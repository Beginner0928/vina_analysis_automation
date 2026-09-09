from __future__ import annotations

import copy
import csv
import json
import sys
import unittest
import uuid
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import batch_execution_v05 as batch_execution_module  # noqa: E402
from batch_execution_v05 import (  # noqa: E402
    CandidateState,
    CandidateExecutionContext,
    ExecutionAdapters,
    REQUIRED_ANALYSIS_FILES,
    V04_FROZEN_ENSEMBLE_RULES,
    build_candidate_identity,
    build_ensemble_spec,
    create_attempt,
    execute_planned_batch,
    execute_batch,
    find_resumable_attempt,
    identity_digest,
    make_attempt_id,
    publish_completion,
    transition_state,
    validate_vina_result,
    validate_completed_attempt,
)


def identity_fixture() -> dict[str, object]:
    return {
        "candidate": {"candidate_id": "X1", "group": "A", "sequence": "ACD"},
        "contracts": {
            "candidate_manifest": {"id": "top40", "sha256": "1" * 64},
            "screening_protocol": {"id": "protocol", "sha256": "2" * 64},
            "receptor_registry": {"id": "registry", "sha256": "3" * 64},
            "ligand_inventory": {"id": "inventory", "sha256": "4" * 64},
        },
        "receptor": {"template_id": "template_a", "receptor_id": "receptor_a"},
        "ligand_preparation": {
            "protocol_id": "prep_v1",
            "sha256": "5" * 64,
        },
        "ligand_generation_provenance": {
            "protocol_id": "generation_v1",
            "sha256": "6" * 64,
        },
        "conformers": [
            {
                "conformer_name": "conf01",
                "source_sdf_sha256": "7" * 64,
                "audit_sha256": "8" * 64,
            }
        ],
        "vina_protocol": {
            "sha256": "9" * 64,
            "expected_model_labels": [1, 2],
        },
    }


class CandidateStateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.attempt = ROOT / "test_output" / f"v05b_state_{uuid.uuid4().hex}"
        self.attempt.mkdir(parents=True)

    def test_valid_state_sequence_is_append_only(self) -> None:
        for state in (
            CandidateState.PENDING,
            CandidateState.PREFLIGHT_PASS,
            CandidateState.RUNNING,
            CandidateState.COMPLETE,
        ):
            transition_state(self.attempt, state, {"source": "test"})
        records = sorted((self.attempt / "states").glob("*.json"))
        self.assertEqual(4, len(records))
        self.assertEqual(
            ["PENDING", "PREFLIGHT_PASS", "RUNNING", "COMPLETE"],
            [path.stem.split("_", 1)[1] for path in records],
        )

    def test_complete_cannot_follow_pending(self) -> None:
        transition_state(self.attempt, CandidateState.PENDING, {})
        with self.assertRaisesRegex(ValueError, "Invalid candidate state transition"):
            transition_state(self.attempt, CandidateState.COMPLETE, {})
        self.assertEqual(1, len(list((self.attempt / "states").glob("*.json"))))

    def test_failure_state_is_terminal(self) -> None:
        transition_state(self.attempt, CandidateState.PENDING, {})
        transition_state(self.attempt, CandidateState.FAILED_PREFLIGHT, {})
        with self.assertRaisesRegex(ValueError, "terminal"):
            transition_state(self.attempt, CandidateState.RUNNING, {})

    def test_complete_remains_terminal_and_cannot_transition_to_failure(self) -> None:
        for state in (
            CandidateState.PENDING,
            CandidateState.PREFLIGHT_PASS,
            CandidateState.RUNNING,
            CandidateState.COMPLETE,
        ):
            transition_state(self.attempt, state, {})
        with self.assertRaisesRegex(ValueError, "terminal"):
            transition_state(self.attempt, CandidateState.FAILED_ANALYSIS, {})


class AttemptIdentityTests(unittest.TestCase):
    def test_attempt_id_is_auditable_and_deterministic_for_supplied_time(self) -> None:
        created = datetime(2026, 9, 8, 12, 34, 56, tzinfo=timezone.utc)
        attempt_id = make_attempt_id(identity_fixture(), created)
        self.assertEqual(
            f"20260908T123456Z_{identity_digest(identity_fixture())[:12]}",
            attempt_id,
        )

    def test_source_sdf_change_changes_identity_digest(self) -> None:
        changed = copy.deepcopy(identity_fixture())
        changed["conformers"][0]["source_sdf_sha256"] = "a" * 64  # type: ignore[index]
        self.assertNotEqual(identity_digest(identity_fixture()), identity_digest(changed))

    def test_inventory_and_ligand_protocol_changes_invalidate_identity(self) -> None:
        for path, value in (
            (("contracts", "ligand_inventory", "sha256"), "a" * 64),
            (("ligand_preparation", "sha256"), "b" * 64),
            (("ligand_generation_provenance", "sha256"), "c" * 64),
        ):
            changed = copy.deepcopy(identity_fixture())
            cursor: dict[str, object] = changed
            for key in path[:-1]:
                cursor = cursor[key]  # type: ignore[assignment]
            cursor[path[-1]] = value
            self.assertNotEqual(
                identity_digest(identity_fixture()),
                identity_digest(changed),
                path,
            )


class ImmutableAttemptTests(unittest.TestCase):
    def setUp(self) -> None:
        self.root = ROOT / "test_output" / f"v05b_attempt_{uuid.uuid4().hex}"
        self.candidate_root = self.root / "candidates" / "X1"
        self.created = datetime(2026, 9, 8, 12, 34, 56, tzinfo=timezone.utc)
        self.identity = identity_fixture()

    def _write_artifacts(self, attempt: Path, analysis_status: str = "PASS") -> dict[str, Path]:
        output = attempt / "workspace" / "outputs" / "X1_conf01_out.pdbqt"
        metrics = attempt / "workspace" / "analysis" / "X1_conf01_pose_metrics.tsv"
        summary = attempt / "ensemble" / "X1" / "candidate_summary.tsv"
        for path, text in (
            (output, "MODEL 1\nREMARK VINA RESULT: -5.0 0 0\nATOM      1  C   LIG A   1       0.0 0.0 0.0  0.0  0.0 C\nENDMDL\n"),
            (metrics, "model\tscore\n1\t-5.0\n"),
        ):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding="utf-8", newline="\n")
        summary.parent.mkdir(parents=True, exist_ok=True)
        with summary.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["candidate", "analysis_status", "warning_count"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerow(
                {
                    "candidate": "X1",
                    "analysis_status": analysis_status,
                    "warning_count": "3",
                }
            )
        return {
            "vina_output.conf01": output,
            "pose_metrics.conf01": metrics,
            "candidate_summary": summary,
        }

    def _complete_attempt(self, analysis_status: str = "PASS") -> tuple[Path, dict[str, Path]]:
        attempt = create_attempt(self.candidate_root, self.identity, self.created)
        transition_state(attempt, CandidateState.PREFLIGHT_PASS, {})
        transition_state(attempt, CandidateState.RUNNING, {})
        artifacts = self._write_artifacts(attempt, analysis_status)
        publish_completion(attempt, self.identity, artifacts, completed_at=self.created)
        return attempt, artifacts

    def test_attempt_collision_refuses_overwrite(self) -> None:
        attempt = create_attempt(self.candidate_root, self.identity, self.created)
        manifest_before = (attempt / "attempt_manifest.json").read_bytes()
        with self.assertRaisesRegex(FileExistsError, "attempt collision"):
            create_attempt(self.candidate_root, self.identity, self.created)
        self.assertEqual(manifest_before, (attempt / "attempt_manifest.json").read_bytes())

    def test_valid_complete_attempt_is_resumable(self) -> None:
        attempt, _ = self._complete_attempt()
        result = validate_completed_attempt(attempt, self.identity)
        self.assertTrue(result.valid, result.message)
        self.assertEqual("VALID_COMPLETE", result.code)
        found = find_resumable_attempt(self.candidate_root, self.identity)
        self.assertTrue(found.valid)
        self.assertEqual(attempt, found.attempt_root)

    def test_incomplete_attempt_is_not_resumable(self) -> None:
        attempt = create_attempt(self.candidate_root, self.identity, self.created)
        transition_state(attempt, CandidateState.PREFLIGHT_PASS, {})
        result = validate_completed_attempt(attempt, self.identity)
        self.assertFalse(result.valid)
        self.assertEqual("COMPLETION_MARKER_MISSING", result.code)
        self.assertFalse(find_resumable_attempt(self.candidate_root, self.identity).valid)

    def test_corrupted_artifact_is_not_resumable(self) -> None:
        attempt, artifacts = self._complete_attempt()
        artifacts["vina_output.conf01"].write_text("corrupt", encoding="utf-8")
        result = validate_completed_attempt(attempt, self.identity)
        self.assertFalse(result.valid)
        self.assertEqual("ARTIFACT_HASH_MISMATCH", result.code)

    def test_malformed_completion_marker_is_not_resumable(self) -> None:
        attempt, _ = self._complete_attempt()
        (attempt / "completion.json").write_text(
            "{malformed\n", encoding="utf-8", newline="\n"
        )
        result = validate_completed_attempt(attempt, self.identity)
        self.assertFalse(result.valid)
        self.assertEqual("COMPLETION_MARKER_INVALID", result.code)
        self.assertFalse(find_resumable_attempt(self.candidate_root, self.identity).valid)

    def test_tampered_state_history_is_not_resumable(self) -> None:
        attempt, _ = self._complete_attempt()
        first_state = sorted((attempt / "states").glob("*.json"))[0]
        first_state.write_text("{}\n", encoding="utf-8", newline="\n")
        result = validate_completed_attempt(attempt, self.identity)
        self.assertFalse(result.valid)
        self.assertEqual("STATE_HISTORY_INVALID", result.code)

    def test_scientifically_weak_pass_summary_can_complete(self) -> None:
        attempt, _ = self._complete_attempt(analysis_status="PASS")
        result = validate_completed_attempt(attempt, self.identity)
        self.assertTrue(result.valid)

    def test_technical_analysis_failure_cannot_complete(self) -> None:
        attempt = create_attempt(self.candidate_root, self.identity, self.created)
        transition_state(attempt, CandidateState.PREFLIGHT_PASS, {})
        transition_state(attempt, CandidateState.RUNNING, {})
        artifacts = self._write_artifacts(attempt, analysis_status="FAIL")
        with self.assertRaisesRegex(ValueError, "candidate_summary analysis_status"):
            publish_completion(attempt, self.identity, artifacts, completed_at=self.created)
        self.assertFalse((attempt / "completion.json").exists())
        self.assertNotEqual(
            CandidateState.COMPLETE.value,
            json.loads(sorted((attempt / "states").glob("*.json"))[-1].read_text(encoding="utf-8"))["state"],
        )

    def test_every_resume_identity_dimension_is_compared(self) -> None:
        attempt, _ = self._complete_attempt()
        mutations = (
            (("candidate", "candidate_id"), "X2"),
            (("candidate", "group"), "B"),
            (("candidate", "sequence"), "AAA"),
            (("contracts", "candidate_manifest", "sha256"), "a" * 64),
            (("contracts", "screening_protocol", "sha256"), "b" * 64),
            (("contracts", "receptor_registry", "sha256"), "c" * 64),
            (("contracts", "ligand_inventory", "id"), "inventory-v2"),
            (("contracts", "ligand_inventory", "sha256"), "d" * 64),
            (("receptor", "template_id"), "template_b"),
            (("receptor", "sha256"), {"receptor_monomer_pdb": "0" * 64}),
            (("ligand_preparation", "sha256"), "e" * 64),
            (("ligand_generation_provenance", "sha256"), "f" * 64),
            (("conformers", 0, "conformer_name"), "conf02"),
            (("conformers", 0, "source_sdf_sha256"), "0" * 64),
            (("vina_protocol", "sha256"), "a" * 64),
            (("vina_protocol", "expected_model_labels"), [1, 2, 3]),
        )
        for path, replacement in mutations:
            expected = copy.deepcopy(self.identity)
            cursor: object = expected
            for key in path[:-1]:
                cursor = cursor[key]  # type: ignore[index]
            cursor[path[-1]] = replacement  # type: ignore[index]
            with self.subTest(path=path):
                result = validate_completed_attempt(attempt, expected)
                self.assertFalse(result.valid)
                self.assertEqual("IDENTITY_MISMATCH", result.code)

    def test_manifest_tampering_and_missing_required_artifact_are_rejected(self) -> None:
        attempt, artifacts = self._complete_attempt()
        (attempt / "attempt_manifest.json").write_text("{}\n", encoding="utf-8")
        self.assertEqual(
            "ATTEMPT_MANIFEST_INVALID",
            validate_completed_attempt(attempt, self.identity).code,
        )

        other_root = ROOT / "test_output" / f"v05b_attempt_missing_{uuid.uuid4().hex}"
        other_candidate = other_root / "candidates" / "X1"
        other = create_attempt(other_candidate, self.identity, self.created)
        transition_state(other, CandidateState.PREFLIGHT_PASS, {})
        transition_state(other, CandidateState.RUNNING, {})
        other_artifacts = self._write_artifacts(other)
        other_artifacts["pose_metrics.conf01"] = (
            other / "workspace" / "analysis" / "missing_pose_metrics.tsv"
        )
        with self.assertRaisesRegex(FileNotFoundError, "Required completion artifact"):
            publish_completion(other, self.identity, other_artifacts, completed_at=self.created)
        self.assertFalse((other / "completion.json").exists())


def planned_candidate(candidate_id: str, *, preflight: str = "PASS") -> dict[str, object]:
    return {
        "candidate_id": candidate_id,
        "group": "A",
        "sequence": "ACD",
        "peptide_length": 3,
        "template_id": "template_a",
        "receptor_id": "receptor_a",
        "comparison_protocol_id": "comparison_a",
        "preflight_status": preflight,
        "conformers": [
            {"conformer_name": "conf01"},
            {"conformer_name": "conf02"},
        ],
        "runtime_spec": {"candidate": candidate_id},
    }


def candidate_identity(candidate_id: str) -> dict[str, object]:
    identity = copy.deepcopy(identity_fixture())
    identity["candidate"]["candidate_id"] = candidate_id  # type: ignore[index]
    identity["conformers"] = [
        {"conformer_name": "conf01", "source_sdf_sha256": "7" * 64, "audit_sha256": "8" * 64},
        {"conformer_name": "conf02", "source_sdf_sha256": "9" * 64, "audit_sha256": "a" * 64},
    ]
    return identity


def preflight_fixture(*entries: dict[str, object], global_failure: bool = False) -> dict[str, object]:
    return {
        "outcome": "DRY_RUN_FAIL_GLOBAL" if global_failure else "DRY_RUN_PASS",
        "selected_candidate_ids": [entry["candidate_id"] for entry in entries],
        "plan": [] if global_failure else list(entries),
        "failures": (
            [{"phase": "global_preflight", "error_code": "GLOBAL_INVALID", "message": "bad contract"}]
            if global_failure
            else []
        ),
    }


class FakeExecution:
    def __init__(
        self,
        *,
        docking_fail: set[str] | None = None,
        analysis_fail: set[str] | None = None,
        bad_models: dict[tuple[str, str], list[int]] | None = None,
    ) -> None:
        self.events: list[str] = []
        self.active_vina = 0
        self.max_active_vina = 0
        self.docking_fail = docking_fail or set()
        self.analysis_fail = analysis_fail or set()
        self.bad_models = bad_models or {}

    @staticmethod
    def _artifact(context, relative: str, text: str) -> Path:
        path = context.attempt_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8", newline="\n")
        return path

    def prepare(self, context, conformer: str) -> dict[str, Path]:
        candidate = context.plan["candidate_id"]
        self.events.append(f"prepare:{candidate}:{conformer}")
        path = self._artifact(context, f"workspace/prepared/{candidate}_{conformer}.json", "{}\n")
        artifacts = {
            f"input_audit.{conformer}": path,
            f"source_sdf.{conformer}": self._artifact(
                context, f"workspace/source_sdf/{candidate}_{conformer}.sdf", "synthetic\n"
            ),
            f"ligand_pdbqt.{conformer}": self._artifact(
                context, f"workspace/ligands/{candidate}_{conformer}.pdbqt", "synthetic\n"
            ),
            f"vina_config.{conformer}": self._artifact(
                context, f"workspace/configs/{candidate}_{conformer}.txt", "synthetic\n"
            ),
        }
        if conformer == "conf01":
            artifacts.update(
                {
                    "protocol": self._artifact(context, "workspace/protocol.json", "{}\n"),
                    "receptor_monomer": self._artifact(
                        context, "workspace/receptors/template_a_TfR1_A.pdb", "synthetic\n"
                    ),
                    "receptor_pdbqt": self._artifact(
                        context, "workspace/receptors/template_a_TfR1_A.pdbqt", "synthetic\n"
                    ),
                    "receptor_dimer": self._artifact(
                        context, "workspace/receptors/template_a_TfR1_AB.pdb", "synthetic\n"
                    ),
                }
            )
        return artifacts

    def audit_input(self, context, conformer: str) -> dict[str, Path]:
        candidate = context.plan["candidate_id"]
        self.events.append(f"input_audit:{candidate}:{conformer}")
        path = self._artifact(context, f"workspace/audit/{candidate}_{conformer}_input.json", "{}\n")
        return {f"pre_docking_gate.{conformer}": path}

    def dock(self, context, conformer: str) -> dict[str, object]:
        candidate = context.plan["candidate_id"]
        self.events.append(f"dock:{candidate}:{conformer}")
        self.active_vina += 1
        self.max_active_vina = max(self.max_active_vina, self.active_vina)
        try:
            if candidate in self.docking_fail:
                raise RuntimeError(f"synthetic docking failure for {candidate}")
            labels = self.bad_models.get(
                (candidate, conformer),
                list(context.identity["vina_protocol"]["expected_model_labels"]),
            )
            output = self._artifact(
                context,
                f"workspace/outputs/{candidate}_{conformer}_out.pdbqt",
                "".join(
                    f"MODEL {label}\nREMARK VINA RESULT: -5.0 0 0\nATOM      1  C   LIG A   1       0.0 0.0 0.0  0.0  0.0 C\nENDMDL\n"
                    for label in labels
                ),
            )
            status = self._artifact(
                context,
                f"workspace/run_status/{candidate}_{conformer}.json",
                json.dumps({"complete": True, "exit_code": 0, "model_labels": labels}) + "\n",
            )
            log = self._artifact(
                context,
                f"workspace/logs/{candidate}_{conformer}.log",
                "EXIT_CODE=0\n",
            )
            return {
                "complete": True,
                "exit_code": 0,
                "model_labels": labels,
                "artifacts": {
                    f"vina_output.{conformer}": output,
                    f"run_status.{conformer}": status,
                    f"vina_log.{conformer}": log,
                },
            }
        finally:
            self.active_vina -= 1

    def metrics(self, context, conformer: str) -> dict[str, Path]:
        candidate = context.plan["candidate_id"]
        self.events.append(f"metrics:{candidate}:{conformer}")
        path = self._artifact(
            context,
            f"workspace/analysis/{candidate}_{conformer}_pose_metrics.tsv",
            "model\tscore\n1\t-5.0\n2\t-4.5\n",
        )
        return {f"pose_metrics.{conformer}": path}

    def audit_output(self, context, conformer: str) -> dict[str, Path]:
        candidate = context.plan["candidate_id"]
        self.events.append(f"output_audit:{candidate}:{conformer}")
        path = self._artifact(context, f"workspace/audit/{candidate}_{conformer}_output.json", "{}\n")
        return {f"post_docking_gate.{conformer}": path}

    def analyze(self, context) -> dict[str, Path]:
        candidate = context.plan["candidate_id"]
        self.events.append(f"analyze:{candidate}")
        if candidate in self.analysis_fail:
            raise RuntimeError(f"synthetic analysis failure for {candidate}")
        artifacts: dict[str, Path] = {}
        analysis_root = context.attempt_root / "ensemble" / candidate
        for name in REQUIRED_ANALYSIS_FILES:
            path = analysis_root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            if name == "candidate_summary.tsv":
                with path.open("x", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=["candidate", "analysis_status", "warning_count"],
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    writer.writerow(
                        {"candidate": candidate, "analysis_status": "PASS", "warning_count": "9"}
                    )
            else:
                path.write_text("synthetic\n", encoding="utf-8", newline="\n")
            artifacts[f"ensemble.{name}"] = path
        artifacts["candidate_summary"] = analysis_root / "candidate_summary.tsv"
        artifacts["ensemble_spec"] = self._artifact(
            context, "workspace/specs/ensemble_spec.json", "{}\n"
        )
        return artifacts

    def adapters(self) -> ExecutionAdapters:
        return ExecutionAdapters(
            prepare_conformer=self.prepare,
            audit_prepared_input=self.audit_input,
            run_vina=self.dock,
            generate_pose_metrics=self.metrics,
            audit_completed_output=self.audit_output,
            analyze_candidate=self.analyze,
        )


class PlannedBatchExecutionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.root = ROOT / "test_output" / f"v05b_batch_{uuid.uuid4().hex}"
        self.created = datetime(2026, 9, 8, 13, 0, 0, tzinfo=timezone.utc)
        self.batch_identity = {
            "candidate_manifest_sha256": "1" * 64,
            "screening_protocol_sha256": "2" * 64,
            "receptor_registry_sha256": "3" * 64,
            "ligand_inventory_sha256": "4" * 64,
        }

    def _run(self, entries, fake: FakeExecution, *, resume: bool = False):
        return execute_planned_batch(
            preflight_fixture(*entries),
            self.root,
            self.batch_identity,
            {entry["candidate_id"]: candidate_identity(str(entry["candidate_id"])) for entry in entries},
            fake.adapters(),
            data_root=ROOT.parent,
            vina_path=ROOT / "test_output" / "never_execute_vina.exe",
            resume=resume,
            started_at=self.created,
        )

    def test_one_candidate_runs_complete_pipeline_and_publishes_last(self) -> None:
        fake = FakeExecution()
        result = self._run([planned_candidate("X1")], fake)
        self.assertEqual("BATCH_COMPLETE", result["outcome"])
        self.assertEqual("COMPLETE", result["candidates"][0]["state"])
        self.assertEqual(
            [
                "prepare:X1:conf01", "input_audit:X1:conf01",
                "prepare:X1:conf02", "input_audit:X1:conf02",
                "dock:X1:conf01", "metrics:X1:conf01", "output_audit:X1:conf01",
                "dock:X1:conf02", "metrics:X1:conf02", "output_audit:X1:conf02",
                "analyze:X1",
            ],
            fake.events,
        )
        attempt = Path(result["candidates"][0]["attempt_root"])
        self.assertTrue((attempt / "completion.json").is_file())
        completion = json.loads((attempt / "completion.json").read_text(encoding="utf-8"))
        self.assertEqual(
            sorted(completion["artifacts"]), completion["required_artifact_names"]
        )
        self.assertEqual(
            ["PENDING", "PREFLIGHT_PASS", "RUNNING", "COMPLETE"],
            [json.loads(path.read_text(encoding="utf-8"))["state"] for path in sorted((attempt / "states").glob("*.json"))],
        )

    def test_multiple_candidates_and_conformers_are_strictly_serial(self) -> None:
        fake = FakeExecution()
        result = self._run([planned_candidate("X1"), planned_candidate("X2")], fake)
        self.assertEqual(["X1", "X2"], [row["candidate_id"] for row in result["candidates"]])
        self.assertEqual(1, fake.max_active_vina)
        self.assertEqual(4, result["vina_docking_process_count"])
        self.assertLess(fake.events.index("analyze:X1"), fake.events.index("prepare:X2:conf01"))

    def test_candidate_preflight_failure_is_isolated(self) -> None:
        fake = FakeExecution()
        entries = [planned_candidate("X1", preflight="FAIL"), planned_candidate("X2")]
        result = self._run(entries, fake)
        self.assertEqual("BATCH_COMPLETE_WITH_FAILURES", result["outcome"])
        self.assertEqual(["FAILED_PREFLIGHT", "COMPLETE"], [row["state"] for row in result["candidates"]])
        self.assertTrue(all("X1" not in event for event in fake.events))

    def test_docking_failure_isolated_and_cannot_publish_completion(self) -> None:
        fake = FakeExecution(docking_fail={"X1"})
        result = self._run([planned_candidate("X1"), planned_candidate("X2")], fake)
        self.assertEqual(["FAILED_DOCKING", "COMPLETE"], [row["state"] for row in result["candidates"]])
        failed = Path(result["candidates"][0]["attempt_root"])
        completed = Path(result["candidates"][1]["attempt_root"])
        self.assertFalse((failed / "completion.json").exists())
        self.assertTrue((completed / "completion.json").exists())

    def test_analysis_failure_isolated(self) -> None:
        fake = FakeExecution(analysis_fail={"X1"})
        result = self._run([planned_candidate("X1"), planned_candidate("X2")], fake)
        self.assertEqual(["FAILED_ANALYSIS", "COMPLETE"], [row["state"] for row in result["candidates"]])

    def test_completion_write_failure_is_candidate_isolated_and_not_resumable(self) -> None:
        original = batch_execution_module._write_json_once

        def fail_first_completion(path, value):
            if path.name == "completion.json" and "X1" in path.parts:
                raise OSError("synthetic completion marker write failure")
            return original(path, value)

        fake = FakeExecution()
        with patch(
            "batch_execution_v05._write_json_once", side_effect=fail_first_completion
        ):
            result = self._run(
                [planned_candidate("X1"), planned_candidate("X2")], fake
            )

        self.assertEqual(
            ["FAILED_ANALYSIS", "COMPLETE"],
            [row["state"] for row in result["candidates"]],
        )
        failed = Path(result["candidates"][0]["attempt_root"])
        self.assertTrue(failed.is_dir())
        self.assertFalse((failed / "completion.json").exists())
        self.assertNotIn(
            "COMPLETE",
            [
                json.loads(path.read_text(encoding="utf-8"))["state"]
                for path in sorted((failed / "states").glob("*.json"))
            ],
        )
        self.assertIn("analyze:X2", fake.events)

        resumed = self._run(
            [planned_candidate("X1"), planned_candidate("X2")],
            FakeExecution(),
            resume=True,
        )
        self.assertEqual("FAILED", resumed["candidates"][0]["action"])
        self.assertEqual(
            "RESUME_NO_VALID_COMPLETE_ATTEMPT",
            resumed["candidates"][0]["error_code"],
        )
        self.assertEqual("SKIPPED_COMPLETE", resumed["candidates"][1]["action"])

    def test_corrupt_completion_after_write_fails_before_complete_state(self) -> None:
        original = batch_execution_module._write_json_once

        def corrupt_first_completion(path, value):
            written = original(path, value)
            if path.name == "completion.json" and "X1" in path.parts:
                path.write_text("{malformed\n", encoding="utf-8", newline="\n")
            return written

        with patch(
            "batch_execution_v05._write_json_once", side_effect=corrupt_first_completion
        ):
            result = self._run(
                [planned_candidate("X1"), planned_candidate("X2")], FakeExecution()
            )

        self.assertEqual(
            ["FAILED_ANALYSIS", "COMPLETE"],
            [row["state"] for row in result["candidates"]],
        )
        failed = Path(result["candidates"][0]["attempt_root"])
        self.assertTrue((failed / "completion.json").is_file())
        self.assertNotIn(
            "COMPLETE",
            [
                json.loads(path.read_text(encoding="utf-8"))["state"]
                for path in sorted((failed / "states").glob("*.json"))
            ],
        )

    def test_model_inventory_audit_failure_is_isolated(self) -> None:
        fake = FakeExecution(bad_models={("X1", "conf01"): [1]})
        result = self._run([planned_candidate("X1"), planned_candidate("X2")], fake)
        self.assertEqual(
            ["FAILED_DOCKING", "COMPLETE"],
            [row["state"] for row in result["candidates"]],
        )
        self.assertIn("MODEL", result["candidates"][0]["message"])

    def test_global_failure_stops_before_output_creation(self) -> None:
        fake = FakeExecution()
        result = execute_planned_batch(
            preflight_fixture(global_failure=True),
            self.root,
            self.batch_identity,
            {},
            fake.adapters(),
            data_root=ROOT.parent,
            vina_path=ROOT / "test_output" / "never_execute_vina.exe",
            started_at=self.created,
        )
        self.assertEqual("BATCH_FAIL_GLOBAL", result["outcome"])
        self.assertFalse(self.root.exists())
        self.assertEqual([], fake.events)

    def test_normal_collision_refuses_existing_batch_without_changes(self) -> None:
        fake = FakeExecution()
        first = self._run([planned_candidate("X1")], fake)
        manifest = self.root / "batch_manifest.json"
        before = manifest.read_bytes()
        second = self._run([planned_candidate("X1")], FakeExecution())
        self.assertEqual("BATCH_FAIL_GLOBAL", second["outcome"])
        self.assertEqual("OUTPUT_PATH_COLLISION", second["failures"][0]["error_code"])
        self.assertEqual(before, manifest.read_bytes())
        self.assertEqual("COMPLETE", first["candidates"][0]["state"])

    def test_resume_skips_only_valid_complete_attempt(self) -> None:
        first_fake = FakeExecution()
        self._run([planned_candidate("X1")], first_fake)
        resume_fake = FakeExecution()
        resumed = self._run([planned_candidate("X1")], resume_fake, resume=True)
        self.assertEqual("SKIPPED_COMPLETE", resumed["candidates"][0]["action"])
        self.assertEqual([], resume_fake.events)

    def test_resume_does_not_reuse_corrupt_attempt_or_create_replacement(self) -> None:
        first = self._run([planned_candidate("X1")], FakeExecution())
        attempt = Path(first["candidates"][0]["attempt_root"])
        output = next((attempt / "workspace" / "outputs").glob("*.pdbqt"))
        output.write_text("corrupt", encoding="utf-8")
        attempt_count = len(list((self.root / "candidates" / "X1" / "attempts").iterdir()))
        resumed = self._run([planned_candidate("X1")], FakeExecution(), resume=True)
        self.assertEqual("FAILED_PREFLIGHT", resumed["candidates"][0]["state"])
        self.assertEqual("RESUME_NO_VALID_COMPLETE_ATTEMPT", resumed["candidates"][0]["error_code"])
        self.assertEqual(attempt_count, len(list((self.root / "candidates" / "X1" / "attempts").iterdir())))

    def test_interrupted_subset_resume_preserves_complete_and_runs_missing_candidate(self) -> None:
        self._run([planned_candidate("X1")], FakeExecution())
        fake = FakeExecution()
        resumed = self._run([planned_candidate("X1"), planned_candidate("X2")], fake, resume=True)
        self.assertEqual(["SKIPPED_COMPLETE", "EXECUTED"], [row["action"] for row in resumed["candidates"]])
        self.assertTrue(all("X1" not in event for event in fake.events))
        self.assertIn("analyze:X2", fake.events)

    def test_later_candidate_failure_cannot_modify_prior_complete_package(self) -> None:
        first = self._run([planned_candidate("X1")], FakeExecution())
        completion = Path(first["candidates"][0]["attempt_root"]) / "completion.json"
        before = completion.read_bytes()
        resumed = self._run(
            [planned_candidate("X1"), planned_candidate("X2")],
            FakeExecution(analysis_fail={"X2"}),
            resume=True,
        )
        self.assertEqual(
            ["SKIPPED_COMPLETE", "FAILED"],
            [row["action"] for row in resumed["candidates"]],
        )
        self.assertEqual(before, completion.read_bytes())

    def test_resume_batch_contract_mismatch_fails_globally(self) -> None:
        self._run([planned_candidate("X1")], FakeExecution())
        changed = dict(self.batch_identity, screening_protocol_sha256="f" * 64)
        result = execute_planned_batch(
            preflight_fixture(planned_candidate("X1")),
            self.root,
            changed,
            {"X1": candidate_identity("X1")},
            FakeExecution().adapters(),
            data_root=ROOT.parent,
            vina_path=ROOT / "test_output" / "never_execute_vina.exe",
            resume=True,
            started_at=self.created,
        )
        self.assertEqual("BATCH_FAIL_GLOBAL", result["outcome"])
        self.assertEqual("BATCH_IDENTITY_MISMATCH", result["failures"][0]["error_code"])


class VinaResultIntegrityTests(unittest.TestCase):
    def test_success_requires_complete_zero_exit_and_exact_model_labels(self) -> None:
        validate_vina_result(
            {"complete": True, "exit_code": 0, "model_labels": [1, 2]},
            [1, 2],
        )

    def test_nonzero_missing_and_unexpected_models_fail(self) -> None:
        invalid = (
            ({"complete": False, "exit_code": 7, "model_labels": []}, "exit"),
            ({"complete": True, "exit_code": 0, "model_labels": []}, "MODEL"),
            ({"complete": True, "exit_code": 0, "model_labels": [1]}, "MODEL"),
            ({"complete": True, "exit_code": 0, "model_labels": [1, 3]}, "MODEL"),
        )
        for result, message in invalid:
            with self.subTest(result=result):
                with self.assertRaisesRegex(RuntimeError, message):
                    validate_vina_result(result, [1, 2])


class ProductionContractConstructionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.root = ROOT / "test_output" / f"v05b_contract_{uuid.uuid4().hex}"
        self.root.mkdir(parents=True)

    def _json(self, name: str, value: dict[str, object]) -> Path:
        path = self.root / name
        path.write_text(json.dumps(value) + "\n", encoding="utf-8", newline="\n")
        return path

    def test_candidate_identity_records_full_lineage_and_unfrozen_generation_dependency(self) -> None:
        manifest = {"manifest_id": "top40"}
        protocol = {
            "protocol_id": "protocol",
            "contracts": {
                "chemistry_contract": {
                    "contract_id": "chemistry",
                    "path": "chemistry.json",
                    "sha256": "c" * 64,
                },
                "analysis_schema": {
                    "path": "analysis.schema.json",
                    "sha256": "d" * 64,
                },
                "dependency_lock": {
                    "path": "requirements-lock.txt",
                    "sha256": "e" * 64,
                },
            },
            "ligand_preparation": {"protocol_id": "prep"},
            "docking_protocol": {"num_modes": 2, "cpu": 8},
            "batch_execution": {"vina_process_concurrency": 1},
        }
        inventory = {"inventory_id": "inventory"}
        manifest_path = self._json("manifest.json", manifest)
        protocol_path = self._json("protocol.json", protocol)
        inventory_path = self._json("inventory.json", inventory)
        registry_path = self._json("registry.json", {"registry_id": "registry"})
        vina = self.root / "vina.exe"
        vina.write_bytes(b"synthetic vina")
        plan = planned_candidate("X1")
        plan["runtime_spec"]["receptor_registry_path"] = "specs/registry.json"
        for index, row in enumerate(plan["conformers"]):
            row.update(
                {
                    "source_sdf_path": f"source/X1_conf0{index + 1}.sdf",
                    "source_sdf_sha256": str(index + 1) * 64,
                    "audit_path": f"audit/X1_conf0{index + 1}.json",
                    "audit_sha256": str(index + 3) * 64,
                    "approval_basis": "formal",
                }
            )
        bundle = {
            "registry_id": "registry",
            "registry_path": registry_path,
            "registry_sha256": "a" * 64,
            "template_id": "template_a",
            "group": "A",
            "source_paths": {"receptor_monomer_pdb": "receptor.pdb"},
            "sha256": {"receptor_monomer_pdb": "b" * 64},
            "target_residues": [150],
            "vina_box": {"center": [0, 0, 0], "size": [10, 10, 10]},
        }
        identity = build_candidate_identity(
            plan,
            manifest=manifest,
            manifest_path=manifest_path,
            protocol=protocol,
            protocol_path=protocol_path,
            inventory=inventory,
            inventory_path=inventory_path,
            receptor_bundle=bundle,
            vina_path=vina,
        )
        self.assertEqual([1, 2], identity["vina_protocol"]["expected_model_labels"])
        self.assertEqual(1, identity["vina_protocol"]["process_concurrency"])
        self.assertEqual(
            "EXTERNAL_READINESS_DEPENDENCY",
            identity["ligand_generation_provenance"]["status"],
        )
        self.assertEqual(2, len(identity["conformers"]))
        self.assertNotIn("A4", json.dumps(identity))
        self.assertNotIn("B3", json.dumps(identity))

    def test_runtime_ensemble_spec_is_generic_and_protocol_driven(self) -> None:
        attempt = self.root / "attempt"
        workspace = attempt / "workspace"
        plan = planned_candidate("X7")
        plan.update(
            {
                "receptor_id": "template_a_TfR1_A",
                "template_id": "template_a",
                "target_residues": [150, 154],
                "vina_box": {"center": [1.0, 2.0, 3.0], "size": [20.0, 21.0, 22.0]},
                "runtime_spec": {
                    "receptor_registry_path": "vina_analysis_automation/specs/registry.json",
                    "receptor_registry_sha256": "a" * 64,
                    "formal_chemistry": {"formal_charge": 1},
                },
            }
        )
        identity = candidate_identity("X7")
        context = CandidateExecutionContext(
            attempt_root=attempt,
            plan=plan,
            identity=identity,
            data_root=ROOT,
            vina_path=self.root / "vina.exe",
        )
        for name in ("template_a_TfR1_A.pdb", "template_a_TfR1_A.pdbqt", "template_a_TfR1_AB.pdb"):
            path = workspace / "receptors" / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("synthetic\n", encoding="utf-8", newline="\n")
        for conformer in ("conf01", "conf02"):
            run_id = f"X7_{conformer}"
            for relative in (
                f"source_sdf/{run_id}.sdf",
                f"ligands/{run_id}.pdbqt",
                f"configs/{run_id}.txt",
                f"logs/{run_id}.log",
                f"outputs/{run_id}_out.pdbqt",
                f"analysis/{run_id}_pose_metrics.tsv",
            ):
                path = workspace / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("synthetic\n", encoding="utf-8", newline="\n")
        protocol = {
            "screening_protocol_status": "formal_standardized",
            "contracts": {
                "chemistry_contract": {
                    "path": "chemistry.json",
                    "sha256": "b" * 64,
                }
            },
            "ligand_preparation": {"protocol_id": "prep"},
            "docking_protocol": {
                "engine": "AutoDock Vina",
                "version": "1.2.7",
                "scoring_function": "vina",
                "exhaustiveness": 32,
                "num_modes": 2,
                "energy_range": 5,
                "seed": 1701,
                "cpu": 8,
            },
        }
        spec = build_ensemble_spec(context, protocol)
        self.assertEqual("X7", spec["candidate"])
        self.assertEqual(4, spec["expected_total_pose_count"])
        self.assertEqual(V04_FROZEN_ENSEMBLE_RULES["contact_cutoffs"], spec["contact_cutoffs"])
        self.assertTrue(all(not Path(row["docking_output_path"]).is_absolute() for row in spec["conformers"]))
        self.assertNotIn("A4", json.dumps(spec))
        self.assertNotIn("B3", json.dumps(spec))

    def test_frozen_v04_analysis_rules_match_both_formal_reference_specs(self) -> None:
        for name in ("A4_formal_ensemble_v04.json", "B3_ensemble_v04.json"):
            value = json.loads((ROOT / "specs" / name).read_text(encoding="utf-8"))
            for field, expected in V04_FROZEN_ENSEMBLE_RULES.items():
                self.assertEqual(expected, value[field], (name, field))

    def test_real_entrypoint_treats_unready_ligand_as_candidate_failure_without_vina(self) -> None:
        output = self.root / "unready_execution"
        fake = FakeExecution()
        with patch(
            "batch_execution_v05.production_execution_adapters",
            return_value=fake.adapters(),
        ):
            result = execute_batch(
                ROOT / "specs" / "top40_candidates_v05.json",
                ROOT / "specs" / "screening_protocol_v05.json",
                ROOT.parent,
                output,
                ROOT.parent / "tools" / "AutoDockVina" / "vina.exe",
                candidate_ids=["C2"],
            )
        self.assertEqual("BATCH_COMPLETE_WITH_FAILURES", result["outcome"])
        self.assertEqual(0, result["vina_docking_process_count"])
        self.assertEqual("FAILED_PREFLIGHT", result["candidates"][0]["state"])
        self.assertEqual([], fake.events)
        self.assertFalse((output / "candidates" / "C2" / "attempts").exists())

    def test_real_entrypoint_executes_ready_plan_only_through_injected_fake(self) -> None:
        output = self.root / "ready_fake_execution"
        fake = FakeExecution()
        with patch(
            "batch_execution_v05.production_execution_adapters",
            return_value=fake.adapters(),
        ):
            result = execute_batch(
                ROOT / "specs" / "top40_candidates_v05.json",
                ROOT / "specs" / "screening_protocol_v05.json",
                ROOT.parent,
                output,
                ROOT.parent / "tools" / "AutoDockVina" / "vina.exe",
                candidate_ids=["A4"],
            )
        self.assertEqual("BATCH_COMPLETE", result["outcome"])
        self.assertEqual(3, result["vina_docking_process_count"])
        self.assertEqual("COMPLETE", result["candidates"][0]["state"])
        self.assertEqual(1, fake.max_active_vina)
        self.assertEqual(
            ["dock:A4:conf01", "dock:A4:conf02", "dock:A4:conf03"],
            [event for event in fake.events if event.startswith("dock:")],
        )


if __name__ == "__main__":
    unittest.main()
