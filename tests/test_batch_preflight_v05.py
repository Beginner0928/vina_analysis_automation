from __future__ import annotations

import copy
import hashlib
import json
import subprocess
import sys
import unittest
import uuid
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT.parent
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from batch_preflight_v05 import (  # noqa: E402
    inspect_vina_identity,
    project_candidate_paths,
    run_dry_run,
)


MANIFEST = ROOT / "specs" / "top40_candidates_v05.json"
PROTOCOL = ROOT / "specs" / "screening_protocol_v05.json"
VINA = DATA_ROOT / "tools" / "AutoDockVina" / "vina.exe"


def unused_output(name: str) -> Path:
    return ROOT / "test_output" / f"v05a_{name}_{uuid.uuid4().hex}"


def write_contract_fixture(
    name: str,
    *,
    manifest_change=None,
    inventory_change=None,
    protocol_change=None,
) -> tuple[Path, Path]:
    folder = unused_output(name)
    folder.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    inventory = json.loads((ROOT / "specs" / "top40_ligand_inventory_v05.json").read_text(encoding="utf-8"))
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if manifest_change is not None:
        manifest_change(manifest)
    manifest_path = folder / "top40_candidates_v05.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n")
    inventory["candidate_manifest_sha256"] = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    if inventory_change is not None:
        inventory_change(inventory)
    inventory_path = folder / "top40_ligand_inventory_v05.json"
    inventory_path.write_text(json.dumps(inventory, indent=2) + "\n", encoding="utf-8", newline="\n")
    protocol["contracts"]["candidate_manifest"]["sha256"] = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    protocol["contracts"]["ligand_inventory"]["sha256"] = hashlib.sha256(inventory_path.read_bytes()).hexdigest()
    if protocol_change is not None:
        protocol_change(protocol)
    protocol_path = folder / "screening_protocol_v05.json"
    protocol_path.write_text(json.dumps(protocol, indent=2) + "\n", encoding="utf-8", newline="\n")
    return manifest_path, protocol_path


class BatchDryRunTests(unittest.TestCase):
    def test_single_ready_candidate_passes_and_creates_no_output(self) -> None:
        output = unused_output("a4")
        with patch.object(subprocess, "run", side_effect=AssertionError("dry-run spawned a process")), patch.object(
            subprocess, "Popen", side_effect=AssertionError("dry-run spawned a process")
        ):
            result = run_dry_run(
                manifest_path=MANIFEST,
                protocol_path=PROTOCOL,
                data_root=DATA_ROOT,
                output_root=output,
                vina_path=VINA,
                candidate_ids=["A4"],
            )
        self.assertEqual("DRY_RUN_PASS", result["outcome"])
        self.assertEqual(["A4"], [row["candidate_id"] for row in result["plan"]])
        self.assertEqual(3, len(result["plan"][0]["conformers"]))
        self.assertEqual("v04_ca42ecda7552d7141bf9", result["plan"][0]["comparison_protocol_id"])
        self.assertEqual("PASS", result["global_checks"]["analysis_schema"]["status"])
        contracts = result["global_checks"]["tracked_contracts"]
        self.assertEqual("PASS", contracts["screening_protocol"]["status"])
        self.assertEqual(
            hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
            contracts["screening_protocol"]["sha256"],
        )
        for name in ("candidate_manifest", "ligand_inventory", "screening_protocol"):
            self.assertRegex(contracts[name]["schema_sha256"], r"^[0-9a-f]{64}$")
            self.assertEqual("PASS", contracts[name]["schema_status"])
        self.assertFalse(output.exists())

    def test_b3_ready_candidate_uses_frozen_group_b_protocol_identity(self) -> None:
        result = run_dry_run(
            MANIFEST,
            PROTOCOL,
            DATA_ROOT,
            unused_output("b3"),
            VINA,
            candidate_ids=["B3"],
        )
        self.assertEqual("DRY_RUN_PASS", result["outcome"])
        self.assertEqual("h26d3_boltz_msa_s2718", result["plan"][0]["template_id"])
        self.assertEqual("v04_b79103f6b9650da06a0a", result["plan"][0]["comparison_protocol_id"])

    def test_ab_and_cd_group_plans_are_canonical_and_report_incomplete_inputs(self) -> None:
        ab = run_dry_run(MANIFEST, PROTOCOL, DATA_ROOT, unused_output("ab"), VINA, groups=["B", "A"])
        cd = run_dry_run(MANIFEST, PROTOCOL, DATA_ROOT, unused_output("cd"), VINA, groups=["D", "C"])
        self.assertEqual("DRY_RUN_FAIL_CANDIDATE", ab["outcome"])
        self.assertEqual("DRY_RUN_FAIL_CANDIDATE", cd["outcome"])
        self.assertEqual([f"A{i}" for i in range(1, 11)] + [f"B{i}" for i in range(1, 11)], ab["selected_candidate_ids"])
        self.assertEqual([f"C{i}" for i in range(1, 11)] + [f"D{i}" for i in range(1, 11)], cd["selected_candidate_ids"])
        self.assertEqual(20, ab["selected_candidate_count"])
        self.assertEqual(20, cd["selected_candidate_count"])
        self.assertEqual(ab["selected_candidate_ids"], [row["candidate_id"] for row in ab["plan"]])
        self.assertEqual(cd["selected_candidate_ids"], [row["candidate_id"] for row in cd["plan"]])
        self.assertEqual(2, ab["planned_ready_candidate_count"])
        self.assertEqual(0, cd["planned_ready_candidate_count"])
        self.assertTrue(any(row["error_code"] == "LIGAND_NOT_PREPARED" for row in ab["failures"]))
        self.assertTrue(any(row["error_code"] == "LIGAND_NOT_APPROVED" for row in cd["failures"]))
        c2 = next(row for row in cd["plan"] if row["candidate_id"] == "C2")
        d8 = next(row for row in cd["plan"] if row["candidate_id"] == "D8")
        self.assertEqual(["PASS", "PASS", "PASS"], [row["artifact_validation_status"] for row in c2["conformers"]])
        self.assertEqual(["PASS", "PASS", "PASS"], [row["artifact_validation_status"] for row in d8["conformers"]])

    def test_projected_conformer_paths_are_deterministic_and_protocol_driven(self) -> None:
        paths = project_candidate_paths("A4", ["shape_a", "shape_b"])
        self.assertEqual("candidates/A4", paths["candidate_root"])
        self.assertEqual(["shape_a", "shape_b"], list(paths["conformers"]))
        self.assertEqual(
            "candidates/A4/outputs/A4_shape_b_out.pdbqt",
            paths["conformers"]["shape_b"]["vina_output"],
        )

    def test_missing_vina_and_wrong_vina_identity_are_global_failures(self) -> None:
        missing = run_dry_run(MANIFEST, PROTOCOL, DATA_ROOT, unused_output("missing_vina"), DATA_ROOT / "missing-vina.exe", candidate_ids=["A4"])
        self.assertEqual("DRY_RUN_FAIL_GLOBAL", missing["outcome"])
        self.assertTrue(any(row["error_code"] == "VINA_EXECUTABLE_MISSING" for row in missing["failures"]))

        fake = unused_output("fake_vina").with_suffix(".exe")
        fake.parent.mkdir(parents=True, exist_ok=True)
        fake.write_bytes(b"not vina")
        wrong = run_dry_run(MANIFEST, PROTOCOL, DATA_ROOT, unused_output("wrong_vina"), fake, candidate_ids=["A4"])
        self.assertEqual("DRY_RUN_FAIL_GLOBAL", wrong["outcome"])
        self.assertTrue(any(row["error_code"] == "VINA_EXECUTABLE_HASH_MISMATCH" for row in wrong["failures"]))

    def test_output_collision_is_global_failure(self) -> None:
        result = run_dry_run(MANIFEST, PROTOCOL, DATA_ROOT, ROOT, VINA, candidate_ids=["A4"])
        self.assertEqual("DRY_RUN_FAIL_GLOBAL", result["outcome"])
        self.assertTrue(any(row["error_code"] == "OUTPUT_PATH_COLLISION" for row in result["failures"]))

    def test_existing_output_root_is_allowed_only_for_explicit_resume_preflight(self) -> None:
        output = unused_output("resume_existing")
        output.mkdir(parents=True)
        default = run_dry_run(
            MANIFEST, PROTOCOL, DATA_ROOT, output, VINA, candidate_ids=["A4"]
        )
        resumed = run_dry_run(
            MANIFEST,
            PROTOCOL,
            DATA_ROOT,
            output,
            VINA,
            candidate_ids=["A4"],
            allow_existing_output_root=True,
        )
        self.assertEqual("DRY_RUN_FAIL_GLOBAL", default["outcome"])
        self.assertEqual("DRY_RUN_PASS", resumed["outcome"])
        self.assertTrue(
            resumed["global_checks"]["output_root"][
                "existing_root_allowed_for_resume"
            ]
        )

    def test_wrong_locked_contract_hash_is_global_failure(self) -> None:
        manifest, altered = write_contract_fixture(
            "bad_registry_contract",
            protocol_change=lambda value: value["contracts"]["receptor_registry"].update(
                {"sha256": "0" * 64}
            ),
        )
        result = run_dry_run(manifest, altered, DATA_ROOT, unused_output("bad_registry"), VINA, candidate_ids=["A4"])
        self.assertEqual("DRY_RUN_FAIL_GLOBAL", result["outcome"])
        self.assertTrue(any(row["error_code"] == "RECEPTOR_REGISTRY_HASH_MISMATCH" for row in result["failures"]))

    def test_missing_ligand_and_ligand_hash_mismatch_are_candidate_failures(self) -> None:
        def missing_change(inventory: dict) -> None:
            slot = next(row for row in inventory["slots"] if row["candidate_id"] == "A4" and row["conformer_name"] == "conf01")
            slot["source_sdf_path"] = "docking_test/missing/A4_conf01.sdf"

        manifest, protocol = write_contract_fixture("missing_ligand", inventory_change=missing_change)
        missing = run_dry_run(manifest, protocol, DATA_ROOT, unused_output("missing_ligand_out"), VINA, candidate_ids=["A4"])
        self.assertEqual("DRY_RUN_FAIL_CANDIDATE", missing["outcome"])
        self.assertTrue(any(row["error_code"] == "LIGAND_MISSING" for row in missing["failures"]))

        def hash_change(inventory: dict) -> None:
            slot = next(row for row in inventory["slots"] if row["candidate_id"] == "A4" and row["conformer_name"] == "conf01")
            slot["source_sdf_sha256"] = "0" * 64

        manifest, protocol = write_contract_fixture("bad_ligand_hash", inventory_change=hash_change)
        mismatch = run_dry_run(manifest, protocol, DATA_ROOT, unused_output("bad_ligand_hash_out"), VINA, candidate_ids=["A4"])
        self.assertEqual("DRY_RUN_FAIL_CANDIDATE", mismatch["outcome"])
        self.assertTrue(any(row["error_code"] == "LIGAND_HASH_MISMATCH" for row in mismatch["failures"]))

    def test_candidate_sequence_chemistry_mismatch_is_candidate_failure(self) -> None:
        def change_sequence(manifest: dict) -> None:
            row = next(item for item in manifest["candidates"] if item["candidate_id"] == "A4")
            row["sequence"] = "SPSNFITMYEW"

        manifest, protocol = write_contract_fixture("chemistry_mismatch", manifest_change=change_sequence)
        result = run_dry_run(manifest, protocol, DATA_ROOT, unused_output("chemistry_mismatch_out"), VINA, candidate_ids=["A4"])
        self.assertEqual("DRY_RUN_FAIL_CANDIDATE", result["outcome"])
        self.assertTrue(any(row["error_code"] == "LIGAND_INVALID" for row in result["failures"]))

    def test_protocol_parameter_or_chemistry_contract_mismatch_is_global_failure(self) -> None:
        manifest, protocol = write_contract_fixture(
            "protocol_mismatch",
            protocol_change=lambda value: value["docking_protocol"].update({"seed": 999}),
        )
        protocol_result = run_dry_run(manifest, protocol, DATA_ROOT, unused_output("protocol_mismatch_out"), VINA, candidate_ids=["A4"])
        self.assertEqual("DRY_RUN_FAIL_GLOBAL", protocol_result["outcome"])

        manifest, protocol = write_contract_fixture(
            "chemistry_contract_mismatch",
            protocol_change=lambda value: value["contracts"]["chemistry_contract"].update(
                {"sha256": "0" * 64}
            ),
        )
        chemistry_result = run_dry_run(manifest, protocol, DATA_ROOT, unused_output("chemistry_contract_mismatch_out"), VINA, candidate_ids=["A4"])
        self.assertEqual("DRY_RUN_FAIL_GLOBAL", chemistry_result["outcome"])
        self.assertTrue(any(row["error_code"] == "CHEMISTRY_CONTRACT_HASH_MISMATCH" for row in chemistry_result["failures"]))

    def test_wrong_declared_vina_version_is_global_failure(self) -> None:
        manifest, protocol = write_contract_fixture(
            "wrong_vina_version",
            protocol_change=lambda value: value["vina_executable"].update(
                {"version": "AutoDock Vina v1.2.6"}
            ),
        )
        result = run_dry_run(manifest, protocol, DATA_ROOT, unused_output("wrong_vina_version_out"), VINA, candidate_ids=["A4"])
        self.assertEqual("DRY_RUN_FAIL_GLOBAL", result["outcome"])
        self.assertEqual("PROTOCOL_OR_SELECTOR_INVALID", result["failures"][0]["error_code"])

    def test_malformed_manifest_is_global_failure(self) -> None:
        malformed = unused_output("malformed") / "top40_candidates_v05.json"
        malformed.parent.mkdir(parents=True, exist_ok=True)
        malformed.write_text("{not-json\n", encoding="utf-8", newline="\n")
        result = run_dry_run(malformed, PROTOCOL, DATA_ROOT, unused_output("malformed_out"), VINA, candidate_ids=["A4"])
        self.assertEqual("DRY_RUN_FAIL_GLOBAL", result["outcome"])
        self.assertEqual("CANDIDATE_MANIFEST_INVALID", result["failures"][0]["error_code"])

    def test_vina_identity_check_is_hash_based_and_never_executes_binary(self) -> None:
        expected = hashlib.sha256(b"binary bytes").hexdigest()
        fake = unused_output("identity").with_suffix(".exe")
        fake.parent.mkdir(parents=True, exist_ok=True)
        fake.write_bytes(b"binary bytes")
        with patch.object(subprocess, "run", side_effect=AssertionError("must not invoke Vina")), patch.object(
            subprocess, "Popen", side_effect=AssertionError("must not invoke Vina")
        ):
            identity = inspect_vina_identity(fake, expected, "AutoDock Vina v1.2.7")
        self.assertEqual(expected, identity["sha256"])
        self.assertEqual("AutoDock Vina v1.2.7", identity["version"])
        self.assertEqual("verified_by_locked_executable_sha256_without_process_execution", identity["version_verification"])


if __name__ == "__main__":
    unittest.main()
