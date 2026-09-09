"""Safe V0.5-B batch execution and immutable resume primitives."""

from __future__ import annotations

import hashlib
import csv
import json
from datetime import datetime, timezone
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Callable


REQUIRED_ANALYSIS_FILES = (
    "pose_table.tsv",
    "pose_target_contacts.tsv",
    "target_contact_frequency.tsv",
    "conformer_summary.tsv",
    "representative_poses.tsv",
    "vina_score_summary.tsv",
    "candidate_summary.tsv",
    "peptide_summary.md",
    "peptide_summary.json",
    "run_manifest.txt",
)


class CandidateState(str, Enum):
    PENDING = "PENDING"
    PREFLIGHT_PASS = "PREFLIGHT_PASS"
    RUNNING = "RUNNING"
    COMPLETE = "COMPLETE"
    FAILED_PREFLIGHT = "FAILED_PREFLIGHT"
    FAILED_DOCKING = "FAILED_DOCKING"
    FAILED_ANALYSIS = "FAILED_ANALYSIS"


TERMINAL_STATES = {
    CandidateState.COMPLETE,
    CandidateState.FAILED_PREFLIGHT,
    CandidateState.FAILED_DOCKING,
    CandidateState.FAILED_ANALYSIS,
}
ALLOWED_TRANSITIONS = {
    None: {CandidateState.PENDING},
    CandidateState.PENDING: {
        CandidateState.PREFLIGHT_PASS,
        CandidateState.FAILED_PREFLIGHT,
    },
    CandidateState.PREFLIGHT_PASS: {
        CandidateState.RUNNING,
        CandidateState.FAILED_PREFLIGHT,
    },
    CandidateState.RUNNING: {
        CandidateState.COMPLETE,
        CandidateState.FAILED_DOCKING,
        CandidateState.FAILED_ANALYSIS,
    },
}


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")


def identity_digest(identity: dict[str, Any]) -> str:
    return hashlib.sha256(canonical_json_bytes(identity)).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def make_attempt_id(identity: dict[str, Any], created_at: datetime) -> str:
    if created_at.tzinfo is None:
        raise ValueError("Attempt creation time must be timezone-aware")
    stamp = created_at.astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{stamp}_{identity_digest(identity)[:12]}"


def _state_records(attempt_root: Path) -> list[Path]:
    state_root = attempt_root / "states"
    return sorted(state_root.glob("[0-9][0-9][0-9][0-9]_*.json")) if state_root.is_dir() else []


def _state_from_record(path: Path) -> CandidateState:
    value = json.loads(path.read_text(encoding="utf-8"))
    return CandidateState(value["state"])


def _prepare_state_transition(
    attempt_root: Path,
    state: CandidateState,
    detail: dict[str, Any],
    recorded_at: datetime | None,
) -> tuple[Path, dict[str, Any]]:
    """Validate and materialize the next state record without publishing it."""
    records = _state_records(attempt_root)
    previous = _state_from_record(records[-1]) if records else None
    if previous in TERMINAL_STATES:
        raise ValueError(f"Candidate state {previous.value} is terminal")
    allowed = ALLOWED_TRANSITIONS.get(previous, set())
    if state not in allowed:
        before = previous.value if previous is not None else "NONE"
        raise ValueError(f"Invalid candidate state transition: {before} -> {state.value}")
    timestamp = recorded_at or datetime.now(timezone.utc)
    if timestamp.tzinfo is None:
        raise ValueError("State timestamp must be timezone-aware")
    output = attempt_root / "states" / f"{len(records) + 1:04d}_{state.value}.json"
    payload = {
        "schema_version": "0.5b",
        "ordinal": len(records) + 1,
        "state": state.value,
        "previous_state": previous.value if previous is not None else None,
        "recorded_at_utc": timestamp.astimezone(timezone.utc).isoformat(),
        "detail": detail,
    }
    return output, payload


def _formatted_json_bytes(value: dict[str, Any]) -> bytes:
    return (json.dumps(value, indent=2, ensure_ascii=False) + "\n").encode("utf-8")


def transition_state(
    attempt_root: Path,
    state: CandidateState,
    detail: dict[str, Any],
    *,
    recorded_at: datetime | None = None,
) -> Path:
    output, payload = _prepare_state_transition(
        attempt_root, state, detail, recorded_at
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return output


def _write_json_once(path: Path, value: dict[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return path


def create_attempt(
    candidate_root: Path,
    identity: dict[str, Any],
    created_at: datetime | None = None,
) -> Path:
    created = created_at or datetime.now(timezone.utc)
    attempt_id = make_attempt_id(identity, created)
    attempts_root = candidate_root / "attempts"
    attempts_root.mkdir(parents=True, exist_ok=True)
    attempt_root = attempts_root / attempt_id
    try:
        attempt_root.mkdir()
    except FileExistsError as exc:
        raise FileExistsError(f"Candidate attempt collision: {attempt_root}") from exc
    manifest = {
        "schema_version": "0.5b",
        "attempt_id": attempt_id,
        "created_at_utc": created.astimezone(timezone.utc).isoformat(),
        "identity_sha256": identity_digest(identity),
        "identity": identity,
        "immutability_policy": "exclusive_create_no_overwrite",
    }
    _write_json_once(attempt_root / "attempt_manifest.json", manifest)
    transition_state(
        attempt_root,
        CandidateState.PENDING,
        {"identity_sha256": manifest["identity_sha256"]},
        recorded_at=created,
    )
    return attempt_root


def _relative_artifact(attempt_root: Path, path: Path) -> str:
    root = attempt_root.resolve()
    resolved = path.resolve()
    if not resolved.is_relative_to(root):
        raise ValueError(f"Completion artifact is outside attempt root: {path}")
    return resolved.relative_to(root).as_posix()


def _validate_candidate_summary(path: Path, candidate_id: str) -> None:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError("candidate_summary must contain exactly one data row")
    row = rows[0]
    if row.get("candidate") != candidate_id:
        raise ValueError("candidate_summary candidate does not match attempt identity")
    if row.get("analysis_status") != "PASS":
        raise ValueError("candidate_summary analysis_status must be PASS")


def _validate_published_completion(
    completion_path: Path, expected_payload: dict[str, Any]
) -> None:
    """Read back the completion marker before the COMPLETE state is appended."""
    try:
        actual_bytes = completion_path.read_bytes()
        actual_payload = json.loads(actual_bytes.decode("utf-8"))
    except Exception as exc:
        raise ValueError(f"Published completion marker is invalid: {exc}") from exc
    if actual_bytes != _formatted_json_bytes(expected_payload) or actual_payload != expected_payload:
        raise ValueError("Published completion marker differs from validated payload")


def publish_completion(
    attempt_root: Path,
    identity: dict[str, Any],
    required_artifacts: dict[str, Path],
    *,
    completed_at: datetime | None = None,
    required_names: set[str] | None = None,
) -> Path:
    if "candidate_summary" not in required_artifacts:
        raise ValueError("candidate_summary is a required completion artifact")
    if required_names is not None:
        missing_names = sorted(required_names - set(required_artifacts))
        if missing_names:
            raise ValueError(
                "Completion artifact inventory is missing required names: "
                + ", ".join(missing_names)
            )
    manifest_path = attempt_root / "attempt_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Attempt manifest is missing: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("identity") != identity or manifest.get("identity_sha256") != identity_digest(identity):
        raise ValueError("Attempt manifest identity differs from completion identity")

    artifact_records: dict[str, dict[str, str]] = {}
    for name, path in sorted(required_artifacts.items()):
        if not path.is_file():
            raise FileNotFoundError(f"Required completion artifact is missing: {name}: {path}")
        relative = _relative_artifact(attempt_root, path)
        artifact_records[name] = {"path": relative, "sha256": file_sha256(path)}
    _validate_candidate_summary(
        required_artifacts["candidate_summary"],
        str(identity["candidate"]["candidate_id"]),
    )

    finished = completed_at or datetime.now(timezone.utc)
    complete_state_path, complete_state_payload = _prepare_state_transition(
        attempt_root,
        CandidateState.COMPLETE,
        {"validated_artifact_count": len(artifact_records)},
        finished,
    )
    complete_state_sha256 = hashlib.sha256(
        _formatted_json_bytes(complete_state_payload)
    ).hexdigest()
    state_history = [
        {
            "path": _relative_artifact(attempt_root, path),
            "sha256": file_sha256(path),
        }
        for path in _state_records(attempt_root)
    ]
    state_history.append(
        {
            "path": _relative_artifact(attempt_root, complete_state_path),
            "sha256": complete_state_sha256,
        }
    )
    payload = {
        "schema_version": "0.5b",
        "status": CandidateState.COMPLETE.value,
        "attempt_id": manifest["attempt_id"],
        "completed_at_utc": finished.astimezone(timezone.utc).isoformat(),
        "identity_sha256": identity_digest(identity),
        "attempt_manifest_sha256": file_sha256(manifest_path),
        "complete_state": {
            "path": _relative_artifact(attempt_root, complete_state_path),
            "sha256": complete_state_sha256,
        },
        "state_history": state_history,
        "artifacts": artifact_records,
        "required_artifact_names": sorted(required_names or required_artifacts),
        "integrity_status": "PASS",
    }
    completion_path = _write_json_once(attempt_root / "completion.json", payload)
    _validate_published_completion(completion_path, payload)
    complete_state = transition_state(
        attempt_root,
        CandidateState.COMPLETE,
        {"validated_artifact_count": len(artifact_records)},
        recorded_at=finished,
    )
    if complete_state != complete_state_path or file_sha256(complete_state) != complete_state_sha256:
        raise RuntimeError("Published COMPLETE state differs from the completion contract")
    return completion_path


@dataclass(frozen=True)
class AttemptValidation:
    valid: bool
    code: str
    message: str
    attempt_root: Path | None = None


def _invalid(code: str, message: str, attempt_root: Path) -> AttemptValidation:
    return AttemptValidation(False, code, message, attempt_root)


def validate_completed_attempt(
    attempt_root: Path, expected_identity: dict[str, Any]
) -> AttemptValidation:
    manifest_path = attempt_root / "attempt_manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if not isinstance(manifest, dict) or not {
            "attempt_id",
            "identity",
            "identity_sha256",
        }.issubset(manifest):
            raise ValueError("required fields are missing")
        if manifest["attempt_id"] != attempt_root.name:
            raise ValueError("attempt_id does not match attempt directory name")
        if manifest["identity_sha256"] != identity_digest(manifest["identity"]):
            raise ValueError("stored identity digest is invalid")
    except Exception as exc:
        return _invalid("ATTEMPT_MANIFEST_INVALID", str(exc), attempt_root)
    if manifest["identity"] != expected_identity or manifest["identity_sha256"] != identity_digest(expected_identity):
        return _invalid(
            "IDENTITY_MISMATCH",
            "Attempt identity does not match current candidate provenance",
            attempt_root,
        )

    completion_path = attempt_root / "completion.json"
    if not completion_path.is_file():
        return _invalid("COMPLETION_MARKER_MISSING", "completion.json is missing", attempt_root)
    try:
        completion = json.loads(completion_path.read_text(encoding="utf-8"))
    except Exception as exc:
        return _invalid("COMPLETION_MARKER_INVALID", str(exc), attempt_root)
    if (
        completion.get("status") != CandidateState.COMPLETE.value
        or completion.get("integrity_status") != "PASS"
        or completion.get("attempt_id") != manifest["attempt_id"]
        or completion.get("identity_sha256") != identity_digest(expected_identity)
    ):
        return _invalid("COMPLETION_MARKER_INVALID", "Completion status or identity is invalid", attempt_root)
    if completion.get("attempt_manifest_sha256") != file_sha256(manifest_path):
        return _invalid("ATTEMPT_MANIFEST_INVALID", "Attempt manifest hash changed", attempt_root)

    complete_state = completion.get("complete_state")
    if not isinstance(complete_state, dict):
        return _invalid("COMPLETE_STATE_INVALID", "Complete state record is missing", attempt_root)
    state_path = (attempt_root / str(complete_state.get("path", ""))).resolve()
    try:
        complete_state_valid = (
            state_path.is_relative_to(attempt_root.resolve())
            and state_path.is_file()
            and file_sha256(state_path) == complete_state.get("sha256")
            and _state_from_record(state_path) is CandidateState.COMPLETE
        )
    except Exception:
        complete_state_valid = False
    if not complete_state_valid:
        return _invalid("COMPLETE_STATE_INVALID", "Complete state record is invalid", attempt_root)
    history = completion.get("state_history")
    records = _state_records(attempt_root)
    if not isinstance(history, list) or len(history) != len(records):
        return _invalid("STATE_HISTORY_INVALID", "State history inventory is incomplete", attempt_root)
    previous: CandidateState | None = None
    for ordinal, (record, path) in enumerate(zip(history, records), start=1):
        if not isinstance(record, dict):
            return _invalid("STATE_HISTORY_INVALID", "State history record is invalid", attempt_root)
        if (
            record.get("path") != _relative_artifact(attempt_root, path)
            or record.get("sha256") != file_sha256(path)
        ):
            return _invalid("STATE_HISTORY_INVALID", "State history hash changed", attempt_root)
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            state = CandidateState(payload["state"])
            if payload.get("ordinal") != ordinal:
                raise ValueError("state ordinal mismatch")
            expected_previous = previous.value if previous is not None else None
            if payload.get("previous_state") != expected_previous:
                raise ValueError("state previous_state mismatch")
            if state not in ALLOWED_TRANSITIONS.get(previous, set()):
                raise ValueError("invalid state transition")
            previous = state
        except Exception as exc:
            return _invalid("STATE_HISTORY_INVALID", str(exc), attempt_root)
    if previous is not CandidateState.COMPLETE:
        return _invalid("STATE_HISTORY_INVALID", "Latest state is not COMPLETE", attempt_root)

    artifacts = completion.get("artifacts")
    if not isinstance(artifacts, dict) or "candidate_summary" not in artifacts:
        return _invalid("ARTIFACT_INVENTORY_INVALID", "Artifact inventory is incomplete", attempt_root)
    for name, record in artifacts.items():
        if not isinstance(record, dict):
            return _invalid("ARTIFACT_INVENTORY_INVALID", f"Invalid artifact record: {name}", attempt_root)
        path = (attempt_root / str(record.get("path", ""))).resolve()
        if not path.is_relative_to(attempt_root.resolve()) or not path.is_file():
            return _invalid("ARTIFACT_MISSING", f"Required artifact is missing: {name}", attempt_root)
        if file_sha256(path) != record.get("sha256"):
            return _invalid("ARTIFACT_HASH_MISMATCH", f"Artifact hash changed: {name}", attempt_root)
    required_names = completion.get("required_artifact_names")
    if (
        not isinstance(required_names, list)
        or any(not isinstance(name, str) for name in required_names)
        or not set(required_names).issubset(artifacts)
    ):
        return _invalid(
            "ARTIFACT_INVENTORY_INVALID",
            "Required artifact-name contract is invalid",
            attempt_root,
        )
    if expected_identity.get("schema_version") == "0.5b":
        expected_names = required_execution_artifact_names(expected_identity)
        if set(required_names) != expected_names:
            return _invalid(
                "ARTIFACT_INVENTORY_INVALID",
                "Completion package contract differs from the V0.5-B required inventory",
                attempt_root,
            )
    try:
        summary_record = artifacts["candidate_summary"]
        _validate_candidate_summary(
            attempt_root / summary_record["path"],
            str(expected_identity["candidate"]["candidate_id"]),
        )
    except Exception as exc:
        return _invalid("CANDIDATE_SUMMARY_INVALID", str(exc), attempt_root)
    return AttemptValidation(True, "VALID_COMPLETE", "Completed attempt is valid", attempt_root)


def find_resumable_attempt(
    candidate_root: Path, expected_identity: dict[str, Any]
) -> AttemptValidation:
    attempts_root = candidate_root / "attempts"
    attempts = sorted(
        (path for path in attempts_root.iterdir() if path.is_dir()), reverse=True
    ) if attempts_root.is_dir() else []
    if not attempts:
        return AttemptValidation(False, "NO_ATTEMPTS", "Candidate has no existing attempts")
    failures: list[str] = []
    for attempt in attempts:
        result = validate_completed_attempt(attempt, expected_identity)
        if result.valid:
            return result
        failures.append(f"{attempt.name}:{result.code}")
    return AttemptValidation(
        False,
        "NO_VALID_COMPLETE_ATTEMPT",
        "; ".join(failures),
    )


@dataclass(frozen=True)
class CandidateExecutionContext:
    attempt_root: Path
    plan: dict[str, Any]
    identity: dict[str, Any]
    data_root: Path
    vina_path: Path


ArtifactFunction = Callable[[CandidateExecutionContext, str], dict[str, Path]]
DockingFunction = Callable[[CandidateExecutionContext, str], dict[str, Any]]
AnalysisFunction = Callable[[CandidateExecutionContext], dict[str, Path]]


@dataclass(frozen=True)
class ExecutionAdapters:
    prepare_conformer: ArtifactFunction
    audit_prepared_input: ArtifactFunction
    run_vina: DockingFunction
    generate_pose_metrics: ArtifactFunction
    audit_completed_output: ArtifactFunction
    analyze_candidate: AnalysisFunction


def validate_vina_result(result: dict[str, Any], expected_model_labels: list[int]) -> None:
    if result.get("complete") is not True or result.get("exit_code") != 0:
        raise RuntimeError(
            f"Vina exit was not successful: complete={result.get('complete')!r}, "
            f"exit_code={result.get('exit_code')!r}"
        )
    labels = result.get("model_labels")
    if labels != expected_model_labels:
        raise RuntimeError(
            f"Vina MODEL labels differ from the frozen protocol: "
            f"expected={expected_model_labels!r}, actual={labels!r}"
        )


def _batch_failure(code: str, message: str) -> dict[str, Any]:
    return {
        "outcome": "BATCH_FAIL_GLOBAL",
        "failures": [
            {
                "phase": "global_preflight",
                "error_code": code,
                "message": message,
            }
        ],
        "candidates": [],
        "vina_docking_process_count": 0,
    }


def _record_candidate_failure(
    attempt_root: Path,
    state: CandidateState,
    candidate_id: str,
    error_code: str,
    exc: Exception,
) -> dict[str, Any]:
    transition_state(
        attempt_root,
        state,
        {"error_code": error_code, "message": str(exc)},
    )
    return {
        "candidate_id": candidate_id,
        "state": state.value,
        "action": "FAILED",
        "error_code": error_code,
        "message": str(exc),
        "attempt_root": str(attempt_root),
    }


def _write_batch_run(output_root: Path, payload: dict[str, Any]) -> Path:
    run_root = output_root / "batch_runs"
    run_root.mkdir(parents=True, exist_ok=True)
    ordinal = 1
    while True:
        path = run_root / f"{ordinal:04d}.json"
        try:
            return _write_json_once(path, payload)
        except FileExistsError:
            ordinal += 1


def _analysis_artifacts(artifacts: dict[str, Path]) -> dict[str, Path]:
    by_basename = {path.name: path for path in artifacts.values()}
    missing = [name for name in REQUIRED_ANALYSIS_FILES if name not in by_basename]
    if missing:
        raise FileNotFoundError(
            "Required candidate analysis artifacts are missing: " + ", ".join(missing)
        )
    result = {f"ensemble.{name}": by_basename[name] for name in REQUIRED_ANALYSIS_FILES}
    result["candidate_summary"] = by_basename["candidate_summary.tsv"]
    return result


def required_execution_artifact_names(identity: dict[str, Any]) -> set[str]:
    names = {
        "runtime_spec",
        "ensemble_spec",
        "protocol",
        "receptor_monomer",
        "receptor_pdbqt",
        "receptor_dimer",
        "candidate_summary",
        *(f"ensemble.{name}" for name in REQUIRED_ANALYSIS_FILES),
    }
    for row in identity["conformers"]:
        conformer = str(row["conformer_name"])
        for kind in (
            "input_audit",
            "pre_docking_gate",
            "source_sdf",
            "ligand_pdbqt",
            "vina_config",
            "vina_output",
            "vina_log",
            "run_status",
            "pose_metrics",
            "post_docking_gate",
        ):
            names.add(f"{kind}.{conformer}")
    return names


def execute_planned_batch(
    preflight: dict[str, Any],
    output_root: Path,
    batch_identity: dict[str, Any],
    identities_by_candidate: dict[str, dict[str, Any]],
    adapters: ExecutionAdapters,
    *,
    data_root: Path,
    vina_path: Path,
    resume: bool = False,
    started_at: datetime | None = None,
) -> dict[str, Any]:
    """Execute a validated plan sequentially with immutable attempt semantics."""
    started = started_at or datetime.now(timezone.utc)
    if preflight.get("outcome") == "DRY_RUN_FAIL_GLOBAL":
        return {
            "outcome": "BATCH_FAIL_GLOBAL",
            "failures": list(preflight.get("failures", [])),
            "candidates": [],
            "vina_docking_process_count": 0,
        }

    manifest_path = output_root / "batch_manifest.json"
    if not resume:
        if output_root.exists():
            return _batch_failure(
                "OUTPUT_PATH_COLLISION",
                f"Batch output root already exists: {output_root}",
            )
        output_root.mkdir(parents=True)
        _write_json_once(
            manifest_path,
            {
                "schema_version": "0.5b",
                "created_at_utc": started.astimezone(timezone.utc).isoformat(),
                "batch_identity_sha256": identity_digest(batch_identity),
                "batch_identity": batch_identity,
                "immutability_policy": "exclusive_create_no_overwrite",
            },
        )
    else:
        if not manifest_path.is_file():
            return _batch_failure(
                "RESUME_BATCH_MANIFEST_MISSING",
                f"Resume requires an existing batch manifest: {manifest_path}",
            )
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception as exc:
            return _batch_failure("RESUME_BATCH_MANIFEST_INVALID", str(exc))
        if (
            manifest.get("batch_identity") != batch_identity
            or manifest.get("batch_identity_sha256") != identity_digest(batch_identity)
        ):
            return _batch_failure(
                "BATCH_IDENTITY_MISMATCH",
                "Existing batch contract identity differs from the requested identity",
            )

    results: list[dict[str, Any]] = []
    vina_count = 0
    for plan_value in preflight.get("plan", []):
        plan = dict(plan_value)
        candidate_id = str(plan["candidate_id"])
        if plan.get("preflight_status") != "PASS":
            results.append(
                {
                    "candidate_id": candidate_id,
                    "state": CandidateState.FAILED_PREFLIGHT.value,
                    "action": "FAILED",
                    "error_code": "CANDIDATE_PREFLIGHT_FAILED",
                    "message": "Candidate did not pass V0.5-A preflight",
                    "attempt_root": None,
                }
            )
            continue
        identity = identities_by_candidate.get(candidate_id)
        if identity is None:
            results.append(
                {
                    "candidate_id": candidate_id,
                    "state": CandidateState.FAILED_PREFLIGHT.value,
                    "action": "FAILED",
                    "error_code": "CANDIDATE_IDENTITY_MISSING",
                    "message": "Candidate execution identity was not supplied",
                    "attempt_root": None,
                }
            )
            continue

        candidate_root = output_root / "candidates" / candidate_id
        if resume and candidate_root.exists():
            existing = find_resumable_attempt(candidate_root, identity)
            if existing.valid:
                results.append(
                    {
                        "candidate_id": candidate_id,
                        "state": CandidateState.COMPLETE.value,
                        "action": "SKIPPED_COMPLETE",
                        "error_code": None,
                        "message": existing.message,
                        "attempt_root": str(existing.attempt_root),
                    }
                )
            else:
                results.append(
                    {
                        "candidate_id": candidate_id,
                        "state": CandidateState.FAILED_PREFLIGHT.value,
                        "action": "FAILED",
                        "error_code": "RESUME_NO_VALID_COMPLETE_ATTEMPT",
                        "message": existing.message,
                        "attempt_root": None,
                    }
                )
            continue

        try:
            attempt_root = create_attempt(candidate_root, identity, started)
        except Exception as exc:
            results.append(
                {
                    "candidate_id": candidate_id,
                    "state": CandidateState.FAILED_PREFLIGHT.value,
                    "action": "FAILED",
                    "error_code": "ATTEMPT_CREATE_FAILED",
                    "message": str(exc),
                    "attempt_root": None,
                }
            )
            continue
        context = CandidateExecutionContext(
            attempt_root=attempt_root,
            plan=plan,
            identity=identity,
            data_root=data_root,
            vina_path=vina_path,
        )
        artifacts: dict[str, Path] = {}
        try:
            runtime_spec_path = _write_json_once(
                attempt_root / "workspace" / "specs" / "candidate_runtime_spec.json",
                dict(plan.get("runtime_spec", {})),
            )
            artifacts["runtime_spec"] = runtime_spec_path
            conformers = [
                str(value["conformer_name"]) for value in plan.get("conformers", [])
            ]
            if not conformers:
                raise ValueError("Candidate plan contains no conformers")
            for conformer in conformers:
                artifacts.update(adapters.prepare_conformer(context, conformer))
                artifacts.update(adapters.audit_prepared_input(context, conformer))
            transition_state(
                attempt_root,
                CandidateState.PREFLIGHT_PASS,
                {"validated_conformers": conformers},
            )
        except Exception as exc:
            results.append(
                _record_candidate_failure(
                    attempt_root,
                    CandidateState.FAILED_PREFLIGHT,
                    candidate_id,
                    "PREPARATION_OR_INPUT_AUDIT_FAILED",
                    exc,
                )
            )
            continue

        transition_state(
            attempt_root,
            CandidateState.RUNNING,
            {"conformer_count": len(conformers), "execution_mode": "strictly_serial"},
        )
        docking_failed = False
        for conformer in conformers:
            try:
                vina_count += 1
                docking = adapters.run_vina(context, conformer)
                validate_vina_result(
                    docking,
                    list(identity["vina_protocol"]["expected_model_labels"]),
                )
                artifacts.update(docking.get("artifacts", {}))
                artifacts.update(adapters.generate_pose_metrics(context, conformer))
                artifacts.update(adapters.audit_completed_output(context, conformer))
            except Exception as exc:
                results.append(
                    _record_candidate_failure(
                        attempt_root,
                        CandidateState.FAILED_DOCKING,
                        candidate_id,
                        "DOCKING_OR_OUTPUT_AUDIT_FAILED",
                        exc,
                    )
                )
                docking_failed = True
                break
        if docking_failed:
            continue

        try:
            analysis = adapters.analyze_candidate(context)
            artifacts.update(analysis)
            artifacts.update(_analysis_artifacts(analysis))
            publish_completion(
                attempt_root,
                identity,
                artifacts,
                required_names=required_execution_artifact_names(identity),
            )
            results.append(
                {
                    "candidate_id": candidate_id,
                    "state": CandidateState.COMPLETE.value,
                    "action": "EXECUTED",
                    "error_code": None,
                    "message": "Candidate execution and analysis completed",
                    "attempt_root": str(attempt_root),
                }
            )
        except Exception as exc:
            results.append(
                _record_candidate_failure(
                    attempt_root,
                    CandidateState.FAILED_ANALYSIS,
                    candidate_id,
                    "ANALYSIS_OR_COMPLETION_FAILED",
                    exc,
                )
            )

    has_failures = any(row["state"].startswith("FAILED_") for row in results)
    outcome = "BATCH_COMPLETE_WITH_FAILURES" if has_failures else "BATCH_COMPLETE"
    payload = {
        "schema_version": "0.5b",
        "started_at_utc": started.astimezone(timezone.utc).isoformat(),
        "resume": resume,
        "outcome": outcome,
        "selected_candidate_ids": list(preflight.get("selected_candidate_ids", [])),
        "vina_docking_process_count": vina_count,
        "candidates": results,
        "failures": [],
    }
    _write_batch_run(output_root, payload)
    return payload


def _load_json_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _resolve_contract_path(
    contract: dict[str, Any], *, protocol_path: Path, data_root: Path
) -> Path:
    bases = {
        "protocol_directory": protocol_path.parent,
        "data_root": data_root,
        "project_root": Path(__file__).resolve().parents[1],
    }
    path_base = contract.get("path_base")
    if path_base not in bases:
        raise ValueError(f"Unsupported contract path_base: {path_base!r}")
    raw = Path(str(contract["path"]))
    if raw.is_absolute():
        raise ValueError("Contract paths must be relative")
    resolved = (bases[path_base] / raw).resolve()
    if not resolved.is_relative_to(bases[path_base].resolve()):
        raise ValueError(f"Contract path escapes {path_base}: {raw}")
    return resolved


def build_candidate_identity(
    plan: dict[str, Any],
    *,
    manifest: dict[str, Any],
    manifest_path: Path,
    protocol: dict[str, Any],
    protocol_path: Path,
    inventory: dict[str, Any],
    inventory_path: Path,
    receptor_bundle: dict[str, Any],
    vina_path: Path,
) -> dict[str, Any]:
    conformers = []
    generation_records = []
    for row in plan["conformers"]:
        conformers.append(
            {
                "conformer_name": row["conformer_name"],
                "source_sdf_path": row["source_sdf_path"],
                "source_sdf_sha256": row["source_sdf_sha256"],
                "audit_path": row["audit_path"],
                "audit_sha256": row["audit_sha256"],
                "approval_basis": row.get("approval_basis"),
            }
        )
        if row.get("generation_provenance") is not None:
            generation_records.append(row["generation_provenance"])
    if generation_records and len(generation_records) != len(conformers):
        raise ValueError(
            "Ligand generation provenance must be present for every conformer or none"
        )
    generation = (
        {
            "status": "PRESENT",
            "records": generation_records,
        }
        if generation_records
        else {
            "status": "EXTERNAL_READINESS_DEPENDENCY",
            "dependency": "UNIFIED_TOP40_LIGAND_GENERATION_LINEAGE_NOT_FROZEN",
        }
    )
    protocol_hash = file_sha256(protocol_path)
    vina_hash = file_sha256(vina_path)
    docking_protocol = protocol["docking_protocol"]
    expected_labels = list(range(1, int(docking_protocol["num_modes"]) + 1))
    return {
        "schema_version": "0.5b",
        "candidate": {
            "candidate_id": plan["candidate_id"],
            "group": plan["group"],
            "sequence": plan["sequence"],
            "peptide_length": plan["peptide_length"],
        },
        "contracts": {
            "candidate_manifest": {
                "id": manifest["manifest_id"],
                "path": manifest_path.name,
                "sha256": file_sha256(manifest_path),
            },
            "screening_protocol": {
                "id": protocol["protocol_id"],
                "path": protocol_path.name,
                "sha256": protocol_hash,
            },
            "ligand_inventory": {
                "id": inventory["inventory_id"],
                "path": inventory_path.name,
                "sha256": file_sha256(inventory_path),
            },
            "receptor_registry": {
                "id": receptor_bundle["registry_id"],
                "path": plan["runtime_spec"]["receptor_registry_path"],
                "sha256": receptor_bundle["registry_sha256"],
            },
            "chemistry_contract": {
                "id": protocol["contracts"]["chemistry_contract"]["contract_id"],
                "path": protocol["contracts"]["chemistry_contract"]["path"],
                "sha256": protocol["contracts"]["chemistry_contract"]["sha256"],
            },
            "analysis_schema": {
                "path": protocol["contracts"]["analysis_schema"]["path"],
                "sha256": protocol["contracts"]["analysis_schema"]["sha256"],
            },
            "dependency_lock": {
                "path": protocol["contracts"]["dependency_lock"]["path"],
                "sha256": protocol["contracts"]["dependency_lock"]["sha256"],
            },
        },
        "receptor": {
            "template_id": receptor_bundle["template_id"],
            "receptor_id": plan["receptor_id"],
            "group": receptor_bundle["group"],
            "source_paths": receptor_bundle["source_paths"],
            "sha256": receptor_bundle["sha256"],
            "target_residues": receptor_bundle["target_residues"],
            "vina_box": receptor_bundle["vina_box"],
        },
        "comparison_protocol_id": plan["comparison_protocol_id"],
        "ligand_preparation": {
            **protocol["ligand_preparation"],
            "contract_sha256": hashlib.sha256(
                canonical_json_bytes(protocol["ligand_preparation"])
            ).hexdigest(),
        },
        "ligand_generation_provenance": generation,
        "conformers": conformers,
        "vina_protocol": {
            **docking_protocol,
            "protocol_sha256": hashlib.sha256(
                canonical_json_bytes(docking_protocol)
            ).hexdigest(),
            "vina_executable_name": vina_path.name,
            "vina_executable_sha256": vina_hash,
            "expected_model_labels": expected_labels,
            "process_concurrency": protocol["batch_execution"][
                "vina_process_concurrency"
            ],
        },
    }


V04_FROZEN_ENSEMBLE_RULES = {
    "contact_cutoffs": {"direct_A": 4.0, "near_A": 5.0},
    "representative_selection_rules": {
        "max_receptor_clash_pairs": 0,
        "max_chain_B_clash_pairs": 0,
        "min_direct_target_contacts": 2,
        "required_any_direct_residues": [],
        "max_per_conformer": 3,
        "diversity_key": "exact_direct_contact_pattern",
    },
    "summary_rules": {"recurrent_contact_frequency_min": 0.25},
}


def _relative_input(path: Path, data_root: Path, label: str) -> str:
    resolved = path.resolve()
    root = data_root.resolve()
    if not resolved.is_relative_to(root):
        raise ValueError(f"{label} must be inside --data-root for provenance: {resolved}")
    return resolved.relative_to(root).as_posix()


def _candidate_workspace(context: CandidateExecutionContext) -> Path:
    return context.attempt_root / "workspace"


def _candidate_spec_path(context: CandidateExecutionContext) -> Path:
    return _candidate_workspace(context) / "specs" / "candidate_runtime_spec.json"


def build_ensemble_spec(
    context: CandidateExecutionContext,
    protocol: dict[str, Any],
) -> dict[str, Any]:
    workspace = _candidate_workspace(context)
    runtime = context.plan["runtime_spec"]
    candidate = str(context.plan["candidate_id"])
    receptor_id = str(context.plan["receptor_id"])
    template_id = str(context.plan["template_id"])
    model_count = int(protocol["docking_protocol"]["num_modes"])
    chemistry_contract = protocol["contracts"]["chemistry_contract"]
    conformers = []
    conformer_hashes: dict[str, dict[str, str]] = {}
    for row in context.plan["conformers"]:
        name = str(row["conformer_name"])
        run_id = f"{candidate}_{name}"
        paths = {
            "starting_sdf_path": workspace / "source_sdf" / f"{run_id}.sdf",
            "ligand_pdbqt_path": workspace / "ligands" / f"{run_id}.pdbqt",
            "vina_config_path": workspace / "configs" / f"{run_id}.txt",
            "vina_log_path": workspace / "logs" / f"{run_id}.log",
            "docking_output_path": workspace / "outputs" / f"{run_id}_out.pdbqt",
            "pose_metrics_path": workspace / "analysis" / f"{run_id}_pose_metrics.tsv",
        }
        conformers.append(
            {
                "conformer_name": name,
                **{
                    field: _relative_input(path, context.data_root, f"{name}.{field}")
                    for field, path in paths.items()
                },
                "expected_model_count": model_count,
            }
        )
        conformer_hashes[name] = {
            field: file_sha256(path) for field, path in paths.items()
        }
    receptor_paths = {
        "receptor_monomer_path": workspace / "receptors" / f"{receptor_id}.pdb",
        "receptor_dimer_path": workspace / "receptors" / f"{template_id}_TfR1_AB.pdb",
        "docking_receptor_pdbqt_path": workspace / "receptors" / f"{receptor_id}.pdbqt",
    }
    docking_protocol = {
        "protocol_id": context.plan["comparison_protocol_id"],
        **protocol["docking_protocol"],
        "score_units": "kcal/mol",
        "score_direction": "lower_is_better",
        "box_center_A": context.plan["vina_box"]["center"],
        "box_size_A": context.plan["vina_box"]["size"],
    }
    return {
        "schema_version": "0.4",
        "analysis_id": f"{candidate}_ensemble_v05b",
        "candidate": candidate,
        "sequence": context.plan["sequence"],
        "peptide_length": context.plan["peptide_length"],
        "group": context.plan["group"],
        "template_id": template_id,
        "receptor_id": receptor_id,
        "receptor_registry_path": runtime["receptor_registry_path"],
        "receptor_registry_sha256": runtime["receptor_registry_sha256"],
        "chemistry_contract_path": chemistry_contract["path"],
        "chemistry_contract_sha256": chemistry_contract["sha256"],
        **{
            field: _relative_input(path, context.data_root, field)
            for field, path in receptor_paths.items()
        },
        "receptor_chain": "A",
        "target_residues": context.plan["target_residues"],
        "screening_protocol_status": protocol["screening_protocol_status"],
        "comparison_protocol_id": context.plan["comparison_protocol_id"],
        "vina_score_primary_source": "docking_output_pdbqt_remark",
        "score_crosscheck_tolerance_kcal_mol": 0.001,
        "ligand_preparation": {
            "protocol_id": protocol["ligand_preparation"]["protocol_id"],
            "formal_charge": runtime["formal_chemistry"]["formal_charge"],
        },
        "docking_protocol": docking_protocol,
        "expected_total_pose_count": model_count * len(conformers),
        "conformers": conformers,
        **V04_FROZEN_ENSEMBLE_RULES,
        "input_sha256": {
            **{field: file_sha256(path) for field, path in receptor_paths.items()},
            "conformers": conformer_hashes,
        },
    }


def production_execution_adapters(protocol: dict[str, Any]) -> ExecutionAdapters:
    """Bind the V0.4 scientific pipeline to one immutable candidate workspace."""
    from analyze_peptide_ensemble import analyze_peptide_ensemble
    from audit_completed_vina_run import audit_completed_run
    from audit_prepared_vina_input import audit_one
    from pose_metrics_v04 import generate_pose_metrics
    from prepare_standardized_vina import prepare_one
    from run_standardized_vina import run_one

    def prepare(context: CandidateExecutionContext, conformer: str) -> dict[str, Path]:
        workspace = _candidate_workspace(context)
        audit = prepare_one(
            _candidate_spec_path(context), context.data_root, workspace, conformer
        )
        candidate = str(context.plan["candidate_id"])
        run_id = f"{candidate}_{conformer}"
        artifacts = {
            f"input_audit.{conformer}": audit,
            f"source_sdf.{conformer}": workspace / "source_sdf" / f"{run_id}.sdf",
            f"ligand_pdbqt.{conformer}": workspace / "ligands" / f"{run_id}.pdbqt",
            f"vina_config.{conformer}": workspace / "configs" / f"{run_id}.txt",
        }
        if conformer == str(context.plan["conformers"][0]["conformer_name"]):
            artifacts.update(
                {
                    "protocol": workspace / "protocol.json",
                    "receptor_monomer": workspace
                    / "receptors"
                    / f"{context.plan['receptor_id']}.pdb",
                    "receptor_pdbqt": workspace
                    / "receptors"
                    / f"{context.plan['receptor_id']}.pdbqt",
                    "receptor_dimer": workspace
                    / "receptors"
                    / f"{context.plan['template_id']}_TfR1_AB.pdb",
                }
            )
        return artifacts

    def audit_input(context: CandidateExecutionContext, conformer: str) -> dict[str, Path]:
        path = audit_one(
            _candidate_spec_path(context),
            context.data_root,
            _candidate_workspace(context),
            conformer,
        )
        return {f"pre_docking_gate.{conformer}": path}

    def dock(context: CandidateExecutionContext, conformer: str) -> dict[str, Any]:
        workspace = _candidate_workspace(context)
        candidate = str(context.plan["candidate_id"])
        status_path, complete = run_one(
            workspace, candidate, conformer, context.vina_path
        )
        status = _load_json_object(status_path)
        artifacts: dict[str, Path] = {f"run_status.{conformer}": status_path}
        for field, name in (("output_path", "vina_output"), ("log_path", "vina_log")):
            raw = status.get(field)
            if isinstance(raw, str):
                artifacts[f"{name}.{conformer}"] = workspace / raw
        return {
            "complete": complete,
            "exit_code": status.get("exit_code"),
            "model_labels": status.get("model_labels", []),
            "artifacts": artifacts,
        }

    def metrics(context: CandidateExecutionContext, conformer: str) -> dict[str, Path]:
        workspace = _candidate_workspace(context)
        candidate = str(context.plan["candidate_id"])
        run_id = f"{candidate}_{conformer}"
        output = workspace / "analysis" / f"{run_id}_pose_metrics.tsv"
        generate_pose_metrics(
            docking_output=workspace / "outputs" / f"{run_id}_out.pdbqt",
            source_sdf=workspace / "source_sdf" / f"{run_id}.sdf",
            ligand_pdbqt=workspace / "ligands" / f"{run_id}.pdbqt",
            receptor_monomer=workspace / "receptors" / f"{context.plan['receptor_id']}.pdb",
            receptor_dimer=workspace / "receptors" / f"{context.plan['template_id']}_TfR1_AB.pdb",
            receptor_chain="A",
            dimer_chain="B",
            sequence=str(context.plan["sequence"]),
            target_residues=set(context.plan["target_residues"]),
            output_tsv=output,
        )
        return {f"pose_metrics.{conformer}": output}

    def audit_output(context: CandidateExecutionContext, conformer: str) -> dict[str, Path]:
        path = audit_completed_run(
            _candidate_workspace(context),
            str(context.plan["candidate_id"]),
            conformer,
        )
        return {f"post_docking_gate.{conformer}": path}

    def analyze(context: CandidateExecutionContext) -> dict[str, Path]:
        spec = build_ensemble_spec(context, protocol)
        spec_path = _candidate_workspace(context) / "specs" / "ensemble_spec.json"
        _write_json_once(spec_path, spec)
        output = analyze_peptide_ensemble(
            spec_path,
            context.data_root,
            context.attempt_root / "ensemble" / str(context.plan["candidate_id"]),
            allowed_output_root=context.attempt_root,
        )
        return {
            "ensemble_spec": spec_path,
            **{
                f"ensemble.{name}": output / name
                for name in REQUIRED_ANALYSIS_FILES
            },
        }

    return ExecutionAdapters(prepare, audit_input, dock, metrics, audit_output, analyze)


def execute_batch(
    manifest_path: Path,
    protocol_path: Path,
    data_root: Path,
    output_root: Path,
    vina_path: Path,
    *,
    groups: list[str] | None = None,
    candidate_ids: list[str] | None = None,
    resume: bool = False,
) -> dict[str, Any]:
    """Production entry point. It runs Vina only after the complete preflight passes."""
    from batch_preflight_v05 import run_dry_run
    from receptor_registry import resolve_receptor_bundle

    data_root = data_root.resolve()
    output_root = output_root.resolve()
    if not output_root.is_relative_to(data_root):
        return _batch_failure(
            "OUTPUT_ROOT_OUTSIDE_DATA_ROOT",
            "Production output root must be inside --data-root for relative provenance",
        )
    preflight = run_dry_run(
        manifest_path,
        protocol_path,
        data_root,
        output_root,
        vina_path,
        groups=groups,
        candidate_ids=candidate_ids,
        allow_existing_output_root=resume,
    )
    if preflight.get("outcome") == "DRY_RUN_FAIL_GLOBAL":
        return {
            "outcome": "BATCH_FAIL_GLOBAL",
            "failures": preflight.get("failures", []),
            "candidates": [],
            "vina_docking_process_count": 0,
        }
    try:
        manifest_path = manifest_path.resolve()
        protocol_path = protocol_path.resolve()
        manifest = _load_json_object(manifest_path)
        protocol = _load_json_object(protocol_path)
        inventory_path = _resolve_contract_path(
            protocol["contracts"]["ligand_inventory"],
            protocol_path=protocol_path,
            data_root=data_root,
        )
        inventory = _load_json_object(inventory_path)
        registry_path = _resolve_contract_path(
            protocol["contracts"]["receptor_registry"],
            protocol_path=protocol_path,
            data_root=data_root,
        )
        bundles = {
            group: resolve_receptor_bundle(registry_path, data_root, group)
            for group in sorted({str(row["group"]) for row in preflight["plan"]})
        }
        identities = {
            str(plan["candidate_id"]): build_candidate_identity(
                plan,
                manifest=manifest,
                manifest_path=manifest_path,
                protocol=protocol,
                protocol_path=protocol_path,
                inventory=inventory,
                inventory_path=inventory_path,
                receptor_bundle=bundles[str(plan["group"])],
                vina_path=vina_path.resolve(),
            )
            for plan in preflight["plan"]
            if plan.get("preflight_status") == "PASS"
        }
        batch_identity = {
            "schema_version": "0.5b",
            "candidate_manifest_sha256": file_sha256(manifest_path),
            "screening_protocol_sha256": file_sha256(protocol_path),
            "ligand_inventory_sha256": file_sha256(inventory_path),
            "receptor_registry_sha256": file_sha256(registry_path),
            "vina_executable_sha256": file_sha256(vina_path.resolve()),
            "chemistry_contract_sha256": protocol["contracts"][
                "chemistry_contract"
            ]["sha256"],
            "analysis_schema_sha256": protocol["contracts"]["analysis_schema"][
                "sha256"
            ],
            "dependency_lock_sha256": protocol["contracts"]["dependency_lock"][
                "sha256"
            ],
            "scientific_semantics_baseline": protocol[
                "scientific_semantics_baseline"
            ],
        }
        adapters = production_execution_adapters(protocol)
    except Exception as exc:
        return _batch_failure("EXECUTION_CONTRACT_CONSTRUCTION_FAILED", str(exc))
    return execute_planned_batch(
        preflight,
        output_root,
        batch_identity,
        identities,
        adapters,
        data_root=data_root,
        vina_path=vina_path.resolve(),
        resume=resume,
    )
