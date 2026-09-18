"""One-off preparation and dry-run for the 23 ready non-Top12 candidates."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

from standardized_vina_inputs import write_vina_config


DATA_ROOT = Path(r"E:\research_related\competition")
LEGACY_DATA_ROOT = DATA_ROOT / "docking_test_v1_top12"
OUTPUT_ROOT = DATA_ROOT / "docking_test" / "common_receptor_remaining23_v1"
TOP12_STAGE = LEGACY_DATA_ROOT / "common_receptor_top12_v1"
MANIFEST_PATH = DATA_ROOT / "vina_analysis_automation" / "specs" / "top40_candidates_v05.json"
INVENTORY_PATH = (
    LEGACY_DATA_ROOT
    / "partial35_exploratory_v05"
    / "inputs"
    / "top40_ligand_inventory_partial35_v05.json"
)
RECEPTOR_PATH = (
    DATA_ROOT
    / "GPT_vina_branch_20260906"
    / "receptor_pdbqt"
    / "h26d3_af3_s2718_m2_TfR1_A.pdbqt"
)

EXPECTED_READY = (
    "A1", "A2", "A5", "A6", "A7", "A10",
    "B2", "B5", "B6", "B8", "B9", "B10",
    "C1", "C3", "C4", "C5", "C6", "C9", "C10",
    "D1", "D2", "D5", "D9",
)
EXPECTED_BLOCKED = ("A8", "B1", "D3", "D7", "D10")
EXPECTED_TOP12 = ("A3", "A9", "A4", "B3", "B7", "B4", "C2", "C7", "C8", "D8", "D4", "D6")
CONFORMERS = ("conf01", "conf02", "conf03")
BOX_CENTER = (-10.827158, 23.888659, -11.674443)
BOX_SIZE = (30.591256, 37.457546, 33.838885)
VINA_SETTINGS = {
    "exhaustiveness": 16,
    "num_modes": 20,
    "energy_range": 5,
    "seed": 1701,
    "cpu": 8,
}
RECEPTOR_SHA256 = "f95e4f3eceb2774a711466c3d5690e14e256122f85d42df6dc459ea820f81ca1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def json_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, ensure_ascii=False) + "\n").encode("utf-8")


def write_json_once(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as handle:
        handle.write(json_bytes(value))


def resolve_moved_path(recorded: str) -> Path:
    normalized = recorded.replace("/", "\\")
    prefix = "docking_test\\"
    if normalized.startswith(prefix):
        normalized = "docking_test_v1_top12\\" + normalized[len(prefix):]
    return (DATA_ROOT / normalized).resolve()


def check_file(path: Path, expected_hash: str, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty {label}: {path}")
    actual = sha256(path)
    if actual != expected_hash:
        raise ValueError(f"{label} SHA256 mismatch: expected {expected_hash}, actual {actual}")


def parse_config(path: Path) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        key, value = line.split("=", 1)
        parsed[key.strip()] = value.strip()
    return parsed


def main() -> int:
    if OUTPUT_ROOT.exists():
        raise FileExistsError(f"Refusing to overwrite existing output root: {OUTPUT_ROOT}")

    top12_protocol_path = TOP12_STAGE / "protocol.json"
    top12_selection_path = TOP12_STAGE / "selection.json"
    top12_protocol = load_json(top12_protocol_path)
    top12_selection = load_json(top12_selection_path)
    manifest = load_json(MANIFEST_PATH)
    inventory = load_json(INVENTORY_PATH)

    if sha256(INVENTORY_PATH) != top12_protocol["ligand_inventory"]["sha256"]:
        raise ValueError("Moved ligand inventory no longer matches the frozen Top12 inventory hash")
    if tuple(row["candidate_id"] for row in top12_selection["candidates"]) != EXPECTED_TOP12:
        raise ValueError("Frozen Top12 candidate identity differs from the approved exclusion list")

    scientific_fields = ("comparison_mode", "receptor", "box", "vina", "model_validation")
    expected_top12 = {
        "comparison_mode": "common_receptor_common_box",
        "receptor": top12_protocol["receptor"],
        "box": {"center_A": list(BOX_CENTER), "size_A": list(BOX_SIZE)},
        "vina": {
            "version": "1.2.7",
            "executable_sha256": "e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5",
            **VINA_SETTINGS,
            "concurrency": 1,
        },
        "model_validation": {
            "minimum_model_count": 1,
            "maximum_model_count": 20,
            "labels": "continuous_1_to_N_without_duplicates",
            "parseable_vina_result_required_per_model": True,
        },
    }
    protocol_match = all(top12_protocol[field] == expected_top12[field] for field in scientific_fields)
    if not protocol_match:
        raise ValueError("Top12 scientific protocol differs from the approved remaining23 protocol")
    check_file(RECEPTOR_PATH, RECEPTOR_SHA256, "common receptor PDBQT")

    canonical = {row["candidate_id"]: row for row in manifest["candidates"]}
    if set(canonical) != set(EXPECTED_READY) | set(EXPECTED_BLOCKED) | set(EXPECTED_TOP12):
        raise ValueError("Top40 manifest partition differs from READY/BLOCKED/Top12 approval")

    slots: dict[tuple[str, str], dict[str, Any]] = {}
    for row in inventory["slots"]:
        key = (row["candidate_id"], row["conformer_name"])
        if key in slots:
            raise ValueError(f"Duplicate inventory slot: {key}")
        slots[key] = row

    readiness_rows: list[dict[str, Any]] = []
    jobs: list[dict[str, Any]] = []
    for candidate_id in EXPECTED_READY:
        candidate = canonical[candidate_id]
        conformer_rows: list[dict[str, Any]] = []
        for conformer in CONFORMERS:
            slot = slots.get((candidate_id, conformer))
            if slot is None:
                raise ValueError(f"Missing ligand inventory slot: {candidate_id}/{conformer}")
            if slot.get("readiness_status") != "approved_for_exploratory_screening":
                raise ValueError(f"Ligand slot is not approved: {candidate_id}/{conformer}")

            sdf_path = (DATA_ROOT / slot["source_sdf_path"]).resolve()
            pdbqt_path = resolve_moved_path(slot["ligand_pdbqt_path"])
            audit_path = resolve_moved_path(slot["audit_path"])
            check_file(sdf_path, slot["source_sdf_sha256"], f"SDF {candidate_id}/{conformer}")
            check_file(pdbqt_path, slot["ligand_pdbqt_sha256"], f"PDBQT {candidate_id}/{conformer}")
            check_file(audit_path, slot["audit_sha256"], f"input audit {candidate_id}/{conformer}")

            audit = load_json(audit_path)
            identity_ok = (
                audit.get("status") == "PASS"
                and audit.get("candidate") == candidate_id
                and audit.get("conformer") == conformer
                and audit.get("source_sdf", {}).get("sequence") == candidate["sequence"]
                and audit.get("source_sdf", {}).get("sha256") == slot["source_sdf_sha256"]
                and audit.get("ligand_pdbqt", {}).get("sha256") == slot["ligand_pdbqt_sha256"]
                and audit.get("ligand_preparation", {}).get("coordinates_regenerated") is False
            )
            if not identity_ok:
                raise ValueError(f"Ligand identity audit failed: {candidate_id}/{conformer}")

            run_id = f"{candidate_id}_{conformer}"
            config_path = OUTPUT_ROOT / "configs" / candidate_id / f"{run_id}.txt"
            log_path = OUTPUT_ROOT / "logs" / candidate_id / f"{run_id}.log"
            output_path = OUTPUT_ROOT / "outputs" / candidate_id / f"{run_id}_out.pdbqt"
            conformer_rows.append(
                {
                    "conformer": conformer,
                    "status": "PASS",
                    "source_sdf": str(sdf_path),
                    "source_sdf_sha256": slot["source_sdf_sha256"],
                    "ligand_pdbqt": str(pdbqt_path),
                    "ligand_pdbqt_sha256": slot["ligand_pdbqt_sha256"],
                    "input_audit": str(audit_path),
                    "input_audit_sha256": slot["audit_sha256"],
                    "recorded_legacy_ligand_path": slot["ligand_pdbqt_path"],
                }
            )
            jobs.append(
                {
                    "job_index": len(jobs) + 1,
                    "candidate_id": candidate_id,
                    "source_group": candidate["group"],
                    "sequence": candidate["sequence"],
                    "conformer_name": conformer,
                    "receptor_pdbqt": str(RECEPTOR_PATH),
                    "receptor_pdbqt_sha256": RECEPTOR_SHA256,
                    "ligand_pdbqt": str(pdbqt_path),
                    "ligand_pdbqt_sha256": slot["ligand_pdbqt_sha256"],
                    "box_center_A": list(BOX_CENTER),
                    "box_size_A": list(BOX_SIZE),
                    "vina_scoring": "vina",
                    "vina_settings": dict(VINA_SETTINGS),
                    "config_path": str(config_path),
                    "log_path": str(log_path),
                    "output_path": str(output_path),
                }
            )
        readiness_rows.append(
            {
                "candidate_id": candidate_id,
                "sequence": candidate["sequence"],
                "status": "READY",
                "conformers": conformer_rows,
            }
        )

    job_keys = [(row["candidate_id"], row["conformer_name"]) for row in jobs]
    candidate_counts = Counter(row["candidate_id"] for row in jobs)
    planned_paths = [Path(row[key]) for row in jobs for key in ("config_path", "log_path", "output_path")]
    preflight = {
        "candidate_count": len(candidate_counts),
        "conformers_per_candidate": sorted(set(candidate_counts.values())),
        "job_count": len(jobs),
        "unique_job_count": len(set(job_keys)),
        "unique_receptor_count": len({row["receptor_pdbqt"] for row in jobs}),
        "unique_receptor_sha256_count": len({row["receptor_pdbqt_sha256"] for row in jobs}),
        "unique_center_count": len({tuple(row["box_center_A"]) for row in jobs}),
        "unique_size_count": len({tuple(row["box_size_A"]) for row in jobs}),
        "unique_vina_parameter_set_count": len(
            {(row["vina_scoring"], tuple(sorted(row["vina_settings"].items()))) for row in jobs}
        ),
        "missing_ligand": sum(not Path(row["ligand_pdbqt"]).is_file() for row in jobs),
        "duplicate_job": len(jobs) - len(set(job_keys)),
        "preexisting_planned_path_collision": sum(path.exists() for path in planned_paths),
    }
    expected_preflight = {
        "candidate_count": 23,
        "conformers_per_candidate": [3],
        "job_count": 69,
        "unique_job_count": 69,
        "unique_receptor_count": 1,
        "unique_receptor_sha256_count": 1,
        "unique_center_count": 1,
        "unique_size_count": 1,
        "unique_vina_parameter_set_count": 1,
        "missing_ligand": 0,
        "duplicate_job": 0,
        "preexisting_planned_path_collision": 0,
    }
    if preflight != expected_preflight:
        raise ValueError(f"Preflight job audit failed: {preflight}")

    selection = {
        "manifest_id": "tfr1_common_receptor_remaining23_selection_v1",
        "candidate_count": 23,
        "conformers": list(CONFORMERS),
        "candidates": [
            {
                "candidate_id": candidate_id,
                "source_group": canonical[candidate_id]["group"],
                "sequence": canonical[candidate_id]["sequence"],
            }
            for candidate_id in EXPECTED_READY
        ],
        "blocked_candidates_not_planned": list(EXPECTED_BLOCKED),
        "excluded_completed_top12": list(EXPECTED_TOP12),
    }
    selection_hash = hashlib.sha256(json_bytes(selection)).hexdigest()
    protocol = {
        "protocol_id": "tfr1_common_receptor_remaining23_v1",
        "derived_from_top12": {
            "protocol_path": str(top12_protocol_path),
            "protocol_sha256": sha256(top12_protocol_path),
            "protocol_match_fields": list(scientific_fields),
            "status": "PASS",
        },
        "comparison_mode": top12_protocol["comparison_mode"],
        "receptor": top12_protocol["receptor"],
        "box": top12_protocol["box"],
        "vina": {**top12_protocol["vina"], "effective_scoring": "vina"},
        "model_validation": top12_protocol["model_validation"],
        "selection": {
            "path": str(OUTPUT_ROOT / "selection.json"),
            "sha256": selection_hash,
        },
        "ligand_inventory": {
            "recorded_legacy_path": top12_protocol["ligand_inventory"]["path"],
            "resolved_current_path": str(INVENTORY_PATH),
            "sha256": sha256(INVENTORY_PATH),
            "path_migration_only": True,
        },
        "output": {"root": str(OUTPUT_ROOT), "no_overwrite": True},
        "dry_run_only": True,
    }

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=False)
    write_json_once(OUTPUT_ROOT / "selection.json", selection)
    write_json_once(OUTPUT_ROOT / "protocol.json", protocol)
    write_json_once(
        OUTPUT_ROOT / "readiness_audit.json",
        {
            "total_remaining": 28,
            "ready": 23,
            "blocked": 5,
            "path_migration_note": "Legacy docking_test paths were resolved to the existing docking_test_v1_top12 tree; frozen metadata was not modified.",
            "ready_candidates": readiness_rows,
            "blocked_candidates_not_modified": list(EXPECTED_BLOCKED),
        },
    )
    for job in jobs:
        write_vina_config(
            Path(job["config_path"]),
            job["receptor_pdbqt"],
            job["ligand_pdbqt"],
            {"center_A": job["box_center_A"], "size_A": job["box_size_A"]},
            job["vina_settings"],
        )
    job_plan_path = OUTPUT_ROOT / "remaining23_job_plan.json"
    write_json_once(job_plan_path, {"jobs": jobs})

    configs = [Path(row["config_path"]) for row in jobs]
    parsed_configs = [parse_config(path) for path in configs]
    config_parameter_keys = (
        "receptor", "center_x", "center_y", "center_z", "size_x", "size_y", "size_z",
        "exhaustiveness", "num_modes", "energy_range", "seed", "cpu",
    )
    postflight = {
        **preflight,
        "config_count": len(configs),
        "missing_config": sum(not path.is_file() for path in configs),
        "unique_config_parameter_set_count": len(
            {tuple((key, config[key]) for key in config_parameter_keys) for config in parsed_configs}
        ),
        "output_collision": sum(Path(row["output_path"]).exists() for row in jobs),
        "log_collision": sum(Path(row["log_path"]).exists() for row in jobs),
    }
    passed = (
        postflight["candidate_count"] == 23
        and postflight["conformers_per_candidate"] == [3]
        and postflight["job_count"] == 69
        and postflight["unique_job_count"] == 69
        and postflight["unique_receptor_count"] == 1
        and postflight["unique_center_count"] == 1
        and postflight["unique_size_count"] == 1
        and postflight["unique_vina_parameter_set_count"] == 1
        and postflight["config_count"] == 69
        and postflight["missing_ligand"] == 0
        and postflight["missing_config"] == 0
        and postflight["duplicate_job"] == 0
        and postflight["output_collision"] == 0
        and postflight["log_collision"] == 0
        and postflight["unique_config_parameter_set_count"] == 1
        and protocol_match
    )
    report = {
        "state": "DRY_RUN_COMPLETE/PREPARED" if passed else "DRY_RUN_FAILED",
        "ready_for_real_docking": passed,
        "candidates": 23,
        "planned_jobs": 69,
        "protocol_match_with_top12": "PASS" if protocol_match else "FAIL",
        "vina_process_count": 0,
        "vina_invoked": False,
        "audit": postflight,
        "hashes": {
            "selection": sha256(OUTPUT_ROOT / "selection.json"),
            "protocol": sha256(OUTPUT_ROOT / "protocol.json"),
            "readiness_audit": sha256(OUTPUT_ROOT / "readiness_audit.json"),
            "job_plan": sha256(job_plan_path),
            "configs": {str(path): sha256(path) for path in configs},
        },
        "blockers": [] if passed else ["One or more dry-run checks failed"],
    }
    write_json_once(OUTPUT_ROOT / "dry_run_report.json", report)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
