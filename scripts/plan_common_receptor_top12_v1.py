"""Prepare the fixed-receptor Top12 comparison package without running Vina."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from receptor_registry import resolve_receptor_bundle
from standardized_vina_inputs import write_vina_config


EXPECTED_CONFORMERS = ("conf01", "conf02", "conf03")
EXPECTED_JOB_COUNT = 36
CONFIG_SETTING_KEYS = ("exhaustiveness", "num_modes", "energy_range", "seed", "cpu")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_contract(path: Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def assert_new_output_root(output_root: Path) -> None:
    if Path(output_root).exists():
        raise FileExistsError(f"Refusing to reuse or overwrite output root: {output_root}")


def resolve_moved_data_file(data_root: Path, recorded_path: str) -> Path:
    """Resolve the approved docking_test -> docking_test_v1_top12 path migration."""
    relative = Path(recorded_path.replace("\\", "/"))
    direct = (data_root / relative).resolve()
    if direct.is_file():
        return direct
    if relative.parts and relative.parts[0] == "docking_test":
        migrated = (data_root / "docking_test_v1_top12" / Path(*relative.parts[1:])).resolve()
        if migrated.is_file():
            return migrated
    return direct


def _resolve_locked_file(data_root: Path, specification: Mapping[str, Any], label: str) -> Path:
    relative = specification.get("path")
    expected = specification.get("sha256")
    if specification.get("path_base") != "data_root" or not isinstance(relative, str):
        raise ValueError(f"{label} must be data-root-relative")
    path = (data_root / relative).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")
    actual = sha256(path)
    if actual != expected:
        raise ValueError(f"{label} SHA-256 mismatch: expected {expected}, actual {actual}")
    return path


def _resolve_receptor(protocol: Mapping[str, Any], data_root: Path) -> dict[str, Any]:
    registry_spec = protocol["receptor_registry"]
    registry_path = _resolve_locked_file(data_root, registry_spec, "receptor registry")
    receptor_spec = protocol["receptor"]
    if receptor_spec.get("registry_group") != "A" or receptor_spec.get("use_for_all_candidates") is not True:
        raise ValueError("Common-receptor protocol must lock registry group A for all candidates")
    bundle = resolve_receptor_bundle(registry_path, data_root, "A")
    if bundle["template_id"] != receptor_spec.get("template_id"):
        raise ValueError("Common receptor template differs from the registry")
    field_pairs = {
        "receptor_pdb": "receptor_monomer_pdb",
        "receptor_pdbqt": "receptor_pdbqt",
        "receptor_dimer_pdb": "receptor_dimer_pdb",
    }
    paths: dict[str, Path] = {}
    for protocol_field, bundle_field in field_pairs.items():
        expected_path = (data_root / receptor_spec[protocol_field]).resolve()
        actual_path = Path(bundle[bundle_field]).resolve()
        if actual_path != expected_path:
            raise ValueError(f"Common receptor {protocol_field} path differs from registry group A")
        expected_hash = receptor_spec["sha256"][protocol_field]
        actual_hash = bundle["sha256"][bundle_field]
        if actual_hash != expected_hash:
            raise ValueError(f"Common receptor {protocol_field} SHA-256 differs from registry group A")
        paths[protocol_field] = actual_path
    return {"bundle": bundle, "paths": paths}


def resolve_jobs(
    *,
    protocol: Mapping[str, Any],
    selection: Mapping[str, Any],
    inventory: Mapping[str, Any],
    data_root: Path,
    output_root: Path,
) -> list[dict[str, Any]]:
    data_root = Path(data_root).resolve()
    output_root = Path(output_root).resolve()
    if protocol.get("protocol_id") != "tfr1_top12_common_receptor_v1":
        raise ValueError("Unexpected common-receptor protocol identity")
    if selection.get("protocol_id") != protocol["protocol_id"]:
        raise ValueError("Selection protocol identity mismatch")
    candidates = selection.get("candidates")
    conformers = selection.get("conformers")
    if not isinstance(candidates, list) or len(candidates) != 12:
        raise ValueError("Selection must contain exactly 12 candidates")
    if conformers != list(EXPECTED_CONFORMERS):
        raise ValueError("Selection must contain exactly conf01, conf02, conf03")

    receptor = _resolve_receptor(protocol, data_root)
    candidate_manifest_path = _resolve_locked_file(
        data_root, protocol["candidate_manifest"], "candidate manifest"
    )
    candidate_manifest = load_contract(candidate_manifest_path)
    canonical = {row["candidate_id"]: row for row in candidate_manifest["candidates"]}

    inventory_slots: dict[tuple[str, str], Mapping[str, Any]] = {}
    for row in inventory.get("slots", []):
        candidate_id = row.get("candidate_id")
        conformer = row.get("conformer_name")
        key = (candidate_id, conformer)
        if key in inventory_slots:
            raise ValueError(f"Duplicate ligand inventory slot: {key}")
        inventory_slots[key] = row

    box = protocol["box"]
    vina = protocol["vina"]
    settings = {key: vina[key] for key in CONFIG_SETTING_KEYS}
    jobs: list[dict[str, Any]] = []
    seen_candidates: set[str] = set()
    for candidate in candidates:
        candidate_id = candidate.get("candidate_id")
        source_group = candidate.get("source_group")
        sequence = candidate.get("sequence")
        if not isinstance(candidate_id, str) or candidate_id in seen_candidates:
            raise ValueError(f"Duplicate or invalid selected candidate: {candidate_id!r}")
        seen_candidates.add(candidate_id)
        canonical_row = canonical.get(candidate_id)
        if canonical_row is None:
            raise ValueError(f"Selected candidate absent from canonical manifest: {candidate_id}")
        if canonical_row.get("group") != source_group or canonical_row.get("sequence") != sequence:
            raise ValueError(f"Selected candidate identity mismatch: {candidate_id}")
        for conformer in conformers:
            slot = inventory_slots.get((candidate_id, conformer))
            if slot is None:
                raise ValueError(f"Missing ligand inventory slot: {candidate_id}/{conformer}")
            if slot.get("readiness_status") != "approved_for_exploratory_screening":
                raise ValueError(f"Ligand slot is not approved: {candidate_id}/{conformer}")
            ligand_path = resolve_moved_data_file(data_root, slot["ligand_pdbqt_path"])
            if not ligand_path.is_file():
                raise FileNotFoundError(f"Missing ligand PDBQT: {ligand_path}")
            ligand_hash = sha256(ligand_path)
            if ligand_hash != slot.get("ligand_pdbqt_sha256"):
                raise ValueError(
                    f"ligand PDBQT SHA-256 mismatch for {candidate_id}/{conformer}: "
                    f"expected {slot.get('ligand_pdbqt_sha256')}, actual {ligand_hash}"
                )
            jobs.append(
                {
                    "job_index": len(jobs) + 1,
                    "candidate_id": candidate_id,
                    "source_group": source_group,
                    "sequence": sequence,
                    "conformer_name": conformer,
                    "receptor_template_id": receptor["bundle"]["template_id"],
                    "receptor_pdb": str(receptor["paths"]["receptor_pdb"]),
                    "receptor_pdbqt": str(receptor["paths"]["receptor_pdbqt"]),
                    "receptor_dimer_pdb": str(receptor["paths"]["receptor_dimer_pdb"]),
                    "receptor_pdbqt_sha256": protocol["receptor"]["sha256"]["receptor_pdbqt"],
                    "ligand_pdbqt": str(ligand_path),
                    "ligand_pdbqt_sha256": ligand_hash,
                    "box_center_A": list(box["center_A"]),
                    "box_size_A": list(box["size_A"]),
                    "vina_settings": dict(settings),
                    "config_path": str(
                        output_root / "configs" / candidate_id / f"{candidate_id}_{conformer}.txt"
                    ),
                }
            )
    return jobs


def validate_jobs(jobs: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    job_keys = [(row["candidate_id"], row["conformer_name"]) for row in jobs]
    candidates = {key[0] for key in job_keys}
    conformers = {key[1] for key in job_keys}
    counts = {
        "job_count": len(jobs),
        "candidate_count": len(candidates),
        "conformer_count": len(conformers),
        "unique_job_count": len(set(job_keys)),
        "unique_receptor_pdbqt_count": len({row["receptor_pdbqt"] for row in jobs}),
        "unique_receptor_pdbqt_sha256_count": len(
            {row["receptor_pdbqt_sha256"] for row in jobs}
        ),
        "unique_center_count": len({tuple(row["box_center_A"]) for row in jobs}),
        "unique_size_count": len({tuple(row["box_size_A"]) for row in jobs}),
        "unique_vina_parameter_set_count": len(
            {tuple(sorted(row["vina_settings"].items())) for row in jobs}
        ),
    }
    coverage = {
        candidate: sorted(row["conformer_name"] for row in jobs if row["candidate_id"] == candidate)
        for candidate in sorted(candidates)
    }
    expected = (
        counts["job_count"] == EXPECTED_JOB_COUNT
        and counts["candidate_count"] == 12
        and counts["conformer_count"] == 3
        and counts["unique_job_count"] == EXPECTED_JOB_COUNT
        and counts["unique_receptor_pdbqt_count"] == 1
        and counts["unique_receptor_pdbqt_sha256_count"] == 1
        and counts["unique_center_count"] == 1
        and counts["unique_size_count"] == 1
        and counts["unique_vina_parameter_set_count"] == 1
        and all(names == list(EXPECTED_CONFORMERS) for names in coverage.values())
    )
    return {**counts, "coverage": coverage, "status": "PASS" if expected else "FAIL"}


def _write_json_once(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def _write_md_summary(path: Path, protocol: Mapping[str, Any], selection: Mapping[str, Any]) -> None:
    receptor = protocol["receptor"]
    box = protocol["box"]
    lines = [
        "# Top12 common-receptor docking inputs",
        "",
        "## Receptor",
        "",
        f"- PDBQT: `{receptor['receptor_pdbqt']}`",
        f"- PDB: `{receptor['receptor_pdb']}`",
        f"- AB dimer PDB: `{receptor['receptor_dimer_pdb']}`",
        f"- Template: `{receptor['template_id']}`",
        "",
        "## Common box",
        "",
        f"- Center: `{box['center_A']}`",
        f"- Size: `{box['size_A']}`",
        "",
        "## Peptides",
        "",
        "| Candidate | Original group | Sequence | Conformers |",
        "|---|---|---|---|",
    ]
    for row in selection["candidates"]:
        lines.append(
            f"| {row['candidate_id']} | {row['source_group']} | `{row['sequence']}` | conf01, conf02, conf03 |"
        )
    lines.extend(
        [
            "",
            "## MD decision note",
            "",
            "This package is a dry-run input package only. MD prioritization should use completed docking poses and post-docking QC; it must not be inferred from preparation readiness alone.",
            "",
        ]
    )
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines))


def materialize_dry_run(
    *,
    protocol_path: Path,
    selection_path: Path,
    inventory_path: Path,
    data_root: Path,
    output_root: Path,
) -> dict[str, Any]:
    output_root = Path(output_root).resolve()
    assert_new_output_root(output_root)
    protocol = load_contract(protocol_path)
    selection = load_contract(selection_path)
    if sha256(Path(selection_path)) != protocol["selection"]["sha256"]:
        raise ValueError("Selection SHA-256 differs from common-receptor protocol")
    if sha256(Path(inventory_path)) != protocol["ligand_inventory"]["sha256"]:
        raise ValueError("Ligand inventory SHA-256 differs from common-receptor protocol")
    inventory = load_contract(inventory_path)
    jobs = resolve_jobs(
        protocol=protocol,
        selection=selection,
        inventory=inventory,
        data_root=data_root,
        output_root=output_root,
    )
    audit = validate_jobs(jobs)
    if audit["status"] != "PASS":
        raise ValueError(f"Projected job audit failed: {audit}")

    output_root.mkdir(parents=True, exist_ok=False)
    for job in jobs:
        write_vina_config(
            Path(job["config_path"]),
            job["receptor_pdbqt"],
            job["ligand_pdbqt"],
            {"center_A": job["box_center_A"], "size_A": job["box_size_A"]},
            job["vina_settings"],
        )
    _write_json_once(output_root / "protocol.json", protocol)
    _write_json_once(output_root / "selection.json", selection)
    _write_json_once(output_root / "top12_job_plan.json", {"jobs": jobs})
    report = {
        "state": "DRY_RUN_COMPLETE/PREPARED",
        "protocol_id": protocol["protocol_id"],
        "protocol_sha256": sha256(Path(protocol_path)),
        "selection_sha256": sha256(Path(selection_path)),
        "ligand_inventory_sha256": sha256(Path(inventory_path)),
        "job_plan_sha256": sha256(output_root / "top12_job_plan.json"),
        "job_audit": audit,
        "vina_docking_process_count": 0,
        "vina_invoked": False,
        "ready_for_real_docking": True,
    }
    _write_json_once(output_root / "dry_run_report.json", report)
    _write_md_summary(output_root / "md_decision_summary.md", protocol, selection)
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--ligand-inventory", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.dry_run:
        raise ValueError("This command supports dry-run only and never executes Vina")
    report = materialize_dry_run(
        protocol_path=args.protocol,
        selection_path=args.selection,
        inventory_path=args.ligand_inventory,
        data_root=args.data_root,
        output_root=args.output_root,
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
