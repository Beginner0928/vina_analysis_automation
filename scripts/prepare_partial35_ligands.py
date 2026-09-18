"""Prepare and publish the immutable PARTIAL35 exploratory ligand PDBQT package."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import meeko

from audit_partial35_conformers import audit_partial35_package
from candidate_manifest_v05 import load_candidate_manifest, validate_ligand_inventory
from competition_conformer_core import sha256_file
from partial35_contract_v05 import EXCLUDED_CANDIDATES, PARTIAL35_CANDIDATE_IDS
from standardized_vina_inputs import (
    audit_ligand_pdbqt,
    audit_source_sdf,
    compare_pdbqt_coordinates,
    prepare_backbone_rigid_ligand,
)


LIGAND_PACKAGE_RELATIVE = Path("docking_test") / "partial35_exploratory_v05" / "inputs"
DERIVED_INVENTORY_FILENAME = "top40_ligand_inventory_partial35_v05.json"
DERIVED_PROTOCOL_FILENAME = "screening_protocol_partial35_v05.json"
PREPARATION_MANIFEST_FILENAME = "ligand_preparation_manifest.json"
LIGAND_PREPARATION_PROTOCOL_ID = "gpt_manual_standardized_ph74_meeko080_backbone_rigid_v1"


def _write_json_exclusive(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def _relative_to_data_root(path: Path, data_root: Path) -> str:
    resolved = path.resolve()
    root = data_root.resolve()
    if not resolved.is_relative_to(root):
        raise ValueError(f"Artifact must be inside --data-root: {resolved}")
    return resolved.relative_to(root).as_posix()


def ligand_preparation_contract_sha256(protocol: dict[str, Any]) -> str:
    payload = json.dumps(
        protocol["ligand_preparation"],
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def build_partial35_inventory(
    canonical_candidates: list[dict[str, Any]],
    artifact_records: list[dict[str, Any]],
    *,
    candidate_manifest_sha256: str,
    partial_package_manifest_sha256: str,
    partial_package_manifest_path: str,
    ligand_preparation_contract_sha256: str,
) -> dict[str, Any]:
    expected_included = {
        (candidate_id, conformer)
        for candidate_id in PARTIAL35_CANDIDATE_IDS
        for conformer in ("conf01", "conf02", "conf03")
    }
    by_key: dict[tuple[str, str], dict[str, Any]] = {}
    for record in artifact_records:
        key = (str(record.get("candidate_id")), str(record.get("conformer_name")))
        if key in by_key:
            raise ValueError(f"Duplicate PARTIAL35 ligand artifact record: {key}")
        by_key[key] = record
    missing = sorted(expected_included - set(by_key))
    extra = sorted(set(by_key) - expected_included)
    if missing or extra:
        raise ValueError(f"PARTIAL35 ligand artifact set mismatch; missing={missing}, extra={extra}")

    excluded_ids = {candidate_id for candidate_id, _reason in EXCLUDED_CANDIDATES}
    slots: list[dict[str, Any]] = []
    for candidate in canonical_candidates:
        candidate_id = str(candidate["candidate_id"])
        for conformer in ("conf01", "conf02", "conf03"):
            if candidate_id in excluded_ids:
                slots.append(
                    {
                        "candidate_id": candidate_id,
                        "conformer_name": conformer,
                        "readiness_status": "conformer_generation_qc_fail_not_docked",
                        "exclusion_reason": "CONFORMER_GENERATION_QC_FAIL",
                        "screening_status": "CONFORMER_GENERATION_QC_FAIL_NOT_DOCKED",
                    }
                )
                continue
            record = by_key[(candidate_id, conformer)]
            slots.append(
                {
                    **record,
                    "readiness_status": "approved_for_exploratory_screening",
                    "approval_basis": "PARTIAL35_SDF_AND_PDBQT_AUDIT_PASS",
                }
            )
    inventory = {
        "schema_version": "0.5",
        "inventory_id": "tfr1_partial35_exploratory_ligands_v05",
        "candidate_manifest_sha256": candidate_manifest_sha256,
        "expected_slot_count": len(canonical_candidates) * 3,
        "approved_exploratory_slot_count": len(expected_included),
        "qc_excluded_slot_count": len(excluded_ids) * 3,
        "partial_package": {
            "package_id": "competition_v1_partial35_mmff_qc_pass",
            "manifest_path": partial_package_manifest_path,
            "manifest_sha256": partial_package_manifest_sha256,
            "completeness_statement": "NOT_A_COMPLETE_TOP40_SCREEN",
        },
        "ligand_preparation_contract_sha256": ligand_preparation_contract_sha256,
        "slots": slots,
    }
    validate_ligand_inventory(
        inventory,
        canonical_candidates,
        {
            "formal_conformers": {
                "required_conformer_names": ["conf01", "conf02", "conf03"]
            }
        },
    )
    return inventory


def build_partial35_screening_protocol(
    base_protocol: dict[str, Any],
    *,
    inventory_path: Path,
    inventory_sha256: str,
    candidate_manifest_path: str | None = None,
) -> dict[str, Any]:
    derived = copy.deepcopy(base_protocol)
    derived["protocol_id"] = "tfr1_partial35_exploratory_screening_v05"
    derived["contracts"]["ligand_inventory"] = {
        "path": inventory_path.as_posix(),
        "path_base": "protocol_directory",
        "sha256": inventory_sha256,
    }
    if candidate_manifest_path is not None:
        derived["contracts"]["candidate_manifest"] = {
            **derived["contracts"]["candidate_manifest"],
            "path": candidate_manifest_path,
            "path_base": "data_root",
        }
    return derived


def _formal_charge(sequence: str) -> int:
    return sequence.count("R") + sequence.count("K") - sequence.count("D") - sequence.count("E")


def prepare_partial35_ligands(
    *,
    package_directory: Path,
    partial_contract_path: Path,
    generation_protocol_path: Path,
    canonical_manifest_path: Path,
    base_screening_protocol_path: Path,
    data_root: Path,
    output_root: Path,
) -> Path:
    """Prepare all 105 ligands, then publish inventory and derived protocol."""
    expected_output = data_root.resolve() / LIGAND_PACKAGE_RELATIVE
    if output_root.resolve() != expected_output:
        raise ValueError(
            f"PARTIAL35 ligand package must use the isolated frozen path: {expected_output}"
        )
    if output_root.exists():
        raise FileExistsError(f"Refusing to overwrite existing ligand package: {output_root}")
    gate = audit_partial35_package(
        package_directory,
        partial_contract_path,
        generation_protocol_path,
        data_root,
    )
    if gate.get("status") != "PASS" or gate.get("sdf_count") != 105:
        raise ValueError("PARTIAL35 conformer package did not pass Gate B")
    canonical_manifest, candidates = load_candidate_manifest(canonical_manifest_path)
    candidate_by_id = {row["candidate_id"]: row for row in candidates}
    core_manifest = json.loads(
        (package_directory / "preparation_manifest.json").read_text(encoding="utf-8")
    )
    partial_manifest_path = package_directory / "partial_package_manifest.json"
    partial_manifest = json.loads(partial_manifest_path.read_text(encoding="utf-8"))
    base_protocol = json.loads(base_screening_protocol_path.read_text(encoding="utf-8"))
    prep_contract_hash = ligand_preparation_contract_sha256(base_protocol)
    output_root.mkdir(parents=True)
    records: list[dict[str, Any]] = []
    for candidate_row in core_manifest["candidates"]:
        candidate_id = str(candidate_row["candidate_id"])
        if candidate_id not in PARTIAL35_CANDIDATE_IDS:
            raise ValueError(f"Unexpected candidate in PARTIAL35 package: {candidate_id}")
        canonical = candidate_by_id[candidate_id]
        for conformer_row in candidate_row["formal_conformers"]:
            conformer = str(conformer_row["formal_conformer_ordinal"])
            run_id = f"{candidate_id}_{conformer}"
            source_sdf = package_directory / conformer_row["path"]
            if sha256_file(source_sdf) != conformer_row["sha256"]:
                raise ValueError(f"Source SDF hash mismatch before PDBQT preparation: {run_id}")
            source_audit = audit_source_sdf(
                source_sdf,
                candidate_id,
                canonical["sequence"],
                _formal_charge(canonical["sequence"]),
            )
            pdbqt = output_root / "ligands" / candidate_id / f"{run_id}.pdbqt"
            preparation = prepare_backbone_rigid_ligand(source_sdf, pdbqt)
            pdbqt_audit = audit_ligand_pdbqt(pdbqt, source_sdf, canonical["sequence"])
            coordinate = compare_pdbqt_coordinates(source_sdf, pdbqt)
            audit_payload = {
                "schema_version": "0.5-partial35-ligand",
                "status": "PASS",
                "candidate": candidate_id,
                "conformer": conformer,
                "run_id": run_id,
                "source_sdf": source_audit,
                "generation_provenance": {
                    "protocol_id": partial_manifest["generation_protocol_id"],
                    "protocol_sha256": partial_manifest["generation_protocol_sha256"],
                    "partial_package_id": partial_manifest["package_id"],
                    "partial_package_manifest_sha256": sha256_file(partial_manifest_path),
                },
                "ligand_preparation": {
                    "protocol": base_protocol["ligand_preparation"],
                    "contract_sha256": prep_contract_hash,
                    "meeko_version": meeko.__version__,
                    **preparation,
                },
                "ligand_pdbqt": {
                    **pdbqt_audit,
                    "coordinate_identity": coordinate,
                },
                "coordinates_regenerated": False,
            }
            audit_path = output_root / "audit" / f"{run_id}_ligand_audit.json"
            _write_json_exclusive(audit_path, audit_payload)
            records.append(
                {
                    "candidate_id": candidate_id,
                    "conformer_name": conformer,
                    "source_sdf_path": _relative_to_data_root(source_sdf, data_root),
                    "source_sdf_sha256": source_audit["sha256"],
                    "audit_path": _relative_to_data_root(audit_path, data_root),
                    "audit_sha256": sha256_file(audit_path),
                    "ligand_pdbqt_path": _relative_to_data_root(pdbqt, data_root),
                    "ligand_pdbqt_sha256": pdbqt_audit["sha256"],
                    "pdbqt_audit_path": _relative_to_data_root(audit_path, data_root),
                    "pdbqt_audit_sha256": sha256_file(audit_path),
                    "generation_provenance": audit_payload["generation_provenance"],
                }
            )
    inventory = build_partial35_inventory(
        candidates,
        records,
        candidate_manifest_sha256=sha256_file(canonical_manifest_path),
        partial_package_manifest_sha256=sha256_file(partial_manifest_path),
        partial_package_manifest_path=_relative_to_data_root(
            partial_manifest_path, data_root
        ),
        ligand_preparation_contract_sha256=prep_contract_hash,
    )
    inventory_path = output_root / DERIVED_INVENTORY_FILENAME
    _write_json_exclusive(inventory_path, inventory)
    derived_protocol = build_partial35_screening_protocol(
        base_protocol,
        inventory_path=Path(DERIVED_INVENTORY_FILENAME),
        inventory_sha256=sha256_file(inventory_path),
        candidate_manifest_path=_relative_to_data_root(
            canonical_manifest_path, data_root
        ),
    )
    derived_protocol_path = output_root / DERIVED_PROTOCOL_FILENAME
    _write_json_exclusive(derived_protocol_path, derived_protocol)
    hash_inventory_path = output_root / "PDBQT_SHA256SUMS.tsv"
    with hash_inventory_path.open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["candidate_id", "conformer_name", "path", "sha256"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(
            {
                "candidate_id": row["candidate_id"],
                "conformer_name": row["conformer_name"],
                "path": row["ligand_pdbqt_path"],
                "sha256": row["ligand_pdbqt_sha256"],
            }
            for row in records
        )
    manifest = {
        "schema_version": "0.5-partial35-ligand-package",
        "package_id": "tfr1_partial35_exploratory_ligand_inputs_v05",
        "package_classification": "PARTIAL_EXPLORATORY_SCREENING",
        "completeness_statement": "NOT_A_COMPLETE_TOP40_SCREEN",
        "partial_conformer_package_id": partial_manifest["package_id"],
        "partial_conformer_package_manifest_sha256": sha256_file(partial_manifest_path),
        "candidate_manifest_id": canonical_manifest["manifest_id"],
        "candidate_manifest_sha256": sha256_file(canonical_manifest_path),
        "ligand_preparation_protocol": base_protocol["ligand_preparation"],
        "ligand_preparation_contract_sha256": prep_contract_hash,
        "meeko_version": meeko.__version__,
        "candidate_count": 35,
        "pdbqt_count": len(records),
        "inventory_path": inventory_path.name,
        "inventory_sha256": sha256_file(inventory_path),
        "derived_screening_protocol_path": derived_protocol_path.name,
        "derived_screening_protocol_sha256": sha256_file(derived_protocol_path),
        "pdbqt_hash_inventory_path": hash_inventory_path.name,
        "pdbqt_hash_inventory_sha256": sha256_file(hash_inventory_path),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS",
    }
    manifest_path = output_root / PREPARATION_MANIFEST_FILENAME
    _write_json_exclusive(manifest_path, manifest)
    return manifest_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partial-package", type=Path, required=True)
    parser.add_argument("--partial-contract", type=Path, required=True)
    parser.add_argument("--generation-protocol", type=Path, required=True)
    parser.add_argument("--candidate-manifest", type=Path, required=True)
    parser.add_argument("--screening-protocol", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        data_root = args.data_root.resolve()
        manifest = prepare_partial35_ligands(
            package_directory=args.partial_package.resolve(),
            partial_contract_path=args.partial_contract.resolve(),
            generation_protocol_path=args.generation_protocol.resolve(),
            canonical_manifest_path=args.candidate_manifest.resolve(),
            base_screening_protocol_path=args.screening_protocol.resolve(),
            data_root=data_root,
            output_root=data_root / LIGAND_PACKAGE_RELATIVE,
        )
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print(f"PARTIAL35 ligand preparation PASS: {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
