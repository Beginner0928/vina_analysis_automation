"""Gate C audit for the immutable PARTIAL35 exploratory PDBQT package."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

from audit_partial35_conformers import audit_partial35_package
from candidate_manifest_v05 import (
    load_candidate_manifest,
    validate_json_schema,
    validate_ligand_inventory,
    validate_screening_protocol,
)
from competition_conformer_core import sha256_file
from partial35_contract_v05 import EXCLUDED_CANDIDATES, PARTIAL35_CANDIDATE_IDS
from prepare_partial35_ligands import (
    DERIVED_INVENTORY_FILENAME,
    DERIVED_PROTOCOL_FILENAME,
    PREPARATION_MANIFEST_FILENAME,
    ligand_preparation_contract_sha256,
)
from standardized_vina_inputs import (
    audit_ligand_pdbqt,
    audit_source_sdf,
    compare_pdbqt_coordinates,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _resolve_under_data_root(data_root: Path, raw: str, label: str) -> Path:
    value = Path(raw)
    if value.is_absolute():
        raise ValueError(f"{label} must be relative to --data-root")
    root = data_root.resolve()
    path = (root / value).resolve()
    if not path.is_relative_to(root):
        raise ValueError(f"{label} escapes --data-root")
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")
    return path


def _expected_formal_charge(sequence: str) -> int:
    return sequence.count("R") + sequence.count("K") - sequence.count("D") - sequence.count("E")


def validate_partial35_inventory_contract(
    inventory: dict[str, Any],
    canonical_candidates: list[dict[str, Any]],
    protocol: dict[str, Any],
    *,
    expected_partial_package_manifest_sha256: str,
) -> dict[str, Any]:
    slots = validate_ligand_inventory(inventory, canonical_candidates, protocol)
    approved = [
        row
        for row in slots.values()
        if row["readiness_status"] == "approved_for_exploratory_screening"
    ]
    excluded = [
        row
        for row in slots.values()
        if row["readiness_status"] == "conformer_generation_qc_fail_not_docked"
    ]
    if len(approved) != 105 or inventory.get("approved_exploratory_slot_count") != 105:
        raise ValueError("PARTIAL35 inventory must contain exactly 105 approved PDBQT slots")
    if {
        row["candidate_id"] for row in approved
    } != set(PARTIAL35_CANDIDATE_IDS):
        raise ValueError("PARTIAL35 approved PDBQT candidate set mismatch")
    expected_excluded = {candidate_id for candidate_id, _reason in EXCLUDED_CANDIDATES}
    if len(excluded) != 15 or {row["candidate_id"] for row in excluded} != expected_excluded:
        raise ValueError("PARTIAL35 QC-excluded slot set must be the exact five candidates x three")
    if inventory.get("qc_excluded_slot_count") != 15:
        raise ValueError("PARTIAL35 qc_excluded_slot_count must equal 15")
    partial = inventory.get("partial_package", {})
    if partial != {
        "package_id": "competition_v1_partial35_mmff_qc_pass",
        "manifest_path": partial.get("manifest_path"),
        "manifest_sha256": expected_partial_package_manifest_sha256,
        "completeness_statement": "NOT_A_COMPLETE_TOP40_SCREEN",
    }:
        raise ValueError("PARTIAL35 inventory conformer-package provenance mismatch")
    if not isinstance(partial.get("manifest_path"), str) or Path(
        partial["manifest_path"]
    ).is_absolute():
        raise ValueError("PARTIAL35 inventory manifest_path must be relative")
    return {
        "status": "PASS",
        "approved_pdbqt_count": len(approved),
        "qc_excluded_slot_count": len(excluded),
    }


def audit_partial35_ligand_package(
    *,
    ligand_package_root: Path,
    conformer_package_root: Path,
    partial_contract_path: Path,
    generation_protocol_path: Path,
    canonical_manifest_path: Path,
    base_screening_protocol_path: Path,
    data_root: Path,
) -> dict[str, Any]:
    conformer_gate = audit_partial35_package(
        conformer_package_root,
        partial_contract_path,
        generation_protocol_path,
        data_root,
    )
    if conformer_gate.get("status") != "PASS":
        raise ValueError("PARTIAL35 conformer Gate B is not PASS")
    partial_manifest_path = conformer_package_root / "partial_package_manifest.json"
    partial_manifest_hash = sha256_file(partial_manifest_path)
    manifest, candidates = load_candidate_manifest(canonical_manifest_path)

    inventory_path = ligand_package_root / DERIVED_INVENTORY_FILENAME
    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    schema = json.loads(
        (PROJECT_ROOT / "specs" / "top40_ligand_inventory_v05.schema.json").read_text(
            encoding="utf-8"
        )
    )
    validate_json_schema(inventory, schema)
    derived_protocol_path = ligand_package_root / DERIVED_PROTOCOL_FILENAME
    derived_protocol = json.loads(derived_protocol_path.read_text(encoding="utf-8"))
    validate_screening_protocol(derived_protocol)
    if derived_protocol["contracts"]["ligand_inventory"] != {
        "path": DERIVED_INVENTORY_FILENAME,
        "path_base": "protocol_directory",
        "sha256": sha256_file(inventory_path),
    }:
        raise ValueError("Derived screening protocol ligand inventory contract mismatch")
    if inventory.get("candidate_manifest_sha256") != sha256_file(canonical_manifest_path):
        raise ValueError("PARTIAL35 inventory canonical manifest SHA-256 mismatch")
    contract_result = validate_partial35_inventory_contract(
        inventory,
        candidates,
        derived_protocol,
        expected_partial_package_manifest_sha256=partial_manifest_hash,
    )
    inventory_partial_manifest = _resolve_under_data_root(
        data_root,
        inventory["partial_package"]["manifest_path"],
        "partial conformer package manifest",
    )
    if (
        inventory_partial_manifest != partial_manifest_path.resolve()
        or sha256_file(inventory_partial_manifest) != partial_manifest_hash
    ):
        raise ValueError("PARTIAL35 inventory package manifest path/hash mismatch")

    candidate_by_id = {row["candidate_id"]: row for row in candidates}
    core_manifest = json.loads(
        (conformer_package_root / "preparation_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    core_sources = {
        (candidate["candidate_id"], formal["formal_conformer_ordinal"]): {
            "path": (conformer_package_root / formal["path"]).resolve(),
            "sha256": formal["sha256"],
        }
        for candidate in core_manifest["candidates"]
        for formal in candidate["formal_conformers"]
    }
    base_protocol = json.loads(base_screening_protocol_path.read_text(encoding="utf-8"))
    preparation_hash = ligand_preparation_contract_sha256(base_protocol)
    audited = 0
    for slot in inventory["slots"]:
        if slot["readiness_status"] != "approved_for_exploratory_screening":
            continue
        candidate_id = slot["candidate_id"]
        conformer = slot["conformer_name"]
        run_id = f"{candidate_id}_{conformer}"
        source = _resolve_under_data_root(data_root, slot["source_sdf_path"], "source SDF")
        pdbqt = _resolve_under_data_root(data_root, slot["ligand_pdbqt_path"], "ligand PDBQT")
        audit_path = _resolve_under_data_root(data_root, slot["pdbqt_audit_path"], "PDBQT audit")
        for path, expected, label in (
            (source, slot["source_sdf_sha256"], "source SDF"),
            (pdbqt, slot["ligand_pdbqt_sha256"], "ligand PDBQT"),
            (audit_path, slot["pdbqt_audit_sha256"], "PDBQT audit"),
        ):
            actual = sha256_file(path)
            if actual != expected:
                raise ValueError(f"{run_id} {label} SHA-256 mismatch")
        expected_source = core_sources.get((candidate_id, conformer))
        if expected_source != {"path": source, "sha256": slot["source_sdf_sha256"]}:
            raise ValueError(f"{run_id} source SDF is not the audited PARTIAL35 conformer")
        candidate = candidate_by_id[candidate_id]
        source_audit = audit_source_sdf(
            source,
            candidate_id,
            candidate["sequence"],
            _expected_formal_charge(candidate["sequence"]),
        )
        pdbqt_audit = audit_ligand_pdbqt(pdbqt, source, candidate["sequence"])
        coordinate_audit = compare_pdbqt_coordinates(source, pdbqt)
        payload = json.loads(audit_path.read_text(encoding="utf-8"))
        if payload.get("status") != "PASS" or payload.get("run_id") != run_id:
            raise ValueError(f"{run_id} PDBQT audit identity/status mismatch")
        if payload.get("source_sdf", {}).get("sha256") != source_audit["sha256"]:
            raise ValueError(f"{run_id} PDBQT audit source SDF provenance mismatch")
        if payload.get("ligand_pdbqt", {}).get("sha256") != pdbqt_audit["sha256"]:
            raise ValueError(f"{run_id} PDBQT audit artifact provenance mismatch")
        if payload.get("ligand_pdbqt", {}).get("coordinate_identity") != coordinate_audit:
            raise ValueError(f"{run_id} PDBQT audit coordinate provenance mismatch")
        generation = payload.get("generation_provenance", {})
        if generation != slot.get("generation_provenance") or generation != {
            "protocol_id": core_manifest["protocol_id"],
            "protocol_sha256": core_manifest["generation_protocol_sha256"],
            "partial_package_id": "competition_v1_partial35_mmff_qc_pass",
            "partial_package_manifest_sha256": partial_manifest_hash,
        }:
            raise ValueError(f"{run_id} partial package provenance mismatch")
        preparation = payload.get("ligand_preparation", {})
        if (
            preparation.get("protocol") != base_protocol["ligand_preparation"]
            or preparation.get("contract_sha256") != preparation_hash
            or preparation.get("meeko_version") != "0.8.0"
        ):
            raise ValueError(f"{run_id} ligand preparation provenance mismatch")
        audited += 1
    if audited != 105:
        raise ValueError(f"PARTIAL35 PDBQT audit count must equal 105; actual {audited}")

    hash_inventory_path = ligand_package_root / "PDBQT_SHA256SUMS.tsv"
    with hash_inventory_path.open("r", encoding="utf-8", newline="") as handle:
        hash_rows = list(csv.DictReader(handle, delimiter="\t"))
    expected_hash_rows = [
        {
            "candidate_id": slot["candidate_id"],
            "conformer_name": slot["conformer_name"],
            "path": slot["ligand_pdbqt_path"],
            "sha256": slot["ligand_pdbqt_sha256"],
        }
        for slot in inventory["slots"]
        if slot["readiness_status"] == "approved_for_exploratory_screening"
    ]
    if hash_rows != expected_hash_rows:
        raise ValueError("PARTIAL35 PDBQT hash inventory differs from ligand inventory")

    package_manifest_path = ligand_package_root / PREPARATION_MANIFEST_FILENAME
    package_manifest = json.loads(package_manifest_path.read_text(encoding="utf-8"))
    expected_manifest = {
        "package_classification": "PARTIAL_EXPLORATORY_SCREENING",
        "completeness_statement": "NOT_A_COMPLETE_TOP40_SCREEN",
        "partial_conformer_package_manifest_sha256": partial_manifest_hash,
        "candidate_manifest_id": manifest["manifest_id"],
        "candidate_manifest_sha256": sha256_file(canonical_manifest_path),
        "ligand_preparation_contract_sha256": preparation_hash,
        "meeko_version": "0.8.0",
        "candidate_count": 35,
        "pdbqt_count": 105,
        "inventory_sha256": sha256_file(inventory_path),
        "derived_screening_protocol_sha256": sha256_file(derived_protocol_path),
        "pdbqt_hash_inventory_path": hash_inventory_path.name,
        "pdbqt_hash_inventory_sha256": sha256_file(hash_inventory_path),
        "status": "PASS",
    }
    mismatches = {
        field: {"expected": value, "actual": package_manifest.get(field)}
        for field, value in expected_manifest.items()
        if package_manifest.get(field) != value
    }
    if mismatches:
        raise ValueError(f"PARTIAL35 ligand package manifest mismatch: {mismatches}")
    return {
        **contract_result,
        "candidate_count": 35,
        "pdbqt_count": audited,
        "inventory_sha256": sha256_file(inventory_path),
        "derived_screening_protocol_sha256": sha256_file(derived_protocol_path),
        "ligand_package_manifest_sha256": sha256_file(package_manifest_path),
        "status": "PASS",
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ligand-package", type=Path, required=True)
    parser.add_argument("--conformer-package", type=Path, required=True)
    parser.add_argument("--partial-contract", type=Path, required=True)
    parser.add_argument("--generation-protocol", type=Path, required=True)
    parser.add_argument("--candidate-manifest", type=Path, required=True)
    parser.add_argument("--screening-protocol", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = audit_partial35_ligand_package(
            ligand_package_root=args.ligand_package.resolve(),
            conformer_package_root=args.conformer_package.resolve(),
            partial_contract_path=args.partial_contract.resolve(),
            generation_protocol_path=args.generation_protocol.resolve(),
            canonical_manifest_path=args.candidate_manifest.resolve(),
            base_screening_protocol_path=args.screening_protocol.resolve(),
            data_root=args.data_root.resolve(),
        )
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print(
        "PARTIAL35 PDBQT Gate C PASS: "
        f"candidates={result['candidate_count']} PDBQTs={result['pdbqt_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
