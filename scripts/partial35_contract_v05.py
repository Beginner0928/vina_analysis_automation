"""Frozen identity and validation helpers for the PARTIAL35 exploratory package."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from candidate_manifest_v05 import validate_json_schema


PACKAGE_ID = "competition_v1_partial35_mmff_qc_pass"
PACKAGE_CLASSIFICATION = "PARTIAL_EXPLORATORY_SCREENING"
COMPLETENESS_STATEMENT = "NOT_A_COMPLETE_TOP40_SCREEN"
PARENT_PROTOCOL_ID = "tfr1_top40_sdfgen_competition_v1_b3derived"
PARENT_PROTOCOL_SHA256 = (
    "f79031dda94e148b851774b33082117ffde70f12b4cbf70b7f90a3266226a4a7"
)
PARTIAL35_CANDIDATE_IDS = (
    "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A9", "A10",
    "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10",
    "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
    "D1", "D2", "D4", "D5", "D6", "D8", "D9",
)
EXCLUDED_CANDIDATES = (
    ("A8", "CONFORMER_GENERATION_QC_FAIL"),
    ("B1", "CONFORMER_GENERATION_QC_FAIL"),
    ("D3", "CONFORMER_GENERATION_QC_FAIL"),
    ("D7", "CONFORMER_GENERATION_QC_FAIL"),
    ("D10", "CONFORMER_GENERATION_QC_FAIL"),
)
OUTPUT_RELATIVE = Path("standardized_ligands") / "partial_packages" / PACKAGE_ID


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def canonical_candidate_set_sha256(candidates: list[dict[str, Any]]) -> str:
    identity = [
        {
            field: row[field]
            for field in ("candidate_id", "group", "sequence", "peptide_length")
        }
        for row in candidates
    ]
    return hashlib.sha256(canonical_json_bytes(identity)).hexdigest()


def load_partial35_contract(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        schema = json.loads(
            path.with_name("partial35_exploratory_package_v05.schema.json").read_text(
                encoding="utf-8"
            )
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot load PARTIAL35 contract: {exc}") from exc
    if not isinstance(value, dict) or not isinstance(schema, dict):
        raise ValueError("PARTIAL35 contract and schema roots must be objects")
    validate_json_schema(value, schema)
    return value


def validate_partial35_contract(
    contract: dict[str, Any],
    canonical_candidates: list[dict[str, Any]],
    protocol: dict[str, Any],
    protocol_sha256: str,
    *,
    canonical_manifest_sha256: str,
) -> list[dict[str, Any]]:
    expected_scalars = {
        "package_id": PACKAGE_ID,
        "package_classification": PACKAGE_CLASSIFICATION,
        "completeness_statement": COMPLETENESS_STATEMENT,
        "included_candidate_count": len(PARTIAL35_CANDIDATE_IDS),
        "expected_formal_sdf_count": len(PARTIAL35_CANDIDATE_IDS) * 3,
    }
    mismatches = {
        field: {"expected": expected, "actual": contract.get(field)}
        for field, expected in expected_scalars.items()
        if contract.get(field) != expected
    }
    if mismatches:
        raise ValueError(f"PARTIAL35 identity/completeness mismatch: {mismatches}")

    included = contract.get("included_candidate_ids")
    if not isinstance(included, list):
        raise ValueError("included_candidate_ids must be an array")
    if len(included) != len(set(included)):
        raise ValueError("included_candidate_ids contains a duplicate")
    known = {row["candidate_id"] for row in canonical_candidates}
    unknown = sorted(set(included) - known)
    if unknown:
        raise ValueError(f"included_candidate_ids contains unknown candidates: {unknown}")
    if included != list(PARTIAL35_CANDIDATE_IDS):
        raise ValueError("included_candidate_ids differs from the frozen canonical PARTIAL35 set")

    excluded = contract.get("excluded_candidates")
    expected_excluded = [
        {
            "candidate_id": candidate_id,
            "generation_exclusion_reason": reason,
            "screening_status": "CONFORMER_GENERATION_QC_FAIL_NOT_DOCKED",
        }
        for candidate_id, reason in EXCLUDED_CANDIDATES
    ]
    if excluded != expected_excluded:
        raise ValueError("excluded_candidates differs from the frozen five-candidate QC set")

    parent = contract.get("parent_generation_protocol", {})
    if protocol.get("protocol_id") != PARENT_PROTOCOL_ID:
        raise ValueError("Loaded generation protocol ID is not competition_v1")
    if protocol_sha256 != PARENT_PROTOCOL_SHA256:
        raise ValueError("Loaded generation protocol SHA-256 differs from the frozen value")
    if parent != {"protocol_id": PARENT_PROTOCOL_ID, "sha256": PARENT_PROTOCOL_SHA256}:
        raise ValueError("PARTIAL35 parent generation protocol contract mismatch")

    authority = contract.get("canonical_candidate_manifest", {})
    if authority.get("manifest_id") != "tfr1_top40_v05":
        raise ValueError("PARTIAL35 canonical manifest ID mismatch")
    if authority.get("sha256") != canonical_manifest_sha256:
        raise ValueError("PARTIAL35 canonical manifest SHA-256 mismatch")
    chemistry = contract.get("chemistry_contract", {})
    expected_chemistry = protocol.get("chemistry_contract")
    if chemistry != {
        "contract_id": expected_chemistry.get("contract_id"),
        "sha256": expected_chemistry.get("sha256"),
    }:
        raise ValueError("PARTIAL35 chemistry contract identity mismatch")

    selected = [row for row in canonical_candidates if row["candidate_id"] in set(included)]
    if [row["candidate_id"] for row in selected] != included:
        raise ValueError("PARTIAL35 selection is not in canonical manifest order")
    digest = canonical_candidate_set_sha256(selected)
    if contract.get("included_candidate_set_sha256") != digest:
        raise ValueError("PARTIAL35 included candidate-set SHA-256 mismatch")
    if contract.get("included_candidate_set_sha256_definition") != (
        "sha256_canonical_json_identity_rows_candidate_id_group_sequence_"
        "peptide_length_in_manifest_order"
    ):
        raise ValueError("PARTIAL35 candidate-set digest definition mismatch")
    return selected


def partial_package_directory(
    data_root: Path, contract: dict[str, Any], protocol: dict[str, Any]
) -> Path:
    if contract.get("package_id") != PACKAGE_ID:
        raise ValueError("Cannot resolve output for a non-PARTIAL35 package identity")
    output = contract.get("output", {}).get("relative_directory")
    if output != OUTPUT_RELATIVE.as_posix():
        raise ValueError("PARTIAL35 output directory differs from the frozen isolated path")
    resolved = data_root.resolve() / OUTPUT_RELATIVE
    reserved = data_root.resolve() / "standardized_ligands" / protocol["protocol_id"]
    if resolved == reserved:
        raise ValueError("Partial package cannot use the reserved complete-Top40 directory")
    return resolved
