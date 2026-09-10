"""Audit the isolated PARTIAL35 competition_v1 conformer package (no docking)."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any

from competition_conformer_core import (
    audit_generation_package_against_contract,
    load_authoritative_generation_inputs,
    sha256_file,
)
from generate_partial35_conformers import PARTIAL_MANIFEST_FILENAME
from partial35_contract_v05 import (
    canonical_json_bytes,
    load_partial35_contract,
    partial_package_directory,
    validate_partial35_contract,
)


def audit_partial35_package(
    package_directory: Path,
    contract_path: Path,
    protocol_path: Path,
    data_root: Path,
) -> dict[str, Any]:
    data_root = data_root.resolve()
    contract = load_partial35_contract(contract_path.resolve())
    (
        protocol,
        protocol_sha256,
        candidates,
        manifest_sha256,
        _chemistry,
        chemistry_sha256,
    ) = load_authoritative_generation_inputs(protocol_path.resolve(), data_root)
    selected = validate_partial35_contract(
        contract,
        candidates,
        protocol,
        protocol_sha256,
        canonical_manifest_sha256=manifest_sha256,
    )
    expected_directory = partial_package_directory(data_root, contract, protocol)
    if package_directory.resolve() != expected_directory:
        raise ValueError(
            f"PARTIAL35 package is not at its frozen isolated path: {package_directory}"
        )
    internal = audit_generation_package_against_contract(
        package_directory, protocol_path, data_root
    )
    if internal.get("candidate_count") != 35 or internal.get("sdf_count") != 105:
        raise ValueError(
            "PARTIAL35 core package must contain exactly 35 candidates and 105 SDFs"
        )
    core_manifest_payload = json.loads(
        (package_directory / "preparation_manifest.json").read_text(encoding="utf-8")
    )
    if [row.get("candidate_id") for row in core_manifest_payload.get("candidates", [])] != [
        row["candidate_id"] for row in selected
    ]:
        raise ValueError("PARTIAL35 core package candidate set/order mismatch")
    manifest_path = package_directory / PARTIAL_MANIFEST_FILENAME
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot read PARTIAL35 package manifest: {exc}") from exc
    if not isinstance(manifest, dict):
        raise ValueError("PARTIAL35 package manifest root must be an object")
    expected = {
        "package_id": contract["package_id"],
        "package_classification": contract["package_classification"],
        "completeness_statement": contract["completeness_statement"],
        "parent_generation_protocol": contract["parent_generation_protocol"],
        "canonical_candidate_manifest": contract["canonical_candidate_manifest"],
        "chemistry_contract": contract["chemistry_contract"],
        "included_candidate_ids": contract["included_candidate_ids"],
        "included_candidate_set_sha256": contract["included_candidate_set_sha256"],
        "excluded_candidates": contract["excluded_candidates"],
        "included_candidate_count": 35,
        "expected_formal_sdf_count": 105,
        "generation_protocol_id": protocol["protocol_id"],
        "generation_protocol_sha256": protocol_sha256,
    }
    mismatches = {
        field: {"expected": value, "actual": manifest.get(field)}
        for field, value in expected.items()
        if manifest.get(field) != value
    }
    if mismatches:
        raise ValueError(f"PARTIAL35 package provenance mismatch: {mismatches}")
    if re.fullmatch(r"[0-9a-f]{40}", str(manifest.get("generator_git_commit"))) is None:
        raise ValueError("PARTIAL35 generator Git commit is invalid")
    core_path = package_directory / str(manifest["core_preparation_manifest_path"])
    core_hash = sha256_file(core_path)
    if core_hash != manifest.get("core_preparation_manifest_sha256"):
        raise ValueError("PARTIAL35 core preparation manifest SHA-256 mismatch")
    identity_fields = (
        "package_id",
        "package_classification",
        "completeness_statement",
        "parent_generation_protocol",
        "canonical_candidate_manifest",
        "chemistry_contract",
        "included_candidate_ids",
        "included_candidate_set_sha256",
        "excluded_candidates",
        "generator_git_commit",
        "core_preparation_manifest_sha256",
    )
    identity = {field: manifest[field] for field in identity_fields}
    identity_hash = hashlib.sha256(canonical_json_bytes(identity)).hexdigest()
    if manifest.get("package_identity_sha256") != identity_hash:
        raise ValueError("PARTIAL35 package identity SHA-256 mismatch")
    rows = manifest.get("candidate_generation")
    if not isinstance(rows, list) or [row.get("candidate_id") for row in rows] != [
        row["candidate_id"] for row in selected
    ]:
        raise ValueError("PARTIAL35 per-candidate generation inventory mismatch")
    if any(row.get("status") != "PASS" for row in rows):
        raise ValueError("PARTIAL35 per-candidate generation status is not all PASS")
    if chemistry_sha256 != contract["chemistry_contract"]["sha256"]:
        raise ValueError("PARTIAL35 chemistry SHA-256 mismatch")
    return {
        **internal,
        "package_id": contract["package_id"],
        "package_identity_sha256": identity_hash,
        "partial_package_manifest_sha256": sha256_file(manifest_path),
        "candidate_count": len(selected),
        "sdf_count": 105,
        "excluded_candidate_count": 5,
        "contract_provenance": "MATCH",
        "status": "PASS",
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--package-directory", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = audit_partial35_package(
            args.package_directory.resolve(),
            args.contract.resolve(),
            args.protocol.resolve(),
            args.data_root.resolve(),
        )
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print(
        "PARTIAL35 conformer audit PASS: "
        f"candidates={result['candidate_count']} SDFs={result['sdf_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
