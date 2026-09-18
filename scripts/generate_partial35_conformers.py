"""Generate the frozen PARTIAL35 competition_v1 conformer package (no docking)."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import rdkit

from batch_preflight_v05 import inspect_git_release
from competition_conformer_core import (
    generate_candidate_conformers,
    load_authoritative_generation_inputs,
    sha256_file,
    write_generation_package,
)
from partial35_contract_v05 import (
    canonical_json_bytes,
    partial_package_directory,
    load_partial35_contract,
    validate_partial35_contract,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PARTIAL_MANIFEST_FILENAME = "partial_package_manifest.json"


def _write_json_exclusive(path: Path, value: dict[str, Any]) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def _result_provenance(result: dict[str, Any]) -> dict[str, Any]:
    selection = result["selection"]
    return {
        "candidate_id": result["candidate"]["candidate_id"],
        "generated_source_conformer_count": result["generated_source_conformer_count"],
        "mmff_converged_conformer_count": result["mmff_converged_conformer_count"],
        "selected_source_conformer_ids": selection["selected_source_conformer_ids"],
        "selected_mmff94s_energies": selection["selected_energies"],
        "selection_rmsd_threshold_used_A": selection["rmsd_threshold_used_A"],
        "pairwise_rmsd_A": selection["pairwise_rmsd_A"],
        "elapsed_seconds": result.get("generation_elapsed_seconds"),
        "status": "PASS",
    }


def generate_partial35_package(
    *,
    contract: dict[str, Any],
    protocol: dict[str, Any],
    protocol_sha256: str,
    canonical_candidates: list[dict[str, Any]],
    canonical_manifest_sha256: str,
    chemistry_contract: dict[str, Any],
    chemistry_contract_sha256: str,
    package_directory: Path,
    generator_git_commit: str,
    generate_fn: Callable[[dict[str, Any], dict[str, Any]], dict[str, Any]] = generate_candidate_conformers,
    write_fn: Callable[..., Path] = write_generation_package,
) -> Path:
    """Generate sequentially and publish an immutable isolated PARTIAL35 package."""
    if package_directory.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing partial package directory: {package_directory}"
        )
    if re.fullmatch(r"[0-9a-f]{40}", generator_git_commit) is None:
        raise ValueError("generator_git_commit must be a lowercase 40-character Git SHA")
    selected = validate_partial35_contract(
        contract,
        canonical_candidates,
        protocol,
        protocol_sha256,
        canonical_manifest_sha256=canonical_manifest_sha256,
    )
    if chemistry_contract.get("contract_id") != contract["chemistry_contract"]["contract_id"]:
        raise ValueError("Loaded chemistry contract identity differs from PARTIAL35 contract")
    if chemistry_contract_sha256 != contract["chemistry_contract"]["sha256"]:
        raise ValueError("Loaded chemistry contract SHA-256 differs from PARTIAL35 contract")

    results: list[dict[str, Any]] = []
    started = time.perf_counter()
    for candidate in selected:
        candidate_started = time.perf_counter()
        result = generate_fn(candidate, protocol)
        result["generation_elapsed_seconds"] = time.perf_counter() - candidate_started
        results.append(result)

    core_manifest = write_fn(
        results,
        protocol,
        protocol_sha256,
        chemistry_contract,
        package_directory,
        candidate_manifest_sha256=canonical_manifest_sha256,
        chemistry_contract_sha256=chemistry_contract_sha256,
    )
    core_manifest_hash = sha256_file(core_manifest)
    identity = {
        "package_id": contract["package_id"],
        "package_classification": contract["package_classification"],
        "completeness_statement": contract["completeness_statement"],
        "parent_generation_protocol": contract["parent_generation_protocol"],
        "canonical_candidate_manifest": contract["canonical_candidate_manifest"],
        "chemistry_contract": contract["chemistry_contract"],
        "included_candidate_ids": contract["included_candidate_ids"],
        "included_candidate_set_sha256": contract["included_candidate_set_sha256"],
        "excluded_candidates": contract["excluded_candidates"],
        "generator_git_commit": generator_git_commit,
        "core_preparation_manifest_sha256": core_manifest_hash,
    }
    payload = {
        "schema_version": "0.5-partial35-package",
        **identity,
        "package_identity_sha256": hashlib.sha256(canonical_json_bytes(identity)).hexdigest(),
        "included_candidate_count": len(selected),
        "expected_formal_sdf_count": len(selected)
        * int(protocol["selection"]["final_conformer_count"]),
        "generation_protocol_id": protocol["protocol_id"],
        "generation_protocol_sha256": protocol_sha256,
        "core_preparation_manifest_path": core_manifest.relative_to(
            package_directory
        ).as_posix(),
        "candidate_generation": [_result_provenance(result) for result in results],
        "software": {
            "python_version": platform.python_version(),
            "python_executable": str(Path(sys.executable).resolve()),
            "rdkit_version": rdkit.__version__,
        },
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "total_generation_elapsed_seconds": time.perf_counter() - started,
        "status": "PASS",
    }
    output = package_directory / PARTIAL_MANIFEST_FILENAME
    _write_json_exclusive(output, payload)
    return output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        data_root = args.data_root.resolve()
        contract = load_partial35_contract(args.contract.resolve())
        (
            protocol,
            protocol_sha256,
            candidates,
            manifest_sha256,
            chemistry_contract,
            chemistry_sha256,
        ) = load_authoritative_generation_inputs(args.protocol.resolve(), data_root)
        git = inspect_git_release(PROJECT_ROOT)
        output = generate_partial35_package(
            contract=contract,
            protocol=protocol,
            protocol_sha256=protocol_sha256,
            canonical_candidates=candidates,
            canonical_manifest_sha256=manifest_sha256,
            chemistry_contract=chemistry_contract,
            chemistry_contract_sha256=chemistry_sha256,
            package_directory=partial_package_directory(data_root, contract, protocol),
            generator_git_commit=git["head_commit"],
        )
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print(f"PARTIAL35 conformer generation PASS: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
