"""Perform the pre-docking gate audit for one prepared standardized Vina input."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

from standardized_vina_inputs import (
    audit_ligand_pdbqt,
    audit_source_sdf,
    compare_pdbqt_coordinates,
    comparison_protocol_id,
    sha256,
)
from receptor_registry import resolve_receptor_bundle, validate_spec_against_bundle


def read_config(path: Path) -> dict[str, str]:
    fields: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        key, separator, value = line.partition("=")
        if not separator or not key.strip() or key.strip() in fields:
            raise ValueError(f"Malformed or duplicate Vina config line: {raw_line!r}")
        fields[key.strip()] = value.strip()
    return fields


def audit_one(spec_path: Path, data_root: Path, root: Path, conformer: str) -> Path:
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    selected = [row for row in spec["conformers"] if row["conformer_name"] == conformer]
    if len(selected) != 1:
        raise ValueError(f"Conformer not found exactly once: {conformer}")
    conf_spec = selected[0]
    run_id = f"{spec['candidate']}_{conformer}"
    source_sdf = root / "source_sdf" / f"{run_id}.sdf"
    ligand_pdbqt = root / "ligands" / f"{run_id}.pdbqt"
    config_path = root / "configs" / f"{run_id}.txt"
    protocol_path = root / "protocol.json"
    for path in (source_sdf, ligand_pdbqt, config_path, protocol_path):
        if not path.is_file():
            raise FileNotFoundError(f"Prepared input is missing: {path}")

    source_audit = audit_source_sdf(
        source_sdf,
        spec["candidate"],
        spec["sequence"],
        int(spec["formal_chemistry"]["formal_charge"]),
    )
    if source_audit["sha256"] != conf_spec["source_sdf_sha256"]:
        raise ValueError("Copied standardized SDF hash differs from approved source hash")
    pdbqt_audit = audit_ligand_pdbqt(ligand_pdbqt, source_sdf, spec["sequence"])
    coordinate_audit = compare_pdbqt_coordinates(source_sdf, ligand_pdbqt)

    registry_raw = Path(spec["receptor_registry_path"])
    if registry_raw.is_absolute():
        raise ValueError("receptor_registry_path must be relative to --data-root")
    resolved_data_root = data_root.resolve()
    registry = (resolved_data_root / registry_raw).resolve()
    if not registry.is_relative_to(resolved_data_root):
        raise ValueError("receptor_registry_path resolves outside --data-root")
    if not registry.is_file():
        raise FileNotFoundError(f"Receptor registry is missing: {registry}")
    if sha256(registry).lower() != spec["receptor_registry_sha256"].lower():
        raise ValueError("Receptor registry hash mismatch")
    bundle = resolve_receptor_bundle(registry, data_root, spec["group"])
    registry_identity = validate_spec_against_bundle(spec, bundle)
    box = {
        "template": bundle["template_id"],
        "center_A": tuple(bundle["vina_box"]["center"]),
        "size_A": tuple(bundle["vina_box"]["size"]),
    }
    config = read_config(config_path)
    expected: dict[str, str] = {
        "receptor": f"receptors/{spec['receptor_id']}.pdbqt",
        "ligand": f"ligands/{run_id}.pdbqt",
        **{
            f"center_{axis}": str(value)
            for axis, value in zip("xyz", box["center_A"])
        },
        **{f"size_{axis}": str(value) for axis, value in zip("xyz", box["size_A"])},
        **{
            field: str(spec["docking_protocol"][field])
            for field in ("exhaustiveness", "num_modes", "energy_range", "seed", "cpu")
        },
    }
    if config != expected:
        raise ValueError(f"Vina config differs from audited sources: expected {expected}, actual {config}")

    receptor_hashes: dict[str, Any] = {}
    for suffix, expected_hash in (
        (f"{spec['receptor_id']}.pdb", bundle["sha256"]["receptor_monomer_pdb"]),
        (f"{spec['receptor_id']}.pdbqt", bundle["sha256"]["receptor_pdbqt"]),
        (f"{spec['template_id']}_TfR1_AB.pdb", bundle["sha256"]["receptor_dimer_pdb"]),
    ):
        path = root / "receptors" / suffix
        actual = sha256(path)
        if actual != expected_hash:
            raise ValueError(f"Prepared receptor hash mismatch: {path}")
        receptor_hashes[suffix] = {"sha256": actual, "status": "MATCH"}

    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    protocol_basis = {
        "ligand_preparation": spec["ligand_preparation"],
        "group": bundle["group"],
        "template_id": bundle["template_id"],
        "receptor_pdbqt_sha256": bundle["sha256"]["receptor_pdbqt"],
        "box_center_A": box["center_A"],
        "box_size_A": box["size_A"],
        "docking_protocol": spec["docking_protocol"],
    }
    expected_protocol_id = comparison_protocol_id(protocol_basis)
    if protocol.get("comparison_protocol_id") != expected_protocol_id:
        raise ValueError(
            "comparison_protocol_id differs from registry-resolved protocol identity"
        )
    result = {
        "audit_datetime": datetime.now().astimezone().isoformat(timespec="seconds"),
        "gate": "PRE_DOCKING_INPUT_AUDIT",
        "gate_status": "PASS",
        "run_id": run_id,
        "screening_protocol_status": spec["screening_protocol_status"],
        "comparison_protocol_id": protocol["comparison_protocol_id"],
        "source_sdf": source_audit,
        "source_sdf_hash_verification": "MATCH",
        "generated_pdbqt": {
            **pdbqt_audit,
            "coordinate_identity": coordinate_audit,
        },
        "formal_chemistry": spec["formal_chemistry"],
        "receptor_registry": {
            **registry_identity,
            "spec_path": spec["receptor_registry_path"],
            "resolved_path": str(registry),
            "sha256": sha256(registry),
            "hash_verification": "MATCH",
            "coordinate_audit": bundle["coordinate_audit"],
        },
        "receptor_files": receptor_hashes,
        "box_registry": box,
        "vina_config": {"fields": config, "sha256": sha256(config_path), "status": "MATCH"},
        "protocol_json": {"sha256": sha256(protocol_path), "content": protocol},
    }
    output = root / "audit" / f"{run_id}_pre_docking_gate.json"
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing gate audit: {output}")
    with output.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(result, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--conformer", required=True)
    args = parser.parse_args(argv)
    try:
        path = audit_one(
            args.spec.resolve(), args.data_root.resolve(), args.root.resolve(), args.conformer
        )
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print(f"Pre-docking gate PASS: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
