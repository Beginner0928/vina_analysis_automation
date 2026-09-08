"""Prepare one standardized peptide conformer for an audited Vina run."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

from standardized_vina_inputs import (
    audit_ligand_pdbqt,
    audit_source_sdf,
    compare_pdbqt_coordinates,
    comparison_protocol_id,
    prepare_backbone_rigid_ligand,
    sha256,
    write_vina_config,
)
from receptor_registry import resolve_receptor_bundle, validate_spec_against_bundle


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def resolve_read_only(data_root: Path, raw_path: str, label: str) -> Path:
    candidate = Path(raw_path)
    if candidate.is_absolute():
        raise ValueError(f"{label} must be relative to --data-root")
    resolved = (data_root.resolve() / candidate).resolve()
    if not resolved.is_relative_to(data_root.resolve()):
        raise ValueError(f"{label} resolves outside --data-root")
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} is missing: {resolved}")
    return resolved


def verify_hash(path: Path, expected: str, label: str) -> str:
    actual = sha256(path)
    if actual.lower() != expected.lower():
        raise ValueError(f"SHA-256 mismatch for {label}: expected {expected}, actual {actual}")
    return actual


def copy_verified_once(source: Path, destination: Path, expected_hash: str) -> str:
    source_hash = verify_hash(source, expected_hash, str(source))
    if destination.exists():
        destination_hash = sha256(destination)
        if destination_hash != source_hash:
            raise FileExistsError(
                f"Existing destination differs; refusing overwrite: {destination}"
            )
        return destination_hash
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as reader, destination.open("xb") as writer:
        for chunk in iter(lambda: reader.read(1024 * 1024), b""):
            writer.write(chunk)
    return verify_hash(destination, source_hash, str(destination))


def write_json_once(path: Path, value: dict[str, Any], allow_identical: bool = False) -> None:
    text = json.dumps(value, indent=2, ensure_ascii=False) + "\n"
    if path.exists():
        if allow_identical and path.read_text(encoding="utf-8") == text:
            return
        raise FileExistsError(f"Refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def prepare_one(
    spec_path: Path, data_root: Path, output_root: Path, conformer_name: str
) -> Path:
    spec = load_json(spec_path)
    if spec.get("schema_version") != "0.4":
        raise ValueError("Standardized docking spec schema_version must be 0.4")
    conformers = {
        item["conformer_name"]: item for item in spec.get("conformers", [])
    }
    if conformer_name not in conformers:
        raise ValueError(f"Unknown conformer {conformer_name!r}")
    conformer = conformers[conformer_name]
    run_id = f"{spec['candidate']}_{conformer_name}"

    source_sdf = resolve_read_only(
        data_root, conformer["source_sdf_path"], f"{conformer_name}.source_sdf_path"
    )
    source_audit = audit_source_sdf(
        source_sdf,
        spec["candidate"],
        spec["sequence"],
        int(spec["formal_chemistry"]["formal_charge"]),
    )
    verify_hash(
        source_sdf, conformer["source_sdf_sha256"], f"{conformer_name}.source_sdf"
    )

    registry_path = resolve_read_only(
        data_root, spec["receptor_registry_path"], "receptor_registry_path"
    )
    registry_hash = verify_hash(
        registry_path,
        spec["receptor_registry_sha256"],
        "receptor_registry_path",
    )
    receptor_bundle = resolve_receptor_bundle(registry_path, data_root, spec["group"])
    registry_identity = validate_spec_against_bundle(spec, receptor_bundle)
    receptor_items = (
        (
            "receptor_monomer_pdb",
            f"receptors/{spec['receptor_id']}.pdb",
        ),
        (
            "receptor_pdbqt",
            f"receptors/{spec['receptor_id']}.pdbqt",
        ),
        (
            "receptor_dimer_pdb",
            f"receptors/{spec['template_id']}_TfR1_AB.pdb",
        ),
    )
    receptor_audits: dict[str, Any] = {}
    for bundle_field, destination_name in receptor_items:
        source = receptor_bundle[bundle_field]
        destination = output_root / destination_name
        expected_hash = receptor_bundle["sha256"][bundle_field]
        copied_hash = copy_verified_once(source, destination, expected_hash)
        receptor_audits[bundle_field] = {
            "source_path": receptor_bundle["source_paths"][bundle_field],
            "source_sha256": expected_hash,
            "copied_path": destination.relative_to(output_root).as_posix(),
            "copied_sha256": copied_hash,
            "hash_verification": "MATCH",
        }

    box = {
        "template": receptor_bundle["template_id"],
        "center_A": tuple(receptor_bundle["vina_box"]["center"]),
        "size_A": tuple(receptor_bundle["vina_box"]["size"]),
    }

    protocol_basis = {
        "ligand_preparation": spec["ligand_preparation"],
        "group": spec["group"],
        "template_id": spec["template_id"],
        "receptor_pdbqt_sha256": receptor_bundle["sha256"]["receptor_pdbqt"],
        "box_center_A": box["center_A"],
        "box_size_A": box["size_A"],
        "docking_protocol": spec["docking_protocol"],
    }
    protocol = {
        "comparison_protocol_id": comparison_protocol_id(protocol_basis),
        "screening_protocol_status": spec["screening_protocol_status"],
        **protocol_basis,
        "receptor_registry_id": receptor_bundle["registry_id"],
        "receptor_registry_path": spec["receptor_registry_path"],
        "receptor_registry_sha256": registry_hash,
    }
    write_json_once(output_root / "protocol.json", protocol, allow_identical=True)

    copied_sdf = output_root / "source_sdf" / f"{run_id}.sdf"
    copy_verified_once(source_sdf, copied_sdf, conformer["source_sdf_sha256"])
    ligand_pdbqt = output_root / "ligands" / f"{run_id}.pdbqt"
    prep_audit = prepare_backbone_rigid_ligand(copied_sdf, ligand_pdbqt)
    pdbqt_audit = audit_ligand_pdbqt(ligand_pdbqt, copied_sdf, spec["sequence"])
    coordinate_audit = compare_pdbqt_coordinates(copied_sdf, ligand_pdbqt)

    config_path = output_root / "configs" / f"{run_id}.txt"
    write_vina_config(
        config_path,
        f"receptors/{spec['receptor_id']}.pdbqt",
        f"ligands/{run_id}.pdbqt",
        box,
        spec["docking_protocol"],
    )
    audit = {
        "audit_datetime": datetime.now().astimezone().isoformat(timespec="seconds"),
        "run_id": run_id,
        "candidate": spec["candidate"],
        "sequence": spec["sequence"],
        "peptide_length": spec["peptide_length"],
        "conformer": conformer_name,
        "screening_protocol_status": spec["screening_protocol_status"],
        "comparison_protocol_id": protocol["comparison_protocol_id"],
        "source_sdf": source_audit,
        "source_sdf_hash_verification": "MATCH",
        "generated_pdbqt": {
            **prep_audit,
            **pdbqt_audit,
            "coordinate_identity": coordinate_audit,
        },
        "formal_chemistry": spec["formal_chemistry"],
        "receptor_registry": {
            **registry_identity,
            "spec_path": spec["receptor_registry_path"],
            "resolved_path": str(registry_path),
            "sha256": registry_hash,
            "hash_verification": "MATCH",
            "coordinate_audit": receptor_bundle["coordinate_audit"],
        },
        "receptors": receptor_audits,
        "box": box,
        "config": {
            "path": config_path.relative_to(output_root).as_posix(),
            "sha256": sha256(config_path),
        },
        "software": {
            "python_version": platform.python_version(),
            "python_executable": str(Path(sys.executable).resolve()),
            "meeko": spec["ligand_preparation"]["software"],
        },
    }
    audit_path = output_root / "audit" / f"{run_id}_input_audit.json"
    write_json_once(audit_path, audit)
    return audit_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--conformer", required=True)
    args = parser.parse_args(argv)
    try:
        audit_path = prepare_one(
            args.spec.resolve(),
            args.data_root.resolve(),
            args.output_root.resolve(),
            args.conformer,
        )
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print(f"Prepared and audited: {audit_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
