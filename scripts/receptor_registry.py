"""Resolve and preflight locked A/B/C/D TfR1 receptor template bundles."""

from __future__ import annotations

import json
from math import dist
from pathlib import Path
from typing import Any

from pose_analysis_core import CoordinateAtom, parse_receptor_pdb_atoms
from standardized_vina_inputs import read_box_registry, sha256


PRODUCTION_GROUPS = {"A", "B", "C", "D"}
PATH_FIELDS = (
    "receptor_monomer_pdb",
    "receptor_pdbqt",
    "receptor_dimer_pdb",
)
PRODUCTION_PROVENANCE_HASH_FIELDS = {
    "normalized_receptor_monomer_pdb_sha256": "receptor_monomer_pdb",
    "docking_receptor_pdbqt_sha256": "receptor_pdbqt",
    "receptor_dimer_pdb_sha256": "receptor_dimer_pdb",
}


def load_registry(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Receptor registry is missing: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("Receptor registry must be a JSON object")
    if value.get("schema_version") != "0.4":
        raise ValueError("Receptor registry schema_version must be 0.4")
    groups = value.get("groups")
    if not isinstance(groups, dict) or set(groups) != PRODUCTION_GROUPS:
        raise ValueError("Receptor registry groups must be exactly A, B, C, and D")
    rules = value.get("production_rules")
    required_rules = {
        "group_is_required": True,
        "allow_cross_group_coordinate_mixing": False,
        "allow_delete_bad_res": False,
        "allow_bad_res": False,
        "allow_serorder_test_files_as_production": False,
        "fail_on_missing_registry_file": True,
        "fail_on_hash_mismatch": True,
    }
    if not isinstance(rules, dict) or any(
        rules.get(key) != expected for key, expected in required_rules.items()
    ):
        raise ValueError("Receptor registry production rules are incomplete or unsafe")
    return value


def _resolve_file(data_root: Path, raw_path: str, label: str) -> Path:
    if "serorder_test" in raw_path.lower():
        raise ValueError(f"SERorder_test diagnostic cannot be production input: {label}")
    candidate = Path(raw_path)
    if candidate.is_absolute():
        raise ValueError(f"Registry path must be relative to data root: {label}")
    root = data_root.resolve()
    resolved = (root / candidate).resolve()
    if not resolved.is_relative_to(root):
        raise ValueError(f"Registry path resolves outside data root: {label}")
    if not resolved.is_file():
        raise FileNotFoundError(f"Registry production file is missing: {resolved}")
    return resolved


def _parse_receptor_pdbqt(path: Path, chain: str = "A") -> list[CoordinateAtom]:
    atoms: list[CoordinateAtom] = []
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[21:22].strip() != chain:
            continue
        fields = line.split()
        atoms.append(
            CoordinateAtom(
                serial=int(line[6:11]),
                name=line[12:16].strip(),
                atom_type=fields[-1],
                xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
                residue_number=int(line[22:26]),
                residue_name=line[17:20].strip(),
                chain=chain,
            )
        )
    if not atoms:
        raise ValueError(f"Receptor PDBQT chain {chain} contains no atoms: {path}")
    return atoms


def _atom_map(atoms: list[CoordinateAtom]) -> dict[tuple[int, str], CoordinateAtom]:
    result: dict[tuple[int, str], CoordinateAtom] = {}
    for atom in atoms:
        if not atom.is_heavy or atom.residue_number is None:
            continue
        key = (atom.residue_number, atom.name)
        if key in result:
            raise ValueError(f"Duplicate receptor atom key: {key}")
        result[key] = atom
    return result


def _coordinate_delta_with_terminal_oxygen_equivalence(
    left: dict[tuple[int, str], CoordinateAtom],
    right: dict[tuple[int, str], CoordinateAtom],
) -> tuple[float, float, bool]:
    raw_deltas = {key: dist(left[key].xyz, right[key].xyz) for key in left}
    raw_maximum = max(raw_deltas.values())
    terminal_residue = max(key[0] for key in left)
    oxygen_keys = ((terminal_residue, "O"), (terminal_residue, "OXT"))
    if not all(key in left and key in right for key in oxygen_keys):
        return raw_maximum, raw_maximum, False
    direct_sum = sum(raw_deltas[key] for key in oxygen_keys)
    swapped_deltas = (
        dist(left[oxygen_keys[0]].xyz, right[oxygen_keys[1]].xyz),
        dist(left[oxygen_keys[1]].xyz, right[oxygen_keys[0]].xyz),
    )
    swapped = sum(swapped_deltas) < direct_sum
    effective = dict(raw_deltas)
    if swapped:
        effective[oxygen_keys[0]], effective[oxygen_keys[1]] = swapped_deltas
    return max(effective.values()), raw_maximum, swapped


def audit_coordinate_bundle(
    monomer: Path,
    receptor_pdbqt: Path,
    dimer: Path,
    target_residues: list[int],
) -> dict[str, Any]:
    monomer_atoms = parse_receptor_pdb_atoms(monomer, "A")
    pdbqt_atoms = _parse_receptor_pdbqt(receptor_pdbqt, "A")
    dimer_a_atoms = parse_receptor_pdb_atoms(dimer, "A")
    dimer_b_atoms = parse_receptor_pdb_atoms(dimer, "B")
    monomer_map = _atom_map(monomer_atoms)
    pdbqt_map = _atom_map(pdbqt_atoms)
    dimer_a_map = _atom_map(dimer_a_atoms)
    if set(monomer_map) != set(pdbqt_map):
        raise ValueError("Receptor monomer/PDBQT heavy-atom identities differ")
    if set(monomer_map) != set(dimer_a_map):
        raise ValueError("Receptor monomer/dimer chain-A heavy-atom identities differ")
    pdbqt_max, pdbqt_raw_max, terminal_oxygen_swap = (
        _coordinate_delta_with_terminal_oxygen_equivalence(
            monomer_map, pdbqt_map
        )
    )
    dimer_max = max(
        dist(monomer_map[key].xyz, dimer_a_map[key].xyz) for key in monomer_map
    )
    if pdbqt_max >= 0.001:
        raise ValueError(
            f"Receptor PDB/PDBQT coordinate mismatch: max delta {pdbqt_max:.9f} A"
        )
    if dimer_max > 1e-9:
        raise ValueError(
            f"Receptor monomer/dimer coordinate mismatch: max delta {dimer_max:.9f} A"
        )
    present = sorted(
        residue
        for residue in target_residues
        if any(key[0] == residue for key in monomer_map)
    )
    if present != target_residues:
        raise ValueError(
            f"Target residues missing from receptor bundle: expected {target_residues}, actual {present}"
        )
    return {
        "monomer_heavy_atom_count": len(monomer_map),
        "pdbqt_heavy_atom_count": len(pdbqt_map),
        "dimer_chain_A_heavy_atom_count": len(dimer_a_map),
        "dimer_chain_B_heavy_atom_count": len(_atom_map(dimer_b_atoms)),
        "target_residues_present": present,
        "monomer_to_pdbqt_max_coordinate_delta_A": pdbqt_max,
        "monomer_to_pdbqt_raw_named_atom_max_coordinate_delta_A": pdbqt_raw_max,
        "terminal_O_OXT_equivalent_name_swap": terminal_oxygen_swap,
        "monomer_to_dimer_chain_A_max_coordinate_delta_A": dimer_max,
        "coordinate_system_status": "MATCH",
    }


def resolve_receptor_bundle(
    registry_path: Path, data_root: Path, group: str
) -> dict[str, Any]:
    registry = load_registry(registry_path.resolve())
    if group not in PRODUCTION_GROUPS:
        raise ValueError(
            f"Unknown receptor group {group!r}; allowed groups are A, B, C, D"
        )
    entry = registry["groups"][group]
    template_id = entry.get("template_id")
    if not isinstance(template_id, str) or not template_id:
        raise ValueError(f"Registry group {group} has no template_id")
    if entry.get("coordinate_system_id") != template_id:
        raise ValueError(f"Registry group {group} has cross-template coordinate_system_id")
    production_status = entry.get("production_status")
    if production_status not in {"READY", "approved"}:
        raise ValueError(f"Registry group {group} is not approved for production")

    raw_hashes = entry.get("sha256")
    if not isinstance(raw_hashes, dict):
        raise ValueError(f"Registry group {group} has no locked hashes")
    paths: dict[str, Path] = {}
    actual_hashes: dict[str, str] = {}
    for field in PATH_FIELDS:
        raw_path = entry.get(field)
        expected_hash = raw_hashes.get(field)
        if not isinstance(raw_path, str) or not raw_path:
            raise ValueError(f"Registry group {group} missing {field}")
        if template_id not in Path(raw_path).name:
            raise ValueError(
                f"Registry group {group} cross-template component in {field}: {raw_path}"
            )
        if not isinstance(expected_hash, str) or len(expected_hash) != 64:
            raise ValueError(f"Registry group {group} missing locked SHA-256 for {field}")
        path = _resolve_file(data_root, raw_path, f"groups.{group}.{field}")
        actual_hash = sha256(path)
        if actual_hash.lower() != expected_hash.lower():
            raise ValueError(
                f"SHA-256 mismatch for registry group {group} {field}: "
                f"expected {expected_hash}, actual {actual_hash}"
            )
        paths[field] = path
        actual_hashes[field] = actual_hash

    normalized_required = entry.get("normalized_receptor_required")
    if normalized_required is True:
        if not paths["receptor_monomer_pdb"].name.endswith("_normalized.pdb"):
            raise ValueError(f"Group {group} requires normalized production receptor")
        if not entry.get("normalization"):
            raise ValueError(f"Group {group} lacks normalization provenance")
        if entry.get("normalization_status") != "PASS":
            raise ValueError(f"Group {group} normalization_status must be PASS")
        if production_status != "approved":
            raise ValueError(f"Group {group} production_status must be approved")
        if entry.get("receptor_preparation_blocker") != "RESOLVED":
            raise ValueError(
                f"Group {group} receptor preparation blocker must be RESOLVED"
            )
        provenance_hashes = entry.get("production_provenance_sha256")
        if not isinstance(provenance_hashes, dict) or set(provenance_hashes) != set(
            PRODUCTION_PROVENANCE_HASH_FIELDS
        ):
            raise ValueError(
                f"Group {group} production provenance SHA-256 fields are incomplete"
            )
        for provenance_field, path_field in PRODUCTION_PROVENANCE_HASH_FIELDS.items():
            provenance_hash = provenance_hashes.get(provenance_field)
            if (
                not isinstance(provenance_hash, str)
                or len(provenance_hash) != 64
                or provenance_hash.lower() != actual_hashes[path_field].lower()
                or provenance_hash.lower() != raw_hashes[path_field].lower()
            ):
                raise ValueError(
                    f"Group {group} production provenance SHA-256 mismatch for "
                    f"{provenance_field}"
                )
    elif normalized_required is not False:
        raise ValueError(f"Group {group} normalized_receptor_required must be boolean")

    source_spec = registry.get("box_registry_source")
    if not isinstance(source_spec, dict):
        raise ValueError("Registry lacks audited box_registry_source")
    box_source = _resolve_file(data_root, source_spec.get("path", ""), "box_registry_source")
    actual_box_hash = sha256(box_source)
    if actual_box_hash != source_spec.get("sha256"):
        raise ValueError(
            f"SHA-256 mismatch for box registry: expected {source_spec.get('sha256')}, actual {actual_box_hash}"
        )
    source_box = read_box_registry(box_source, template_id)
    registry_box = entry.get("vina_box")
    if not isinstance(registry_box, dict):
        raise ValueError(f"Registry group {group} has no Vina box")
    center = tuple(float(value) for value in registry_box.get("center", []))
    size = tuple(float(value) for value in registry_box.get("size", []))
    if center != source_box["center_A"] or size != source_box["size_A"]:
        raise ValueError(f"Registry group {group} Vina box differs from audited source")

    targets = registry.get("target_residues")
    if not isinstance(targets, list) or not targets:
        raise ValueError("Registry target_residues must be a non-empty list")
    coordinate_audit = audit_coordinate_bundle(
        paths["receptor_monomer_pdb"],
        paths["receptor_pdbqt"],
        paths["receptor_dimer_pdb"],
        targets,
    )
    return {
        "registry_id": registry["registry_id"],
        "registry_path": registry_path.resolve(),
        "registry_sha256": sha256(registry_path.resolve()),
        "group": group,
        "template_id": template_id,
        "coordinate_system_id": entry["coordinate_system_id"],
        **paths,
        "source_paths": {field: entry[field] for field in PATH_FIELDS},
        "sha256": actual_hashes,
        "vina_box": {"center": list(center), "size": list(size)},
        "target_residues": list(targets),
        "normalized_receptor_required": normalized_required,
        "normalization": entry.get("normalization"),
        "normalization_status": entry.get("normalization_status"),
        "production_status": production_status,
        "receptor_preparation_blocker": entry.get("receptor_preparation_blocker"),
        "production_provenance_sha256": entry.get("production_provenance_sha256"),
        "coordinate_audit": coordinate_audit,
        "preflight_status": "PASS",
    }


def validate_spec_against_bundle(
    spec: dict[str, Any], bundle: dict[str, Any]
) -> dict[str, Any]:
    """Reject any spec identity, box, target, or locked hash outside its group bundle."""
    expected_receptor_id = f"{bundle['template_id']}_TfR1_A"
    identity_expectations = {
        "group": bundle["group"],
        "template_id": bundle["template_id"],
        "receptor_id": expected_receptor_id,
    }
    for field, expected in identity_expectations.items():
        if spec.get(field) != expected:
            raise ValueError(
                f"{field} differs from receptor registry group {bundle['group']}: "
                f"expected {expected!r}, actual {spec.get(field)!r}"
            )

    if "target_residues" in spec and spec["target_residues"] != bundle["target_residues"]:
        raise ValueError("target_residues differ from receptor registry")

    protocol = spec.get("docking_protocol")
    if isinstance(protocol, dict):
        for spec_field, bundle_field in (
            ("box_center_A", "center"),
            ("box_size_A", "size"),
        ):
            if spec_field in protocol and list(protocol[spec_field]) != list(
                bundle["vina_box"][bundle_field]
            ):
                raise ValueError(f"docking_protocol.{spec_field} differs from receptor registry")

    hashes = spec.get("input_sha256")
    if isinstance(hashes, dict):
        hash_mapping = {
            "receptor_monomer_path": "receptor_monomer_pdb",
            "receptor_dimer_path": "receptor_dimer_pdb",
            "docking_receptor_pdbqt_path": "receptor_pdbqt",
        }
        for spec_field, bundle_field in hash_mapping.items():
            if spec_field in hashes and hashes[spec_field].lower() != bundle["sha256"][
                bundle_field
            ].lower():
                raise ValueError(f"{spec_field} hash differs from receptor registry")

    return {
        "registry_id": bundle["registry_id"],
        "registry_sha256": bundle["registry_sha256"],
        "group": bundle["group"],
        "template_id": bundle["template_id"],
        "receptor_id": expected_receptor_id,
        "preflight_status": "PASS",
    }
