"""Prepare and audit standardized peptide SDF inputs for Vina without regenerating coordinates."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from math import dist
from pathlib import Path
from typing import Any


PEPTIDE_BACKBONE_SMARTS = "[N;X3,X4+]-[C;X4]-[C;X3](=[O;X1])"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_single_sdf(path: Path):
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    molecules = [molecule for molecule in supplier if molecule is not None]
    if len(molecules) != 1:
        raise ValueError(f"Expected exactly one readable molecule in {path}; found {len(molecules)}")
    return molecules[0]


def audit_source_sdf(
    path: Path, candidate: str, sequence: str, expected_formal_charge: int
) -> dict[str, Any]:
    from rdkit import Chem

    molecule = read_single_sdf(path)
    title = molecule.GetProp("_Name") if molecule.HasProp("_Name") else ""
    if title != path.stem or not title.startswith(f"{candidate}_conf"):
        raise ValueError(f"SDF identity mismatch: title {title!r}, file {path.name!r}")
    recorded_sequence = molecule.GetProp("Sequence") if molecule.HasProp("Sequence") else ""
    if recorded_sequence != sequence:
        raise ValueError(
            f"SDF sequence mismatch for {path.name}: expected {sequence}, actual {recorded_sequence!r}"
        )
    computed_charge = sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())
    recorded_charge = (
        int(molecule.GetProp("FormalCharge"))
        if molecule.HasProp("FormalCharge")
        else computed_charge
    )
    if computed_charge != expected_formal_charge or recorded_charge != expected_formal_charge:
        raise ValueError(
            f"SDF formal charge mismatch for {path.name}: expected {expected_formal_charge}, "
            f"computed {computed_charge}, recorded {recorded_charge}"
        )
    template = Chem.MolFromSequence(sequence)
    if template is None:
        raise ValueError(f"Cannot construct peptide template for sequence {sequence!r}")
    source_heavy = [atom for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
    template_heavy = [atom for atom in template.GetAtoms() if atom.GetAtomicNum() > 1]
    if [atom.GetSymbol() for atom in source_heavy] != [
        atom.GetSymbol() for atom in template_heavy
    ]:
        raise ValueError("SDF heavy-atom identity/order differs from sequence template")
    sites: dict[tuple[int, str], Any] = {}
    for source_atom, template_atom in zip(source_heavy, template_heavy):
        info = template_atom.GetPDBResidueInfo()
        if info is not None:
            sites[(info.GetResidueNumber(), info.GetName().strip())] = source_atom
    n_terminal = sites[(1, "N")]
    c_terminal = sites[(len(sequence), "OXT")]
    arg_nh1 = next(
        (
            source_atom
            for source_atom, template_atom in zip(source_heavy, template_heavy)
            if template_atom.GetPDBResidueInfo() is not None
            and template_atom.GetPDBResidueInfo().GetResidueName().strip() == "ARG"
            and template_atom.GetPDBResidueInfo().GetName().strip() == "NH1"
        ),
        None,
    )
    chemistry_checks = {
        "n_terminus_formal_charge": n_terminal.GetFormalCharge(),
        "n_terminus_hydrogen_count": sum(
            neighbor.GetAtomicNum() == 1 for neighbor in n_terminal.GetNeighbors()
        ),
        "c_terminus_oxt_formal_charge": c_terminal.GetFormalCharge(),
        "c_terminus_oxt_hydrogen_count": sum(
            neighbor.GetAtomicNum() == 1 for neighbor in c_terminal.GetNeighbors()
        ),
        "arginine_nh1_formal_charge": (
            arg_nh1.GetFormalCharge() if arg_nh1 is not None else None
        ),
        "sulfur_sulfur_bond_count": sum(
            bond.GetBeginAtom().GetAtomicNum() == 16
            and bond.GetEndAtom().GetAtomicNum() == 16
            for bond in molecule.GetBonds()
        ),
        "heavy_atom_chirality_matches_L_template": all(
            source_atom.GetChiralTag() == template_atom.GetChiralTag()
            for source_atom, template_atom in zip(source_heavy, template_heavy)
        ),
    }
    expected_checks = {
        "n_terminus_formal_charge": 1,
        "n_terminus_hydrogen_count": 3,
        "c_terminus_oxt_formal_charge": -1,
        "c_terminus_oxt_hydrogen_count": 0,
        "sulfur_sulfur_bond_count": 0,
        "heavy_atom_chirality_matches_L_template": True,
    }
    if "R" in sequence:
        expected_checks["arginine_nh1_formal_charge"] = 1
    mismatches = {
        field: {"expected": expected, "actual": chemistry_checks[field]}
        for field, expected in expected_checks.items()
        if chemistry_checks[field] != expected
    }
    if mismatches:
        raise ValueError(f"SDF formal chemistry checks failed for {path.name}: {mismatches}")
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "title": title,
        "candidate": candidate,
        "sequence": recorded_sequence,
        "formal_charge": computed_charge,
        "recorded_formal_charge": recorded_charge,
        "atom_count": molecule.GetNumAtoms(),
        "heavy_atom_count": sum(atom.GetAtomicNum() > 1 for atom in molecule.GetAtoms()),
        "conformer_count": molecule.GetNumConformers(),
        "source_conformer_id": (
            int(molecule.GetProp("SourceConformerID"))
            if molecule.HasProp("SourceConformerID")
            else None
        ),
        "mmff94s_energy": (
            float(molecule.GetProp("MMFF94s_Energy"))
            if molecule.HasProp("MMFF94s_Energy")
            else None
        ),
        "coordinates_regenerated": False,
        "chemistry_checks": chemistry_checks,
    }


def prepare_backbone_rigid_setup(source_sdf: Path):
    from meeko import MoleculePreparation

    molecule = read_single_sdf(source_sdf)
    preparation = MoleculePreparation(
        rigidify_bonds_smarts=[PEPTIDE_BACKBONE_SMARTS, PEPTIDE_BACKBONE_SMARTS],
        rigidify_bonds_indices=[(0, 1), (1, 2)],
        add_index_map=True,
    )
    setups = preparation.prepare(molecule)
    if len(setups) != 1:
        raise RuntimeError(f"Expected one Meeko ligand setup for {source_sdf}; got {len(setups)}")
    return setups[0]


def prepare_backbone_rigid_ligand(source_sdf: Path, output_pdbqt: Path) -> dict[str, Any]:
    from meeko import PDBQTWriterLegacy

    if output_pdbqt.exists():
        raise FileExistsError(f"Refusing to overwrite existing PDBQT: {output_pdbqt}")
    source_hash = sha256(source_sdf)
    setup = prepare_backbone_rigid_setup(source_sdf)
    pdbqt_text, ok, error = PDBQTWriterLegacy.write_string(
        setup, add_index_map=True, remove_smiles=False
    )
    if not ok:
        raise RuntimeError(f"Meeko could not write {source_sdf}: {error}")
    torsion_lines = [
        line for line in pdbqt_text.splitlines() if line.startswith("TORSDOF ")
    ]
    if len(torsion_lines) != 1:
        raise RuntimeError(f"Expected one TORSDOF record in generated PDBQT for {source_sdf}")
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    with output_pdbqt.open("x", encoding="ascii", newline="\n") as handle:
        handle.write(pdbqt_text)
    if sha256(source_sdf) != source_hash:
        raise RuntimeError(f"Source SDF changed during preparation: {source_sdf}")
    return {
        "source_sdf_sha256": source_hash,
        "generated_pdbqt_sha256": sha256(output_pdbqt),
        "torsdof": int(torsion_lines[0].split()[1]),
        "coordinates_regenerated": False,
        "meeko_index_map_present": "REMARK INDEX MAP" in pdbqt_text,
    }


def parse_pdbqt_index_map(text: str) -> dict[int, int]:
    mapping: dict[int, int] = {}
    for line in text.splitlines():
        if not line.startswith("REMARK INDEX MAP"):
            continue
        values = [int(value) for value in line.split()[3:]]
        if len(values) % 2:
            raise ValueError("Malformed REMARK INDEX MAP record")
        for input_index, pdbqt_serial in zip(values[0::2], values[1::2]):
            if input_index in mapping:
                raise ValueError(f"Duplicate SDF atom index in PDBQT map: {input_index}")
            mapping[input_index] = pdbqt_serial
    if not mapping:
        raise ValueError("Generated PDBQT has no REMARK INDEX MAP")
    return mapping


def compare_pdbqt_coordinates(source_sdf: Path, ligand_pdbqt: Path) -> dict[str, Any]:
    molecule = read_single_sdf(source_sdf)
    conformer = molecule.GetConformer()
    text = ligand_pdbqt.read_text(encoding="ascii", errors="strict")
    index_map = parse_pdbqt_index_map(text)
    pdbqt_coordinates: dict[int, tuple[float, float, float]] = {}
    for line in text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            pdbqt_coordinates[int(line[6:11])] = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
    deltas: list[float] = []
    for input_index, pdbqt_serial in index_map.items():
        if pdbqt_serial not in pdbqt_coordinates:
            raise ValueError(f"Mapped PDBQT atom {pdbqt_serial} has no coordinate record")
        position = conformer.GetAtomPosition(input_index - 1)
        deltas.append(
            dist(
                (position.x, position.y, position.z),
                pdbqt_coordinates[pdbqt_serial],
            )
        )
    maximum = max(deltas)
    status = "MATCH_WITHIN_PDBQT_PRECISION" if maximum <= 0.001 else "MISMATCH"
    if status != "MATCH_WITHIN_PDBQT_PRECISION":
        raise ValueError(
            f"PDBQT coordinates differ from source SDF beyond 0.001 A: max {maximum:.9f}"
        )
    return {
        "mapped_atom_count": len(deltas),
        "maximum_coordinate_delta_A": maximum,
        "mean_coordinate_delta_A": sum(deltas) / len(deltas),
        "comparison_tolerance_A": 0.001,
        "status": status,
        "note": "PDBQT coordinates are serialized to 0.001 A precision",
    }


def audit_ligand_pdbqt(
    path: Path, source_sdf: Path | None = None, sequence: str | None = None
) -> dict[str, Any]:
    text = path.read_text(encoding="ascii", errors="strict")
    atom_lines = [
        line for line in text.splitlines() if line.startswith(("ATOM  ", "HETATM"))
    ]
    if not atom_lines:
        raise ValueError(f"PDBQT contains no atom records: {path}")
    atom_types = [line.split()[-1] for line in atom_lines]
    heavy_count = sum(not atom_type.upper().startswith("H") for atom_type in atom_types)
    torsion_lines = [line for line in text.splitlines() if line.startswith("TORSDOF ")]
    if len(torsion_lines) != 1:
        raise ValueError(f"PDBQT must contain exactly one TORSDOF record: {path}")
    index_map = parse_pdbqt_index_map(text)
    result: dict[str, Any] = {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "atom_record_count": len(atom_lines),
        "heavy_atom_count": heavy_count,
        "torsdof": int(torsion_lines[0].split()[1]),
        "index_map_entry_count": len(index_map),
    }
    if source_sdf is not None:
        from rdkit import Chem

        if not sequence:
            raise ValueError("sequence is required when auditing PDBQT residue mapping")
        molecule = read_single_sdf(source_sdf)
        template = Chem.MolFromSequence(sequence)
        if template is None:
            raise ValueError(f"Cannot construct peptide template for sequence {sequence!r}")
        source_heavy = [atom for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
        template_heavy = [atom for atom in template.GetAtoms() if atom.GetAtomicNum() > 1]
        if [atom.GetSymbol() for atom in source_heavy] != [
            atom.GetSymbol() for atom in template_heavy
        ]:
            raise ValueError("Standardized SDF heavy-atom order differs from sequence template")
        residue_by_input_index = {
            source_atom.GetIdx() + 1: template_atom.GetPDBResidueInfo().GetResidueNumber()
            for source_atom, template_atom in zip(source_heavy, template_heavy)
        }
        mapped_heavy_serials: dict[int, int] = {}
        for input_index, pdbqt_serial in index_map.items():
            atom = molecule.GetAtomWithIdx(input_index - 1)
            if atom.GetAtomicNum() <= 1:
                continue
            mapped_heavy_serials[pdbqt_serial] = residue_by_input_index[input_index]
        result.update(
            {
                "mapped_heavy_atom_count": len(mapped_heavy_serials),
                "mapped_residue_numbers": sorted(set(mapped_heavy_serials.values())),
            }
        )
    return result


def read_box_registry(path: Path, template_id: str) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["template"] == template_id]
    if len(rows) != 1:
        raise ValueError(f"Expected one box registry row for {template_id}; found {len(rows)}")
    row = rows[0]
    return {
        "template_id": template_id,
        "center_A": tuple(float(row[f"center_{axis}"]) for axis in "xyz"),
        "size_A": tuple(float(row[f"size_{axis}"]) for axis in "xyz"),
        "definition": row["definition"],
        "source_path": str(path.resolve()),
        "source_sha256": sha256(path),
    }


def comparison_protocol_id(protocol: dict[str, Any]) -> str:
    canonical = json.dumps(protocol, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return "v04_" + hashlib.sha256(canonical.encode("ascii")).hexdigest()[:20]


def write_vina_config(
    path: Path,
    receptor_relative: str,
    ligand_relative: str,
    box: dict[str, Any],
    settings: dict[str, int],
) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing Vina config: {path}")
    center = box["center_A"]
    size = box["size_A"]
    fields: list[tuple[str, Any]] = [
        ("receptor", receptor_relative),
        ("ligand", ligand_relative),
        ("center_x", center[0]),
        ("center_y", center[1]),
        ("center_z", center[2]),
        ("size_x", size[0]),
        ("size_y", size[1]),
        ("size_z", size[2]),
        ("exhaustiveness", settings["exhaustiveness"]),
        ("num_modes", settings["num_modes"]),
        ("energy_range", settings["energy_range"]),
        ("seed", settings["seed"]),
        ("cpu", settings["cpu"]),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write("".join(f"{key} = {value}\n" for key, value in fields))
