"""Generate explicitly provenance-labelled legacy pose metrics for standardized runs."""

from __future__ import annotations

import csv
from math import dist
from pathlib import Path
from typing import Any

from peptide_score_summary import parse_vina_scores
from pose_analysis_core import (
    CoordinateAtom,
    parse_pdbqt_atoms,
    parse_receptor_pdb_atoms,
    parse_vina_model_blocks,
)
from standardized_vina_inputs import parse_pdbqt_index_map, read_single_sdf


METRIC_COLUMNS = (
    "model",
    "score",
    "contact_fraction",
    "contacted_peptide_residues",
    "receptor_contact_residues",
    "target_contacts",
    "target_contact_count",
    "receptor_clash_pairs",
    "min_receptor_distance",
    "chain_b_clash_pairs",
    "min_chain_b_distance",
    "basic_geometry_pass",
    "contact_definition",
    "contact_fraction_provenance",
    "basic_geometry_pass_provenance",
    "receptor_clash_provenance",
    "chain_b_clash_provenance",
    "vina_score_provenance",
)


def ligand_serial_to_residue(
    source_sdf: Path, ligand_pdbqt: Path, sequence: str
) -> dict[int, int]:
    from rdkit import Chem

    source = read_single_sdf(source_sdf)
    template = Chem.MolFromSequence(sequence)
    if template is None:
        raise ValueError(f"Cannot construct sequence template: {sequence!r}")
    source_heavy = [atom for atom in source.GetAtoms() if atom.GetAtomicNum() > 1]
    template_heavy = [atom for atom in template.GetAtoms() if atom.GetAtomicNum() > 1]
    if [atom.GetSymbol() for atom in source_heavy] != [
        atom.GetSymbol() for atom in template_heavy
    ]:
        raise ValueError(
            "Standardized SDF heavy-atom order differs from sequence-derived template"
        )
    residue_by_input_index = {
        source_atom.GetIdx() + 1: template_atom.GetPDBResidueInfo().GetResidueNumber()
        for source_atom, template_atom in zip(source_heavy, template_heavy)
    }
    index_map = parse_pdbqt_index_map(
        ligand_pdbqt.read_text(encoding="ascii", errors="strict")
    )
    mapping = {
        pdbqt_serial: residue_by_input_index[input_index]
        for input_index, pdbqt_serial in index_map.items()
        if input_index in residue_by_input_index
    }
    if len(mapping) != len(source_heavy):
        raise ValueError(
            f"Mapped {len(mapping)} heavy atoms, expected {len(source_heavy)}"
        )
    return mapping


def compute_legacy_pose_metrics(
    *,
    model: int,
    score: float,
    ligand_atoms: list[CoordinateAtom],
    receptor_atoms: list[CoordinateAtom],
    chain_b_atoms: list[CoordinateAtom],
    ligand_residue_by_serial: dict[int, int],
    peptide_length: int,
    target_residues: set[int],
) -> dict[str, Any]:
    ligand_heavy = [atom for atom in ligand_atoms if atom.is_heavy]
    receptor_heavy = [atom for atom in receptor_atoms if atom.is_heavy]
    chain_b_heavy = [atom for atom in chain_b_atoms if atom.is_heavy]
    if not ligand_heavy or not receptor_heavy or not chain_b_heavy:
        raise ValueError("Ligand, receptor chain A, and receptor chain B need heavy atoms")

    contacted_peptide: set[int] = set()
    contacted_receptor: set[int] = set()
    contacted_targets: set[int] = set()
    receptor_clashes = 0
    chain_b_clashes = 0
    min_receptor = float("inf")
    min_chain_b = float("inf")
    for ligand in ligand_heavy:
        ligand_residue = ligand_residue_by_serial.get(ligand.serial)
        if ligand_residue is None:
            raise ValueError(f"No peptide residue mapping for ligand atom {ligand.serial}")
        for receptor in receptor_heavy:
            distance = dist(ligand.xyz, receptor.xyz)
            min_receptor = min(min_receptor, distance)
            receptor_clashes += int(distance < 2.0)
            if distance <= 5.0:
                contacted_peptide.add(ligand_residue)
                assert receptor.residue_number is not None
                contacted_receptor.add(receptor.residue_number)
                if receptor.residue_number in target_residues:
                    contacted_targets.add(receptor.residue_number)
        for receptor in chain_b_heavy:
            distance = dist(ligand.xyz, receptor.xyz)
            min_chain_b = min(min_chain_b, distance)
            chain_b_clashes += int(distance < 2.0)

    metrics: dict[str, Any] = {
        "model": model,
        "score": score,
        "contact_fraction": len(contacted_peptide) / peptide_length,
        "contacted_peptide_residues": ";".join(map(str, sorted(contacted_peptide))),
        "receptor_contact_residues": ";".join(map(str, sorted(contacted_receptor))),
        "target_contacts": ";".join(map(str, sorted(contacted_targets))),
        "target_contact_count": len(contacted_targets),
        "receptor_clash_pairs": receptor_clashes,
        "min_receptor_distance": min_receptor,
        "chain_b_clash_pairs": chain_b_clashes,
        "min_chain_b_distance": min_chain_b,
        "contact_definition": "legacy_heavy_atoms_distance_le_5A",
        "contact_fraction_provenance": "V0.4 recalculation of legacy 5A whole-receptor metric",
        "basic_geometry_pass_provenance": "V0.4 recalculation of legacy composite metric",
        "receptor_clash_provenance": "V0.4 heavy-atom distance_lt_2A versus receptor monomer chain A",
        "chain_b_clash_provenance": "V0.4 heavy-atom distance_lt_2A versus receptor dimer chain B",
        "vina_score_provenance": "docking_output_pdbqt_remark",
    }
    metrics["basic_geometry_pass"] = int(
        metrics["contact_fraction"] >= 0.5
        and metrics["target_contact_count"] >= 2
        and receptor_clashes == 0
        and chain_b_clashes == 0
    )
    return metrics


def generate_pose_metrics(
    *,
    docking_output: Path,
    source_sdf: Path,
    ligand_pdbqt: Path,
    receptor_monomer: Path,
    receptor_dimer: Path,
    receptor_chain: str,
    dimer_chain: str,
    sequence: str,
    target_residues: set[int],
    output_tsv: Path,
) -> list[dict[str, Any]]:
    if output_tsv.exists():
        raise FileExistsError(f"Refusing to overwrite pose metrics: {output_tsv}")
    blocks = parse_vina_model_blocks(docking_output)
    scores = parse_vina_scores(blocks)
    receptor_atoms = parse_receptor_pdb_atoms(receptor_monomer, receptor_chain)
    chain_b_atoms = parse_receptor_pdb_atoms(receptor_dimer, dimer_chain)
    mapping = ligand_serial_to_residue(source_sdf, ligand_pdbqt, sequence)
    rows = [
        compute_legacy_pose_metrics(
            model=model,
            score=scores[model],
            ligand_atoms=parse_pdbqt_atoms(block),
            receptor_atoms=receptor_atoms,
            chain_b_atoms=chain_b_atoms,
            ligand_residue_by_serial=mapping,
            peptide_length=len(sequence),
            target_residues=target_residues,
        )
        for model, block in blocks.items()
    ]
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=METRIC_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    return rows
