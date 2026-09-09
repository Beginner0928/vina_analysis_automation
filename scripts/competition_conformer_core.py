"""Deterministic Top40 competition conformer generation primitives."""

from __future__ import annotations

import hashlib
import csv
import json
import platform
import sys
from pathlib import Path
from typing import Any

from candidate_manifest_v05 import (
    load_candidate_manifest,
    validate_json_schema,
)


PROTOCOL_ID = "tfr1_top40_sdfgen_competition_v1_b3derived"
PROTOCOL_CLASSIFICATION = "COMPETITION_EXPLORATORY_SCREENING"
FROZEN_EMBEDDING = {
    "method": "ETKDGv3",
    "requested_source_conformer_count": 30,
    "use_random_coords": True,
    "prune_rms_threshold_A": 0.25,
    "enforce_chirality": True,
    "num_threads": 1,
}
FROZEN_OPTIMIZATION = {
    "force_field": "MMFF94s",
    "max_iterations": 1000,
    "num_threads": 1,
    "eligible_convergence_status": 0,
    "fallback": "none",
}
FROZEN_THRESHOLDS = [2.0, 1.5, 1.0, 0.5]
FROZEN_CANDIDATE_AUTHORITY = {
    "manifest_id": "tfr1_top40_v05",
    "path": "vina_analysis_automation/specs/top40_candidates_v05.json",
    "path_base": "data_root",
    "sha256": "297c36027fdab031c70f799ca589635bf10349202df065ec0edfb6ed804a112e",
}
FROZEN_CHEMISTRY_CONTRACT = {
    "contract_id": "v04_formal_standardized_peptide_chemistry",
    "path": "docking_test/vina_standardized_v04_20260907/chemistry_contract_v04.json",
    "path_base": "data_root",
    "sha256": "dd849afa0baa203fbbd4e2c6104b051edf7489b8da8a08fac9305d5696b04a7b",
}
FROZEN_SEED = {
    "algorithm": "sha256_protocol_colon_candidate_first8_mod2147483646_plus1",
    "material": "protocol_id + ':' + candidate_id",
    "minimum": 1,
    "maximum": 2147483646,
}
FROZEN_RMSD = {
    "atoms": "heavy_atoms_only",
    "correspondence": "fixed_identical_atom_indices",
    "alignment": "optimal_rigid_body_kabsch",
    "rdkit_api": "rdMolAlign.GetAlignmentTransform",
    "hydrogens_included": False,
}
FROZEN_SELECTION = {
    "algorithm": "energy_ordered_greedy_all_selected_pairwise_threshold",
    "energy_order": "MMFF94s_energy_ascending_then_source_conformer_id",
    "final_conformer_count": 3,
    "formal_conformer_names": ["conf01", "conf02", "conf03"],
    "rmsd_threshold_ladder_A": FROZEN_THRESHOLDS,
    "threshold_policy": "highest_full_set_success",
    "source_conformer_id_semantics": "original_rdkit_conformer_id",
}
FROZEN_FAIL_CLOSED = {
    "fewer_than_three_converged": True,
    "no_three_member_set_at_0_5_A": True,
    "chemistry_mismatch": True,
    "existing_production_protocol_directory": True,
    "force_overwrite_supported": False,
}
FROZEN_OUTPUT = {
    "root_template": "standardized_ligands/{protocol_id}",
    "candidate_directory": "{candidate_id}",
    "sdf_filename": "{candidate_id}_{formal_conformer_name}.sdf",
    "manifest_filename": "preparation_manifest.json",
    "summary_filename": "generation_summary.tsv",
    "hash_filename": "SHA256SUMS.tsv",
}
FROZEN_SDF_PROVENANCE_FIELDS = [
    "CandidateID",
    "Group",
    "Sequence",
    "PeptideLength",
    "FormalCharge",
    "AtomCount",
    "HeavyAtomCount",
    "GenerationProtocolID",
    "GenerationProtocolSHA256",
    "RDKitVersion",
    "PythonVersion",
    "Platform",
    "GenerationSeed",
    "RequestedSourceConformerCount",
    "GeneratedSourceConformerCount",
    "MMFFConvergedConformerCount",
    "MMFFVariant",
    "MMFFMaxIterations",
    "MMFFConvergenceStatus",
    "SourceConformerID",
    "FormalConformerOrdinal",
    "MMFF94sEnergy",
    "SelectionAlgorithm",
    "SelectionRMSDThresholdUsed",
    "PairwiseRMSD_conf01_conf02_A",
    "PairwiseRMSD_conf01_conf03_A",
    "PairwiseRMSD_conf02_conf03_A",
    "ChemistryContractID",
    "ChemistryContractSHA256",
    "BondGraphHash",
    "StereochemistryHash",
]
SUMMARY_FIELDS = [
    "candidate_id",
    "group",
    "sequence",
    "peptide_length",
    "formal_charge",
    "generation_seed",
    "requested_source_conformer_count",
    "generated_source_conformer_count",
    "mmff_converged_conformer_count",
    "selected_source_conformer_ids",
    "selected_mmff94s_energies",
    "selection_rmsd_threshold_used_A",
    "rmsd_conf01_conf02_A",
    "rmsd_conf01_conf03_A",
    "rmsd_conf02_conf03_A",
    "status",
]
STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def _write_json_exclusive(path: Path, value: dict[str, Any]) -> None:
    text = json.dumps(value, indent=2, ensure_ascii=False) + "\n"
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def _write_tsv_exclusive(
    path: Path, fieldnames: list[str], rows: list[dict[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def protocol_identity_sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def derive_candidate_seed(protocol_id: str, candidate_id: str) -> int:
    material = f"{protocol_id}:{candidate_id}".encode("utf-8")
    prefix = hashlib.sha256(material).hexdigest()[:8]
    return 1 + (int(prefix, 16) % 2147483646)


def expected_formal_charge(sequence: str) -> int:
    invalid = sorted(set(sequence) - STANDARD_AMINO_ACIDS)
    if not sequence or invalid:
        raise ValueError(f"Sequence must contain only canonical amino acids; invalid={invalid}")
    # The protonated N and deprotonated C termini cancel.
    return sequence.count("K") + sequence.count("R") - sequence.count("D") - sequence.count("E")


def _site_map(molecule) -> dict[tuple[int, str], Any]:
    sites: dict[tuple[int, str], Any] = {}
    for atom in molecule.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None or atom.GetAtomicNum() == 1:
            continue
        key = (info.GetResidueNumber(), info.GetName().strip())
        if key in sites:
            raise ValueError(f"Duplicate peptide atom site {key}")
        sites[key] = atom
    return sites


def _hydrogen_neighbor_count(atom) -> int:
    return sum(neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors())


def build_standardized_peptide(sequence: str):
    """Build one canonical peptide using the frozen project protonation convention."""
    from rdkit import Chem

    expected_formal_charge(sequence)
    molecule = Chem.MolFromSequence(sequence)
    if molecule is None:
        raise ValueError(f"RDKit could not construct peptide sequence {sequence!r}")
    editable = Chem.RWMol(molecule)
    residue_count = len(sequence)
    for atom in editable.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            continue
        residue_number = info.GetResidueNumber()
        residue_name = info.GetResidueName().strip()
        atom_name = info.GetName().strip()
        if residue_number == 1 and atom_name == "N":
            heavy_bond_order = int(
                sum(
                    bond.GetBondTypeAsDouble()
                    for bond in atom.GetBonds()
                    if bond.GetOtherAtom(atom).GetAtomicNum() > 1
                )
            )
            terminal_hydrogens = 4 - heavy_bond_order
            if terminal_hydrogens not in (2, 3):
                raise ValueError(
                    f"Unsupported N-terminal bonding for sequence {sequence!r}: "
                    f"heavy bond order {heavy_bond_order}"
                )
            atom.SetFormalCharge(1)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(terminal_hydrogens)
        elif residue_number == residue_count and atom_name == "OXT":
            atom.SetFormalCharge(-1)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)
        elif (residue_name, atom_name) in {("ASP", "OD2"), ("GLU", "OE2")}:
            atom.SetFormalCharge(-1)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)
        elif (residue_name, atom_name) == ("LYS", "NZ"):
            atom.SetFormalCharge(1)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(3)
        elif (residue_name, atom_name) == ("ARG", "NH1"):
            atom.SetFormalCharge(1)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(2)
    standardized = editable.GetMol()
    Chem.SanitizeMol(standardized)
    standardized = Chem.AddHs(standardized)
    Chem.AssignStereochemistry(standardized, cleanIt=True, force=True)
    return standardized


def bond_graph_hash(molecule) -> str:
    heavy_indices = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
    normalized_index = {atom_index: ordinal for ordinal, atom_index in enumerate(heavy_indices)}
    atoms = []
    for atom_index in heavy_indices:
        atom = molecule.GetAtomWithIdx(atom_index)
        atoms.append(
            [
                atom.GetAtomicNum(),
                atom.GetFormalCharge(),
            ]
        )
    bonds = []
    for bond in molecule.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin in normalized_index and end in normalized_index:
            left, right = sorted((normalized_index[begin], normalized_index[end]))
            bonds.append([left, right, str(bond.GetBondType())])
    value = {"atoms": atoms, "bonds": sorted(bonds)}
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def stereochemistry_hash(molecule) -> str:
    from rdkit import Chem

    copy = Chem.Mol(molecule)
    Chem.AssignStereochemistry(copy, cleanIt=True, force=True)
    values = []
    for atom in copy.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        values.append(
            [
                atom.GetIdx(),
                str(atom.GetChiralTag()),
                atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None,
            ]
        )
    encoded = json.dumps(values, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def audit_standardized_peptide(molecule, candidate: dict[str, Any]) -> dict[str, Any]:
    from rdkit import Chem

    candidate_id = str(candidate["candidate_id"])
    sequence = str(candidate["sequence"])
    peptide_length = candidate["peptide_length"]
    if peptide_length != len(sequence):
        raise ValueError(f"Peptide length mismatch for {candidate_id}")
    sulfur_sulfur_bonds = sum(
        bond.GetBeginAtom().GetAtomicNum() == 16 and bond.GetEndAtom().GetAtomicNum() == 16
        for bond in molecule.GetBonds()
    )
    if sulfur_sulfur_bonds:
        raise ValueError(f"Unintended S-S bond detected for {candidate_id}")

    expected = build_standardized_peptide(sequence)
    actual_heavy = [atom for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
    expected_heavy = [atom for atom in expected.GetAtoms() if atom.GetAtomicNum() > 1]
    actual_identity = [atom.GetSymbol() for atom in actual_heavy]
    expected_identity = [atom.GetSymbol() for atom in expected_heavy]
    if actual_identity != expected_identity:
        raise ValueError(f"Peptide heavy-atom identity/order mismatch for {candidate_id}")
    if molecule.GetNumAtoms() != expected.GetNumAtoms():
        raise ValueError(f"Peptide atom count mismatch for {candidate_id}")
    if bond_graph_hash(molecule) != bond_graph_hash(expected):
        raise ValueError(f"Peptide bond graph mismatch for {candidate_id}")
    formal_charge = Chem.GetFormalCharge(molecule)
    expected_charge = expected_formal_charge(sequence)
    if formal_charge != expected_charge:
        raise ValueError(
            f"Peptide formal charge mismatch for {candidate_id}: {formal_charge} != {expected_charge}"
        )

    expected_sites = _site_map(expected)
    actual_sites = {}
    for actual_atom, expected_atom in zip(actual_heavy, expected_heavy):
        info = expected_atom.GetPDBResidueInfo()
        if info is None:
            raise ValueError(f"Expected peptide template lacks residue metadata for {candidate_id}")
        actual_sites[(info.GetResidueNumber(), info.GetName().strip())] = actual_atom
    charge_and_hydrogen_mismatches = {}
    for key, expected_atom in expected_sites.items():
        actual_atom = actual_sites[key]
        actual_state = (actual_atom.GetFormalCharge(), _hydrogen_neighbor_count(actual_atom))
        expected_state = (expected_atom.GetFormalCharge(), _hydrogen_neighbor_count(expected_atom))
        if actual_state != expected_state:
            charge_and_hydrogen_mismatches[str(key)] = {
                "expected": expected_state,
                "actual": actual_state,
            }
    if charge_and_hydrogen_mismatches:
        raise ValueError(
            f"Peptide protonation mismatch for {candidate_id}: {charge_and_hydrogen_mismatches}"
        )

    actual_chiral = [str(atom.GetChiralTag()) for atom in actual_heavy]
    expected_chiral = [str(atom.GetChiralTag()) for atom in expected_heavy]
    standard_l = actual_chiral == expected_chiral
    if not standard_l:
        raise ValueError(f"Peptide standard-L stereochemistry mismatch for {candidate_id}")

    return {
        "status": "PASS",
        "candidate_id": candidate_id,
        "group": candidate["group"],
        "sequence": sequence,
        "peptide_length": peptide_length,
        "formal_charge": formal_charge,
        "atom_count": molecule.GetNumAtoms(),
        "heavy_atom_count": len(actual_heavy),
        "bond_graph_hash": bond_graph_hash(molecule),
        "stereochemistry_hash": stereochemistry_hash(molecule),
        "chemistry_checks": {
            "n_terminus_protonated": (
                actual_sites[(1, "N")].GetFormalCharge() == 1
                and _hydrogen_neighbor_count(actual_sites[(1, "N")])
                == _hydrogen_neighbor_count(expected_sites[(1, "N")])
            ),
            "c_terminus_deprotonated": (
                actual_sites[(len(sequence), "OXT")].GetFormalCharge() == -1
                and _hydrogen_neighbor_count(actual_sites[(len(sequence), "OXT")]) == 0
            ),
            "histidine_HIE_like": all(
                _hydrogen_neighbor_count(actual_sites[(index, "ND1")]) == 0
                and _hydrogen_neighbor_count(actual_sites[(index, "NE2")]) == 1
                for index, code in enumerate(sequence, start=1)
                if code == "H"
            ),
            "cysteine_free_thiol": all(
                actual_sites[(index, "SG")].GetFormalCharge() == 0
                and _hydrogen_neighbor_count(actual_sites[(index, "SG")]) == 1
                for index, code in enumerate(sequence, start=1)
                if code == "C"
            ),
            "tyrosine_neutral": all(
                actual_sites[(index, "OH")].GetFormalCharge() == 0
                and _hydrogen_neighbor_count(actual_sites[(index, "OH")]) == 1
                for index, code in enumerate(sequence, start=1)
                if code == "Y"
            ),
            "standard_L_stereochemistry": standard_l,
            "sulfur_sulfur_bond_count": sulfur_sulfur_bonds,
        },
    }


def build_etkdg_parameters(protocol: dict[str, Any], seed: int):
    from rdkit.Chem import AllChem

    validate_generation_protocol(protocol)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = seed
    parameters.useRandomCoords = protocol["embedding"]["use_random_coords"]
    parameters.pruneRmsThresh = protocol["embedding"]["prune_rms_threshold_A"]
    parameters.enforceChirality = protocol["embedding"]["enforce_chirality"]
    parameters.numThreads = protocol["embedding"]["num_threads"]
    return parameters


def heavy_atom_aligned_rmsd(molecule, left_conformer_id: int, right_conformer_id: int) -> float:
    """Return Kabsch-aligned RMSD using fixed corresponding heavy-atom indices."""
    from rdkit.Chem import rdMolAlign

    atom_map = [
        (atom.GetIdx(), atom.GetIdx())
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() > 1
    ]
    rmsd, _transform = rdMolAlign.GetAlignmentTransform(
        molecule,
        molecule,
        prbCid=left_conformer_id,
        refCid=right_conformer_id,
        atomMap=atom_map,
        reflect=False,
        maxIters=50,
    )
    return float(rmsd)


def select_formal_conformers(
    optimization_rows: list[dict[str, Any]],
    rmsd_function,
    protocol: dict[str, Any],
) -> dict[str, Any]:
    validate_generation_protocol(protocol)
    eligible_status = protocol["optimization"]["eligible_convergence_status"]
    eligible = sorted(
        (
            row
            for row in optimization_rows
            if row["convergence_status"] == eligible_status
        ),
        key=lambda row: (float(row["energy"]), int(row["source_conformer_id"])),
    )
    required_count = protocol["selection"]["final_conformer_count"]
    if len(eligible) < required_count:
        raise RuntimeError(
            f"Candidate has fewer than 3 converged MMFF94s conformers: {len(eligible)}"
        )

    rmsd_cache: dict[tuple[int, int], float] = {}

    def distance(left: int, right: int) -> float:
        key = tuple(sorted((left, right)))
        if key not in rmsd_cache:
            rmsd_cache[key] = float(rmsd_function(left, right))
        return rmsd_cache[key]

    selected_rows: list[dict[str, Any]] | None = None
    threshold_used: float | None = None
    for threshold in protocol["selection"]["rmsd_threshold_ladder_A"]:
        trial = [eligible[0]]
        for row in eligible[1:]:
            candidate_id = int(row["source_conformer_id"])
            if all(
                distance(candidate_id, int(selected["source_conformer_id"])) >= threshold
                for selected in trial
            ):
                trial.append(row)
            if len(trial) == required_count:
                selected_rows = trial
                threshold_used = float(threshold)
                break
        if selected_rows is not None:
            break
    if selected_rows is None or threshold_used is None:
        raise RuntimeError("Could not select three conformers at the minimum 0.5 A threshold")

    selected_ids = [int(row["source_conformer_id"]) for row in selected_rows]
    formal_names = protocol["selection"]["formal_conformer_names"]
    pairwise = {}
    for left_index in range(required_count):
        for right_index in range(left_index + 1, required_count):
            pairwise[f"{formal_names[left_index]}-{formal_names[right_index]}"] = distance(
                selected_ids[left_index], selected_ids[right_index]
            )
    return {
        "selected_source_conformer_ids": selected_ids,
        "selected_energies": [float(row["energy"]) for row in selected_rows],
        "rmsd_threshold_used_A": threshold_used,
        "pairwise_rmsd_A": pairwise,
        "eligible_conformer_count": len(eligible),
    }


def generate_candidate_conformers(
    candidate: dict[str, Any], protocol: dict[str, Any]
) -> dict[str, Any]:
    import rdkit
    from rdkit.Chem import AllChem

    validate_generation_protocol(protocol)
    expected_version = protocol["software"]["rdkit_version"]
    if rdkit.__version__ != expected_version:
        raise RuntimeError(
            f"RDKit version mismatch: expected {expected_version}, actual {rdkit.__version__}"
        )
    molecule = build_standardized_peptide(candidate["sequence"])
    chemistry_audit = audit_standardized_peptide(molecule, candidate)
    seed = derive_candidate_seed(protocol["protocol_id"], candidate["candidate_id"])
    parameters = build_etkdg_parameters(protocol, seed)
    requested = protocol["embedding"]["requested_source_conformer_count"]
    conformer_ids = list(
        AllChem.EmbedMultipleConfs(molecule, numConfs=requested, params=parameters)
    )
    if len(conformer_ids) < 3:
        raise RuntimeError(
            f"Candidate generated fewer than 3 source conformers: {len(conformer_ids)}"
        )
    optimization = protocol["optimization"]
    results = list(
        AllChem.MMFFOptimizeMoleculeConfs(
            molecule,
            numThreads=optimization["num_threads"],
            mmffVariant=optimization["force_field"],
            maxIters=optimization["max_iterations"],
        )
    )
    if len(results) != len(conformer_ids):
        raise RuntimeError(
            "RDKit MMFF result count differs from the generated conformer inventory"
        )
    rows = [
        {
            "source_conformer_id": int(conformer_id),
            "convergence_status": int(status),
            "energy": float(energy),
        }
        for conformer_id, (status, energy) in zip(conformer_ids, results)
    ]
    selection = select_formal_conformers(
        rows,
        lambda left, right: heavy_atom_aligned_rmsd(molecule, left, right),
        protocol,
    )
    return {
        "molecule": molecule,
        "candidate": candidate,
        "chemistry_audit": chemistry_audit,
        "generation_seed": seed,
        "requested_source_conformer_count": requested,
        "generated_source_conformer_count": len(conformer_ids),
        "mmff_converged_conformer_count": sum(
            row["convergence_status"] == optimization["eligible_convergence_status"]
            for row in rows
        ),
        "optimization_rows": rows,
        "selection": selection,
    }


def _format_float(value: float) -> str:
    return f"{value:.12f}"


def _sdf_properties(
    result: dict[str, Any],
    protocol: dict[str, Any],
    protocol_sha256: str,
    chemistry_contract: dict[str, Any],
    chemistry_contract_sha256: str,
    formal_name: str,
    source_conformer_id: int,
    energy: float,
) -> list[tuple[str, str]]:
    audit = result["chemistry_audit"]
    pairwise = result["selection"]["pairwise_rmsd_A"]
    return [
        ("CandidateID", result["candidate"]["candidate_id"]),
        ("Group", result["candidate"]["group"]),
        ("Sequence", result["candidate"]["sequence"]),
        ("PeptideLength", str(result["candidate"]["peptide_length"])),
        ("FormalCharge", str(audit["formal_charge"])),
        ("AtomCount", str(audit["atom_count"])),
        ("HeavyAtomCount", str(audit["heavy_atom_count"])),
        ("GenerationProtocolID", protocol["protocol_id"]),
        ("GenerationProtocolSHA256", protocol_sha256),
        ("RDKitVersion", protocol["software"]["rdkit_version"]),
        ("PythonVersion", platform.python_version()),
        ("Platform", platform.platform()),
        ("GenerationSeed", str(result["generation_seed"])),
        (
            "RequestedSourceConformerCount",
            str(result["requested_source_conformer_count"]),
        ),
        ("GeneratedSourceConformerCount", str(result["generated_source_conformer_count"])),
        ("MMFFConvergedConformerCount", str(result["mmff_converged_conformer_count"])),
        ("MMFFVariant", protocol["optimization"]["force_field"]),
        ("MMFFMaxIterations", str(protocol["optimization"]["max_iterations"])),
        ("MMFFConvergenceStatus", "0"),
        ("SourceConformerID", str(source_conformer_id)),
        ("FormalConformerOrdinal", formal_name),
        ("MMFF94sEnergy", _format_float(energy)),
        ("SelectionAlgorithm", protocol["selection"]["algorithm"]),
        (
            "SelectionRMSDThresholdUsed",
            _format_float(result["selection"]["rmsd_threshold_used_A"]),
        ),
        ("PairwiseRMSD_conf01_conf02_A", _format_float(pairwise["conf01-conf02"])),
        ("PairwiseRMSD_conf01_conf03_A", _format_float(pairwise["conf01-conf03"])),
        ("PairwiseRMSD_conf02_conf03_A", _format_float(pairwise["conf02-conf03"])),
        ("ChemistryContractID", chemistry_contract["contract_id"]),
        ("ChemistryContractSHA256", chemistry_contract_sha256),
        ("BondGraphHash", audit["bond_graph_hash"]),
        ("StereochemistryHash", audit["stereochemistry_hash"]),
    ]


def _write_one_sdf(
    path: Path,
    result: dict[str, Any],
    protocol: dict[str, Any],
    protocol_sha256: str,
    chemistry_contract: dict[str, Any],
    chemistry_contract_sha256: str,
    formal_name: str,
    source_conformer_id: int,
    energy: float,
) -> None:
    from rdkit import Chem

    molecule = Chem.Mol(result["molecule"])
    for name in list(molecule.GetPropNames(includePrivate=False, includeComputed=False)):
        molecule.ClearProp(name)
    molecule.SetProp("_Name", f"{result['candidate']['candidate_id']}_{formal_name}")
    for key, value in _sdf_properties(
        result,
        protocol,
        protocol_sha256,
        chemistry_contract,
        chemistry_contract_sha256,
        formal_name,
        source_conformer_id,
        energy,
    ):
        molecule.SetProp(key, value)
    writer = Chem.SDWriter(str(path))
    if writer is None:
        raise OSError(f"Could not create SDF writer for {path}")
    try:
        writer.write(molecule, confId=source_conformer_id)
    finally:
        writer.close()
    if not path.is_file() or path.stat().st_size == 0:
        raise OSError(f"RDKit did not publish SDF {path}")


def _read_single_sdf(path: Path):
    from rdkit import Chem

    molecules = [molecule for molecule in Chem.SDMolSupplier(str(path), removeHs=False) if molecule]
    if len(molecules) != 1:
        raise ValueError(f"Expected one readable molecule in {path}; found {len(molecules)}")
    return molecules[0]


def _audit_sdf_properties(
    molecule,
    candidate: dict[str, Any],
    formal_name: str,
    protocol: dict[str, Any],
    protocol_sha256: str,
    chemistry_contract_sha256: str,
) -> dict[str, Any]:
    required = set(protocol["sdf_provenance_fields"])
    actual = set(molecule.GetPropNames(includePrivate=False, includeComputed=False))
    if actual != required:
        raise ValueError(
            f"SDF provenance fields differ for {candidate['candidate_id']}_{formal_name}: "
            f"missing={sorted(required - actual)}, extra={sorted(actual - required)}"
        )
    expected_properties = {
        "CandidateID": candidate["candidate_id"],
        "Group": candidate["group"],
        "Sequence": candidate["sequence"],
        "PeptideLength": str(candidate["peptide_length"]),
        "GenerationProtocolID": protocol["protocol_id"],
        "GenerationProtocolSHA256": protocol_sha256,
        "FormalConformerOrdinal": formal_name,
        "ChemistryContractID": protocol["chemistry_contract"]["contract_id"],
        "ChemistryContractSHA256": chemistry_contract_sha256,
    }
    mismatches = {
        key: {"expected": value, "actual": molecule.GetProp(key)}
        for key, value in expected_properties.items()
        if molecule.GetProp(key) != value
    }
    if mismatches:
        raise ValueError(f"SDF provenance mismatch: {mismatches}")
    if molecule.HasProp("SDFSHA256"):
        raise ValueError("SDF raw self-hash must not be embedded in the SDF")
    chemistry_audit = audit_standardized_peptide(molecule, candidate)
    if molecule.GetProp("BondGraphHash") != chemistry_audit["bond_graph_hash"]:
        raise ValueError("SDF BondGraphHash differs from read-back chemistry audit")
    if molecule.GetProp("StereochemistryHash") != chemistry_audit["stereochemistry_hash"]:
        raise ValueError("SDF StereochemistryHash differs from read-back chemistry audit")
    return chemistry_audit


def write_generation_package(
    results: list[dict[str, Any]],
    protocol: dict[str, Any],
    protocol_sha256: str,
    chemistry_contract: dict[str, Any],
    protocol_directory: Path,
    *,
    candidate_manifest_sha256: str,
    chemistry_contract_sha256: str,
) -> Path:
    """Publish a complete package into a new directory; existing targets are immutable."""
    validate_generation_protocol(protocol)
    if protocol_directory.exists():
        raise FileExistsError(f"Refusing to overwrite existing protocol directory: {protocol_directory}")
    protocol_directory.parent.mkdir(parents=True, exist_ok=True)
    protocol_directory.mkdir()
    manifest_candidates = []
    summary_rows = []
    hash_rows = []
    formal_names = protocol["selection"]["formal_conformer_names"]
    for result in results:
        candidate = result["candidate"]
        candidate_directory = protocol_directory / candidate["candidate_id"]
        candidate_directory.mkdir()
        selected_ids = result["selection"]["selected_source_conformer_ids"]
        selected_energies = result["selection"]["selected_energies"]
        if len(selected_ids) != len(formal_names) or len(selected_energies) != len(formal_names):
            raise ValueError("Selection result does not contain exactly three formal conformers")
        formal_rows = []
        for formal_name, source_id, energy in zip(
            formal_names, selected_ids, selected_energies
        ):
            relative_path = Path(candidate["candidate_id"]) / f"{candidate['candidate_id']}_{formal_name}.sdf"
            sdf_path = protocol_directory / relative_path
            _write_one_sdf(
                sdf_path,
                result,
                protocol,
                protocol_sha256,
                chemistry_contract,
                chemistry_contract_sha256,
                formal_name,
                int(source_id),
                float(energy),
            )
            read_back = _read_single_sdf(sdf_path)
            read_back_audit = _audit_sdf_properties(
                read_back,
                candidate,
                formal_name,
                protocol,
                protocol_sha256,
                chemistry_contract_sha256,
            )
            sdf_sha256 = sha256_file(sdf_path)
            formal_rows.append(
                {
                    "formal_conformer_ordinal": formal_name,
                    "source_conformer_id": int(source_id),
                    "mmff94s_energy": float(energy),
                    "mmff_convergence_status": 0,
                    "path": relative_path.as_posix(),
                    "sha256": sdf_sha256,
                    "read_back_audit_status": read_back_audit["status"],
                }
            )
            hash_rows.append({"path": relative_path.as_posix(), "sha256": sdf_sha256})
        audit = result["chemistry_audit"]
        selection = result["selection"]
        manifest_candidates.append(
            {
                "candidate_id": candidate["candidate_id"],
                "group": candidate["group"],
                "sequence": candidate["sequence"],
                "peptide_length": candidate["peptide_length"],
                "formal_charge": audit["formal_charge"],
                "atom_count": audit["atom_count"],
                "heavy_atom_count": audit["heavy_atom_count"],
                "bond_graph_hash": audit["bond_graph_hash"],
                "stereochemistry_hash": audit["stereochemistry_hash"],
                "generation_seed": result["generation_seed"],
                "requested_source_conformer_count": result[
                    "requested_source_conformer_count"
                ],
                "generated_source_conformer_count": result[
                    "generated_source_conformer_count"
                ],
                "mmff_converged_conformer_count": result[
                    "mmff_converged_conformer_count"
                ],
                "selection_rmsd_threshold_used_A": selection["rmsd_threshold_used_A"],
                "pairwise_rmsd_A": selection["pairwise_rmsd_A"],
                "chemistry_audit": audit,
                "formal_conformers": formal_rows,
            }
        )
        summary_rows.append(
            {
                "candidate_id": candidate["candidate_id"],
                "group": candidate["group"],
                "sequence": candidate["sequence"],
                "peptide_length": candidate["peptide_length"],
                "formal_charge": audit["formal_charge"],
                "generation_seed": result["generation_seed"],
                "requested_source_conformer_count": result["requested_source_conformer_count"],
                "generated_source_conformer_count": result["generated_source_conformer_count"],
                "mmff_converged_conformer_count": result["mmff_converged_conformer_count"],
                "selected_source_conformer_ids": ";".join(map(str, selected_ids)),
                "selected_mmff94s_energies": ";".join(_format_float(value) for value in selected_energies),
                "selection_rmsd_threshold_used_A": _format_float(selection["rmsd_threshold_used_A"]),
                "rmsd_conf01_conf02_A": _format_float(selection["pairwise_rmsd_A"]["conf01-conf02"]),
                "rmsd_conf01_conf03_A": _format_float(selection["pairwise_rmsd_A"]["conf01-conf03"]),
                "rmsd_conf02_conf03_A": _format_float(selection["pairwise_rmsd_A"]["conf02-conf03"]),
                "status": "PASS",
            }
        )

    manifest = {
        "schema_version": "1.0",
        "protocol_id": protocol["protocol_id"],
        "classification": protocol["classification"],
        "generation_protocol_sha256": protocol_sha256,
        "candidate_manifest_sha256": candidate_manifest_sha256,
        "chemistry_contract_id": chemistry_contract["contract_id"],
        "chemistry_contract_sha256": chemistry_contract_sha256,
        "software": {
            "python_version": platform.python_version(),
            "python_executable": str(Path(sys.executable).resolve()),
            "rdkit_version": protocol["software"]["rdkit_version"],
            "platform": platform.platform(),
        },
        "candidate_count": len(manifest_candidates),
        "formal_sdf_count": sum(len(row["formal_conformers"]) for row in manifest_candidates),
        "sdf_provenance_fields": protocol["sdf_provenance_fields"],
        "candidates": manifest_candidates,
        "status": "PASS",
    }
    _write_tsv_exclusive(protocol_directory / "generation_summary.tsv", SUMMARY_FIELDS, summary_rows)
    _write_tsv_exclusive(protocol_directory / "SHA256SUMS.tsv", ["path", "sha256"], hash_rows)
    manifest_path = protocol_directory / "preparation_manifest.json"
    _write_json_exclusive(manifest_path, manifest)
    audit_generation_package(protocol_directory)
    return manifest_path


def audit_generation_package(protocol_directory: Path) -> dict[str, Any]:
    manifest_path = protocol_directory / "preparation_manifest.json"
    hashes_path = protocol_directory / "SHA256SUMS.tsv"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot load preparation manifest: {exc}") from exc
    if manifest.get("protocol_id") != PROTOCOL_ID or manifest.get("status") != "PASS":
        raise ValueError("Preparation manifest protocol/status mismatch")
    with hashes_path.open("r", encoding="utf-8", newline="") as handle:
        hash_rows = list(csv.DictReader(handle, delimiter="\t"))
    expected_hashes = {row["path"]: row["sha256"] for row in hash_rows}
    if len(expected_hashes) != len(hash_rows):
        raise ValueError("SHA256SUMS.tsv contains duplicate paths")
    audited = 0
    protocol_stub = {
        "protocol_id": manifest["protocol_id"],
        "chemistry_contract": {"contract_id": manifest["chemistry_contract_id"]},
        "sdf_provenance_fields": manifest["sdf_provenance_fields"],
    }
    expected_summary_rows = []
    for candidate_row in manifest.get("candidates", []):
        candidate = {
            field: candidate_row[field]
            for field in ("candidate_id", "group", "sequence", "peptide_length")
        }
        for formal in candidate_row["formal_conformers"]:
            relative = formal["path"]
            if relative not in expected_hashes:
                raise ValueError(f"Missing SHA256SUMS entry for {relative}")
            path = (protocol_directory / relative).resolve()
            try:
                path.relative_to(protocol_directory.resolve())
            except ValueError as exc:
                raise ValueError(f"SDF path escapes protocol directory: {relative}") from exc
            actual_hash = sha256_file(path)
            if actual_hash != expected_hashes[relative] or actual_hash != formal["sha256"]:
                raise ValueError(f"SDF SHA-256 mismatch for {relative}")
            molecule = _read_single_sdf(path)
            chemistry_audit = _audit_sdf_properties(
                molecule,
                candidate,
                formal["formal_conformer_ordinal"],
                protocol_stub,
                manifest["generation_protocol_sha256"],
                manifest["chemistry_contract_sha256"],
            )
            expected_sdf_values = {
                "GenerationSeed": str(candidate_row["generation_seed"]),
                "RequestedSourceConformerCount": str(
                    candidate_row["requested_source_conformer_count"]
                ),
                "GeneratedSourceConformerCount": str(
                    candidate_row["generated_source_conformer_count"]
                ),
                "MMFFConvergedConformerCount": str(
                    candidate_row["mmff_converged_conformer_count"]
                ),
                "MMFFConvergenceStatus": str(formal["mmff_convergence_status"]),
                "SourceConformerID": str(formal["source_conformer_id"]),
                "MMFF94sEnergy": _format_float(formal["mmff94s_energy"]),
                "SelectionRMSDThresholdUsed": _format_float(
                    candidate_row["selection_rmsd_threshold_used_A"]
                ),
                "PairwiseRMSD_conf01_conf02_A": _format_float(
                    candidate_row["pairwise_rmsd_A"]["conf01-conf02"]
                ),
                "PairwiseRMSD_conf01_conf03_A": _format_float(
                    candidate_row["pairwise_rmsd_A"]["conf01-conf03"]
                ),
                "PairwiseRMSD_conf02_conf03_A": _format_float(
                    candidate_row["pairwise_rmsd_A"]["conf02-conf03"]
                ),
                "FormalCharge": str(candidate_row["formal_charge"]),
                "AtomCount": str(candidate_row["atom_count"]),
                "HeavyAtomCount": str(candidate_row["heavy_atom_count"]),
                "BondGraphHash": candidate_row["bond_graph_hash"],
                "StereochemistryHash": candidate_row["stereochemistry_hash"],
            }
            sdf_mismatches = {
                field: {"expected": expected, "actual": molecule.GetProp(field)}
                for field, expected in expected_sdf_values.items()
                if molecule.GetProp(field) != expected
            }
            if chemistry_audit != candidate_row["chemistry_audit"]:
                sdf_mismatches["chemistry_audit"] = "manifest differs from SDF read-back"
            if sdf_mismatches:
                raise ValueError(
                    f"Package manifest/SDF provenance mismatch for {relative}: {sdf_mismatches}"
                )
            audited += 1
        formal_rows = candidate_row["formal_conformers"]
        expected_summary_rows.append(
            {
                "candidate_id": candidate_row["candidate_id"],
                "group": candidate_row["group"],
                "sequence": candidate_row["sequence"],
                "peptide_length": str(candidate_row["peptide_length"]),
                "formal_charge": str(candidate_row["formal_charge"]),
                "generation_seed": str(candidate_row["generation_seed"]),
                "requested_source_conformer_count": str(
                    candidate_row["requested_source_conformer_count"]
                ),
                "generated_source_conformer_count": str(
                    candidate_row["generated_source_conformer_count"]
                ),
                "mmff_converged_conformer_count": str(
                    candidate_row["mmff_converged_conformer_count"]
                ),
                "selected_source_conformer_ids": ";".join(
                    str(row["source_conformer_id"]) for row in formal_rows
                ),
                "selected_mmff94s_energies": ";".join(
                    _format_float(row["mmff94s_energy"]) for row in formal_rows
                ),
                "selection_rmsd_threshold_used_A": _format_float(
                    candidate_row["selection_rmsd_threshold_used_A"]
                ),
                "rmsd_conf01_conf02_A": _format_float(
                    candidate_row["pairwise_rmsd_A"]["conf01-conf02"]
                ),
                "rmsd_conf01_conf03_A": _format_float(
                    candidate_row["pairwise_rmsd_A"]["conf01-conf03"]
                ),
                "rmsd_conf02_conf03_A": _format_float(
                    candidate_row["pairwise_rmsd_A"]["conf02-conf03"]
                ),
                "status": "PASS",
            }
        )
    if audited != manifest.get("formal_sdf_count") or audited != len(expected_hashes):
        raise ValueError("Manifest/SHA256SUMS/formal SDF counts differ")
    summary_path = protocol_directory / "generation_summary.tsv"
    with summary_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != SUMMARY_FIELDS:
            raise ValueError("generation_summary.tsv columns differ from the frozen contract")
        actual_summary_rows = list(reader)
    if actual_summary_rows != expected_summary_rows:
        raise ValueError("generation_summary.tsv differs from preparation manifest provenance")
    return {"status": "PASS", "candidate_count": len(manifest["candidates"]), "sdf_count": audited}


def audit_generation_package_against_contract(
    protocol_directory: Path, protocol_path: Path, data_root: Path
) -> dict[str, Any]:
    (
        protocol,
        protocol_sha256,
        canonical_candidates,
        manifest_sha256,
        _chemistry_contract,
        chemistry_sha256,
    ) = load_authoritative_generation_inputs(protocol_path, data_root)
    internal = audit_generation_package(protocol_directory)
    manifest = json.loads(
        (protocol_directory / protocol["output"]["manifest_filename"]).read_text(
            encoding="utf-8"
        )
    )
    expected_provenance = {
        "generation_protocol_sha256": protocol_sha256,
        "candidate_manifest_sha256": manifest_sha256,
        "chemistry_contract_id": protocol["chemistry_contract"]["contract_id"],
        "chemistry_contract_sha256": chemistry_sha256,
    }
    mismatches = {
        field: {"expected": expected, "actual": manifest.get(field)}
        for field, expected in expected_provenance.items()
        if manifest.get(field) != expected
    }
    if mismatches:
        raise ValueError(f"Package provenance differs from frozen contracts: {mismatches}")
    canonical_by_id = {row["candidate_id"]: row for row in canonical_candidates}
    seen = set()
    for row in manifest["candidates"]:
        candidate_id = row["candidate_id"]
        if candidate_id in seen:
            raise ValueError(f"Duplicate candidate in package manifest: {candidate_id}")
        seen.add(candidate_id)
        if candidate_id not in canonical_by_id:
            raise ValueError(f"Package contains non-canonical candidate: {candidate_id}")
        canonical = canonical_by_id[candidate_id]
        actual_identity = {
            field: row[field]
            for field in ("candidate_id", "group", "sequence", "peptide_length")
        }
        if actual_identity != canonical:
            raise ValueError(
                f"Package candidate identity differs from canonical manifest: {candidate_id}"
            )
        if len(row["formal_conformers"]) != protocol["selection"]["final_conformer_count"]:
            raise ValueError(f"Package candidate does not have exactly three conformers: {candidate_id}")
    return {**internal, "contract_provenance": "MATCH"}


def generate_protocol_package(
    candidates: list[dict[str, Any]],
    protocol: dict[str, Any],
    protocol_sha256: str,
    candidate_manifest_sha256: str,
    chemistry_contract: dict[str, Any],
    chemistry_contract_sha256: str,
    protocol_directory: Path,
) -> Path:
    if protocol_directory.exists():
        raise FileExistsError(f"Refusing to overwrite existing protocol directory: {protocol_directory}")
    results = [generate_candidate_conformers(candidate, protocol) for candidate in candidates]
    return write_generation_package(
        results,
        protocol,
        protocol_sha256,
        chemistry_contract,
        protocol_directory,
        candidate_manifest_sha256=candidate_manifest_sha256,
        chemistry_contract_sha256=chemistry_contract_sha256,
    )


def validate_generation_protocol(protocol: dict[str, Any]) -> dict[str, Any]:
    if protocol.get("schema_version") != "1.0":
        raise ValueError("schema_version differs from the frozen competition_v1 contract")
    if protocol.get("protocol_id") != PROTOCOL_ID:
        raise ValueError("protocol_id differs from the frozen competition_v1 identity")
    if protocol.get("classification") != PROTOCOL_CLASSIFICATION:
        raise ValueError("classification differs from COMPETITION_EXPLORATORY_SCREENING")
    if protocol.get("software") != {
        "python_version": "3.12.x",
        "rdkit_version": "2026.03.6",
    }:
        raise ValueError("software versions differ from the frozen protocol")
    if protocol.get("embedding") != FROZEN_EMBEDDING:
        raise ValueError("embedding parameters differ from the frozen protocol")
    if protocol.get("optimization") != FROZEN_OPTIMIZATION:
        raise ValueError("optimization parameters differ from the frozen protocol")
    if protocol.get("candidate_authority") != FROZEN_CANDIDATE_AUTHORITY:
        raise ValueError("candidate authority differs from the frozen canonical Top40 contract")
    if protocol.get("chemistry_contract") != FROZEN_CHEMISTRY_CONTRACT:
        raise ValueError("chemistry contract differs from the frozen project contract")
    if protocol.get("seed") != FROZEN_SEED:
        raise ValueError("seed derivation differs from the frozen protocol")
    if protocol.get("rmsd") != FROZEN_RMSD:
        raise ValueError("rmsd definition differs from the frozen protocol")
    if protocol.get("selection") != FROZEN_SELECTION:
        raise ValueError("selection definition differs from the frozen protocol")
    if protocol.get("fail_closed") != FROZEN_FAIL_CLOSED:
        raise ValueError("fail_closed definition differs from the frozen protocol")
    if protocol.get("output") != FROZEN_OUTPUT:
        raise ValueError("output definition differs from the frozen protocol")
    if protocol.get("sdf_provenance_fields") != FROZEN_SDF_PROVENANCE_FIELDS:
        raise ValueError("SDF provenance fields differ from the frozen protocol")
    return protocol


def load_generation_protocol(path: Path) -> dict[str, Any]:
    try:
        protocol = json.loads(path.read_text(encoding="utf-8"))
        schema_path = path.with_name("competition_conformer_protocol_v1.schema.json")
        schema = json.loads(schema_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot load competition conformer protocol: {exc}") from exc
    if not isinstance(protocol, dict) or not isinstance(schema, dict):
        raise ValueError("Protocol and schema roots must be JSON objects")
    validate_json_schema(protocol, schema)
    return validate_generation_protocol(protocol)


def _resolve_under(root: Path, relative_path: str, label: str) -> Path:
    path = Path(relative_path)
    if path.is_absolute():
        raise ValueError(f"{label} must be relative to --data-root")
    resolved_root = root.resolve()
    resolved = (resolved_root / path).resolve()
    try:
        resolved.relative_to(resolved_root)
    except ValueError as exc:
        raise ValueError(f"{label} escapes --data-root") from exc
    if not resolved.is_file():
        raise FileNotFoundError(f"Missing {label}: {resolved}")
    return resolved


def load_authoritative_generation_inputs(
    protocol_path: Path, data_root: Path
) -> tuple[dict[str, Any], str, list[dict[str, Any]], str, dict[str, Any], str]:
    protocol = load_generation_protocol(protocol_path)
    protocol_sha256 = sha256_file(protocol_path)
    authority = protocol["candidate_authority"]
    manifest_path = _resolve_under(data_root, authority["path"], "candidate manifest")
    manifest_sha256 = sha256_file(manifest_path)
    if manifest_sha256 != authority["sha256"]:
        raise ValueError(
            f"Canonical candidate manifest SHA-256 mismatch: {manifest_sha256} != {authority['sha256']}"
        )
    manifest, candidates = load_candidate_manifest(manifest_path)
    if manifest["manifest_id"] != authority["manifest_id"]:
        raise ValueError("Canonical candidate manifest ID differs from the protocol contract")

    chemistry_ref = protocol["chemistry_contract"]
    chemistry_path = _resolve_under(data_root, chemistry_ref["path"], "chemistry contract")
    chemistry_sha256 = sha256_file(chemistry_path)
    if chemistry_sha256 != chemistry_ref["sha256"]:
        raise ValueError(
            f"Chemistry contract SHA-256 mismatch: {chemistry_sha256} != {chemistry_ref['sha256']}"
        )
    try:
        chemistry_contract = json.loads(chemistry_path.read_text(encoding="utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot read chemistry contract: {exc}") from exc
    if chemistry_contract.get("contract_id") != chemistry_ref["contract_id"]:
        raise ValueError("Chemistry contract ID differs from the protocol contract")
    if chemistry_contract.get("status") != "FROZEN":
        raise ValueError("Chemistry contract is not FROZEN")
    return (
        protocol,
        protocol_sha256,
        candidates,
        manifest_sha256,
        chemistry_contract,
        chemistry_sha256,
    )
