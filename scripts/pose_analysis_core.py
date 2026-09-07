"""Reusable core utilities for V0.2 parameterized single-pose analysis."""

from __future__ import annotations

import csv
import hashlib
import re
from dataclasses import dataclass
from math import dist
from pathlib import Path
from typing import Any, Callable


INPUT_PATH_FIELDS = (
    "docking_output_path",
    "pose_metrics_path",
    "receptor_monomer_path",
    "receptor_dimer_path",
)

REQUIRED_TEXT_FIELDS = (
    "analysis_id",
    "candidate",
    "sequence",
    "group",
    "conformer",
    "receptor_chain",
)


@dataclass(frozen=True)
class CoordinateAtom:
    serial: int
    name: str
    atom_type: str
    xyz: tuple[float, float, float]
    residue_number: int | None
    residue_name: str
    chain: str

    @property
    def is_heavy(self) -> bool:
        return not self.atom_type.upper().startswith("H")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_spec(spec: dict[str, Any]) -> None:
    if spec.get("schema_version") != "0.2":
        raise ValueError("schema_version must be '0.2'")
    for field in REQUIRED_TEXT_FIELDS:
        if not isinstance(spec.get(field), str) or not spec[field].strip():
            raise ValueError(f"{field} must be a non-empty string")
    _require_safe_pymol_token(spec["receptor_chain"], "receptor_chain")
    model = spec.get("model")
    if isinstance(model, bool) or not isinstance(model, int) or model < 1:
        raise ValueError("model must be a positive integer Vina MODEL label")
    for field in INPUT_PATH_FIELDS:
        raw_path = spec.get(field)
        if not isinstance(raw_path, str) or not raw_path.strip():
            raise ValueError(f"{field} must be a non-empty path string")
        if Path(raw_path).is_absolute():
            raise ValueError(f"{field} must be relative to --data-root: {raw_path}")
    residues = spec.get("target_residues")
    if (
        not isinstance(residues, list)
        or not residues
        or any(
            isinstance(value, bool) or not isinstance(value, int) or value < 1
            for value in residues
        )
        or len(set(residues)) != len(residues)
    ):
        raise ValueError(
            "target_residues must be a non-empty list of unique positive integers"
        )
    direct_cutoff = spec.get("direct_contact_cutoff_A")
    near_cutoff = spec.get("near_contact_cutoff_A")
    for field, cutoff in (
        ("direct_contact_cutoff_A", direct_cutoff),
        ("near_contact_cutoff_A", near_cutoff),
    ):
        if (
            isinstance(cutoff, bool)
            or not isinstance(cutoff, (int, float))
            or cutoff <= 0
        ):
            raise ValueError(f"{field} must be a positive number")
    if float(direct_cutoff) != 4.0 or float(near_cutoff) != 5.0:
        raise ValueError(
            "V0.2 contact cutoffs are fixed at direct_contact_cutoff_A=4.0 "
            "and near_contact_cutoff_A=5.0"
        )
    interactions = spec.get("manual_interactions")
    if not isinstance(interactions, list):
        raise ValueError("manual_interactions must be a list")
    for index, interaction in enumerate(interactions):
        _validate_manual_interaction(interaction, index)

    expected_hashes = spec.get("input_sha256", {})
    if not isinstance(expected_hashes, dict):
        raise ValueError("input_sha256 must be an object when provided")
    unknown_hashes = sorted(set(expected_hashes) - set(INPUT_PATH_FIELDS))
    if unknown_hashes:
        raise ValueError(f"input_sha256 contains unknown fields: {unknown_hashes}")
    for field, value in expected_hashes.items():
        if (
            not isinstance(value, str)
            or re.fullmatch(r"[0-9a-fA-F]{64}", value) is None
        ):
            raise ValueError(
                f"input_sha256.{field} must be exactly 64 hexadecimal characters"
            )


def _validate_manual_interaction(interaction: Any, index: int) -> None:
    prefix = f"manual_interactions[{index}]"
    if not isinstance(interaction, dict):
        raise ValueError(f"{prefix} must be an object")
    if not isinstance(interaction.get("interaction_id"), str) or not interaction[
        "interaction_id"
    ].strip():
        raise ValueError(f"{prefix}.interaction_id must be a non-empty string")
    distance = interaction.get("distance")
    if not isinstance(distance, dict):
        raise ValueError(f"{prefix}.distance must be an object")
    _require_positive_integer(distance, "ligand_atom_id", f"{prefix}.distance")
    _require_positive_integer(distance, "receptor_residue", f"{prefix}.distance")
    if not isinstance(distance.get("receptor_atom"), str) or not distance[
        "receptor_atom"
    ].strip():
        raise ValueError(f"{prefix}.distance.receptor_atom must be a non-empty string")
    _require_safe_pymol_token(
        distance["receptor_atom"], f"{prefix}.distance.receptor_atom"
    )
    angle = interaction.get("angle")
    if angle is not None:
        if not isinstance(angle, dict):
            raise ValueError(f"{prefix}.angle must be an object or null")
        _require_positive_integer(angle, "ligand_donor_atom_id", f"{prefix}.angle")
        _require_positive_integer(angle, "ligand_hydrogen_atom_id", f"{prefix}.angle")
        acceptor = angle.get("receptor_acceptor")
        if not isinstance(acceptor, dict):
            raise ValueError(f"{prefix}.angle.receptor_acceptor must be an object")
        _require_positive_integer(
            acceptor, "residue", f"{prefix}.angle.receptor_acceptor"
        )
        if not isinstance(acceptor.get("atom"), str) or not acceptor["atom"].strip():
            raise ValueError(
                f"{prefix}.angle.receptor_acceptor.atom must be a non-empty string"
            )
        _require_safe_pymol_token(
            acceptor["atom"], f"{prefix}.angle.receptor_acceptor.atom"
        )


def _require_positive_integer(mapping: dict[str, Any], field: str, prefix: str) -> None:
    value = mapping.get(field)
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{prefix}.{field} must be a positive integer")


def _require_safe_pymol_token(value: str, field: str) -> None:
    if re.fullmatch(r"[A-Za-z0-9_]+", value) is None:
        raise ValueError(f"{field} must be a safe PyMOL token: {value!r}")


def parse_vina_model_blocks(path: Path) -> dict[int, str]:
    """Return complete Vina MODEL blocks keyed by their integer labels."""
    text = path.read_text(encoding="utf-8")
    blocks = re.findall(r"(?ms)^MODEL\s+\d+\s*$.*?^ENDMDL\s*$", text)
    result: dict[int, str] = {}
    for block in blocks:
        match = re.match(r"MODEL\s+(\d+)\s*$", block.splitlines()[0])
        if match is None:
            raise ValueError(f"Cannot parse MODEL label in {path}")
        label = int(match.group(1))
        if label in result:
            raise ValueError(f"Duplicate MODEL label {label} in Vina output: {path}")
        result[label] = block + "\n"
    if not result:
        raise ValueError(f"No MODEL/ENDMDL blocks found in Vina output: {path}")
    return result


def require_model_block(blocks: dict[int, str], model_label: int) -> str:
    try:
        return blocks[model_label]
    except KeyError as exc:
        available = ", ".join(str(label) for label in blocks)
        raise ValueError(
            f"MODEL {model_label} is absent from Vina output; "
            f"available MODEL labels: {available}"
        ) from exc


def parse_pdbqt_atoms(model_block: str) -> list[CoordinateAtom]:
    atoms: list[CoordinateAtom] = []
    for line in model_block.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        fields = line.split()
        atoms.append(
            CoordinateAtom(
                serial=int(line[6:11]),
                name=line[12:16].strip(),
                atom_type=fields[-1],
                xyz=(
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ),
                residue_number=(
                    int(line[22:26]) if line[22:26].strip() else None
                ),
                residue_name=line[17:20].strip(),
                chain=line[21:22].strip(),
            )
        )
    if not atoms:
        raise ValueError("Selected MODEL block contains no ligand atoms")
    return atoms


def parse_receptor_pdb_atoms(path: Path, chain: str) -> list[CoordinateAtom]:
    atoms: list[CoordinateAtom] = []
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom_chain = line[21:22].strip()
        if atom_chain != chain:
            continue
        element = line[76:78].strip()
        if not element:
            element = next(
                (character for character in line[12:16] if character.isalpha()), ""
            )
        atoms.append(
            CoordinateAtom(
                serial=int(line[6:11]),
                name=line[12:16].strip(),
                atom_type=element,
                xyz=(
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ),
                residue_number=int(line[22:26]),
                residue_name=line[17:20].strip(),
                chain=atom_chain,
            )
        )
    if not atoms:
        raise ValueError(f"Receptor chain {chain!r} contains no atoms in {path}")
    return atoms


def compute_target_residue_contacts(
    ligand_atoms: list[CoordinateAtom],
    receptor_atoms: list[CoordinateAtom],
    target_residues: list[int],
    direct_cutoff_A: float,
    near_cutoff_A: float,
) -> list[dict[str, Any]]:
    ligand_heavy = [atom for atom in ligand_atoms if atom.is_heavy]
    receptor_heavy = [atom for atom in receptor_atoms if atom.is_heavy]
    if not ligand_heavy:
        raise ValueError("Selected ligand MODEL contains no heavy atoms")

    by_residue: dict[int, list[CoordinateAtom]] = {}
    for atom in receptor_heavy:
        if atom.residue_number is not None:
            by_residue.setdefault(atom.residue_number, []).append(atom)

    rows: list[dict[str, Any]] = []
    for residue_number in target_residues:
        residue_atoms = by_residue.get(residue_number, [])
        if not residue_atoms:
            raise ValueError(
                f"Target residue {residue_number} has no heavy atoms in receptor chain"
            )
        minimum, ligand_atom, receptor_atom = min(
            (
                (dist(ligand.xyz, receptor.xyz), ligand, receptor)
                for ligand in ligand_heavy
                for receptor in residue_atoms
            ),
            key=lambda item: item[0],
        )
        rows.append(
            {
                "target_residue_number": residue_number,
                "target_residue": f"{receptor_atom.residue_name}{residue_number}",
                "classification": classify_contact(
                    minimum, direct_cutoff_A, near_cutoff_A
                ),
                "min_distance_A": minimum,
                "ligand_atom_id": ligand_atom.serial,
                "ligand_atom_name": ligand_atom.name,
                "ligand_atom_type": ligand_atom.atom_type,
                "receptor_atom": receptor_atom.name,
                "direct_cutoff_A": direct_cutoff_A,
                "near_cutoff_A": near_cutoff_A,
                "provenance": "v02_raw_coordinate_heavy_atom_recomputed",
            }
        )
    return rows


def resolve_inputs(
    spec: dict[str, Any], data_root: Path
) -> tuple[dict[str, Path], dict[str, str]]:
    root = data_root.resolve()
    inputs: dict[str, Path] = {}
    hashes: dict[str, str] = {}
    expected_hashes = spec.get("input_sha256", {})
    for field in INPUT_PATH_FIELDS:
        raw_path = Path(str(spec[field]))
        path = raw_path.resolve() if raw_path.is_absolute() else (root / raw_path).resolve()
        if not raw_path.is_absolute() and not path.is_relative_to(root):
            raise ValueError(f"{field} resolves outside --data-root: {raw_path}")
        if not path.is_file():
            raise FileNotFoundError(f"Required input is missing ({field}): {path}")
        actual_hash = sha256(path)
        expected_hash = expected_hashes.get(field)
        if expected_hash is not None and actual_hash.lower() != str(expected_hash).lower():
            raise ValueError(
                f"SHA-256 mismatch for {field}: expected {expected_hash}, "
                f"actual {actual_hash}, path {path}"
            )
        inputs[field] = path
        hashes[field] = actual_hash
    return inputs, hashes


def parse_metric(
    row: dict[str, str], header: set[str], field: str, converter: Callable[[str], Any]
) -> Any:
    if field not in header:
        return "NA"
    raw = row.get(field, "").strip()
    if raw == "":
        return "NA"
    try:
        return converter(raw)
    except (TypeError, ValueError):
        return "NA"


def parse_semicolon_integers(raw: str) -> list[int]:
    if not raw.strip():
        return []
    return [int(value) for value in raw.split(";") if value.strip()]


def read_model_metrics(path: Path, model_label: int) -> tuple[list[str], dict[str, Any]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames
        if not header or len(header) < 2:
            raise ValueError(f"Pose metrics is not a valid TSV: {path}")
        if "model" not in header:
            raise ValueError(f"Pose metrics header has no 'model' field: {header}")
        model_row = next(
            (row for row in reader if row.get("model", "").strip() == str(model_label)),
            None,
        )
    if model_row is None:
        raise ValueError(f"MODEL {model_label} is absent from pose metrics: {path}")

    fields = set(header)
    metrics = {
        "model": parse_metric(model_row, fields, "model", int),
        "score": parse_metric(model_row, fields, "score", float),
        "contact_fraction": parse_metric(model_row, fields, "contact_fraction", float),
        "contacted_peptide_residues": parse_metric(
            model_row, fields, "contacted_peptide_residues", parse_semicolon_integers
        ),
        "receptor_contact_residues": parse_metric(
            model_row, fields, "receptor_contact_residues", parse_semicolon_integers
        ),
        "target_contacts": parse_metric(
            model_row, fields, "target_contacts", parse_semicolon_integers
        ),
        "target_contact_count": parse_metric(
            model_row, fields, "target_contact_count", int
        ),
        "receptor_clash_pairs": parse_metric(
            model_row, fields, "receptor_clash_pairs", int
        ),
        "min_receptor_distance": parse_metric(
            model_row, fields, "min_receptor_distance", float
        ),
        "chain_b_clash_pairs": parse_metric(
            model_row, fields, "chain_b_clash_pairs", int
        ),
        "min_chain_b_distance": parse_metric(
            model_row, fields, "min_chain_b_distance", float
        ),
        "basic_geometry_pass": parse_metric(
            model_row, fields, "basic_geometry_pass", int
        ),
        "raw_model_row": {name: model_row.get(name, "") for name in header},
    }
    return list(header), metrics


def compare_target_contacts(
    recomputed: list[int], pose_metrics: list[int] | str
) -> dict[str, Any]:
    recomputed_set = set(recomputed)
    if not isinstance(pose_metrics, list):
        return {
            "status": "NOT_AVAILABLE",
            "v02_recomputed": sorted(recomputed_set),
            "pose_metrics": None,
            "only_in_v02_recomputed": [],
            "only_in_pose_metrics": [],
        }
    pose_metrics_set = set(pose_metrics)
    only_recomputed = sorted(recomputed_set - pose_metrics_set)
    only_metrics = sorted(pose_metrics_set - recomputed_set)
    return {
        "status": "MATCH" if not only_recomputed and not only_metrics else "WARNING_MISMATCH",
        "v02_recomputed": sorted(recomputed_set),
        "pose_metrics": sorted(pose_metrics_set),
        "only_in_v02_recomputed": only_recomputed,
        "only_in_pose_metrics": only_metrics,
    }


def classify_contact(
    distance_A: float, direct_cutoff_A: float, near_cutoff_A: float
) -> str:
    """Classify an unrounded heavy-atom minimum distance."""
    if direct_cutoff_A <= 0 or near_cutoff_A <= direct_cutoff_A:
        raise ValueError("Contact cutoffs must satisfy 0 < direct < near")
    if distance_A <= direct_cutoff_A:
        return "direct"
    if distance_A <= near_cutoff_A:
        return "near"
    return "none"


def compare_legacy_contacts(
    direct_contacts: list[int],
    near_contacts: list[int],
    legacy_5A_contacts: list[int] | str,
) -> dict[str, Any]:
    """Compare formal V0.2 contacts with the legacy 5 Å set."""
    direct = set(direct_contacts)
    near = set(near_contacts)
    if not isinstance(legacy_5A_contacts, list):
        return {
            "status": "NOT_AVAILABLE",
            "v02_direct_contacts": sorted(direct),
            "v02_near_contacts": sorted(near),
            "legacy_5A_target_contacts": None,
            "legacy_only_explained_by_near": [],
            "unexplained_residues": [],
        }
    legacy = set(legacy_5A_contacts)
    expected_legacy = direct | near
    unexplained = sorted(legacy ^ expected_legacy)
    explained = sorted((legacy - direct) & near)
    if unexplained:
        status = "WARNING_MISMATCH"
    elif explained:
        status = "EXPECTED_DEFINITION_DIFFERENCE"
    else:
        status = "MATCH"
    return {
        "status": status,
        "v02_direct_contacts": sorted(direct),
        "v02_near_contacts": sorted(near),
        "legacy_5A_target_contacts": sorted(legacy),
        "legacy_only_explained_by_near": explained,
        "unexplained_residues": unexplained,
    }
