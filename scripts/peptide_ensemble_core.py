"""Core utilities for V0.3/V0.4 peptide-level pose ensemble analysis."""

from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path
from typing import Any

from pose_analysis_core import parse_metric, parse_semicolon_integers


ENSEMBLE_PATH_FIELDS = ("receptor_monomer_path", "receptor_dimer_path")
METRIC_FIELDS: dict[str, tuple[str, Any]] = {
    "model": ("model", int),
    "score": ("score", float),
    "contact_fraction": ("contact_fraction", float),
    "contacted_peptide_residues": (
        "contacted_peptide_residues",
        parse_semicolon_integers,
    ),
    "receptor_contact_residues": (
        "receptor_contact_residues",
        parse_semicolon_integers,
    ),
    "target_contacts": ("target_contacts", parse_semicolon_integers),
    "target_contact_count": ("target_contact_count", int),
    "receptor_clash_pairs": ("receptor_clash_pairs", int),
    "min_receptor_distance": ("min_receptor_distance", float),
    "chain_b_clash_pairs": ("chain_b_clash_pairs", int),
    "min_chain_b_distance": ("min_chain_b_distance", float),
    "basic_geometry_pass": ("basic_geometry_pass", int),
}
REQUIRED_IMPORTED_METRICS = (
    "score",
    "contact_fraction",
    "receptor_clash_pairs",
    "min_receptor_distance",
    "chain_b_clash_pairs",
    "min_chain_b_distance",
    "basic_geometry_pass",
)


def validate_ensemble_spec(spec: dict[str, Any]) -> None:
    schema_version = spec.get("schema_version")
    if schema_version not in {"0.3", "0.4"}:
        raise ValueError("schema_version must be '0.3' or '0.4'")
    for field in ("analysis_id", "candidate", "sequence", "group", "receptor_chain"):
        if not isinstance(spec.get(field), str) or not spec[field].strip():
            raise ValueError(f"{field} must be a non-empty string")
    if re.fullmatch(r"[A-Za-z0-9_]+", spec["receptor_chain"]) is None:
        raise ValueError("receptor_chain must be a safe PyMOL token")
    path_fields = list(ENSEMBLE_PATH_FIELDS)
    if schema_version == "0.4":
        path_fields.extend(
            (
                "docking_receptor_pdbqt_path",
                "receptor_registry_path",
                "chemistry_contract_path",
            )
        )
        _validate_v04_provenance(spec)
    for field in path_fields:
        _require_relative_path(spec.get(field), field)

    residues = spec.get("target_residues")
    if (
        not isinstance(residues, list)
        or not residues
        or any(isinstance(value, bool) or not isinstance(value, int) or value < 1 for value in residues)
        or len(residues) != len(set(residues))
    ):
        raise ValueError("target_residues must be unique positive integers")

    conformers = spec.get("conformers")
    if not isinstance(conformers, list) or not conformers:
        raise ValueError("conformers must be a non-empty list")
    names: list[str] = []
    expected_sum = 0
    for index, conformer in enumerate(conformers):
        prefix = f"conformers[{index}]"
        if not isinstance(conformer, dict):
            raise ValueError(f"{prefix} must be an object")
        name = conformer.get("conformer_name")
        if not isinstance(name, str) or not name.strip():
            raise ValueError(f"{prefix}.conformer_name must be a non-empty string")
        if re.fullmatch(r"[A-Za-z0-9_]+", name) is None:
            raise ValueError(f"{prefix}.conformer_name must be a safe token")
        names.append(name)
        conformer_path_fields = ["docking_output_path", "pose_metrics_path"]
        if schema_version == "0.4":
            conformer_path_fields.extend(
                (
                    "starting_sdf_path",
                    "ligand_pdbqt_path",
                    "vina_config_path",
                    "vina_log_path",
                )
            )
        for field in conformer_path_fields:
            _require_relative_path(conformer.get(field), f"{prefix}.{field}")
        count = conformer.get("expected_model_count")
        if isinstance(count, bool) or not isinstance(count, int) or count < 1:
            raise ValueError(f"{prefix}.expected_model_count must be a positive integer")
        expected_sum += count
    if len(names) != len(set(names)):
        raise ValueError("conformer_name values must be unique")

    expected_total = spec.get("expected_total_pose_count")
    if isinstance(expected_total, bool) or not isinstance(expected_total, int) or expected_total < 1:
        raise ValueError("expected_total_pose_count must be a positive integer")
    if expected_total != expected_sum:
        raise ValueError(
            f"expected_total_pose_count {expected_total} does not equal "
            f"the conformer expected_model_count sum {expected_sum}"
        )

    cutoffs = spec.get("contact_cutoffs")
    if not isinstance(cutoffs, dict) or cutoffs.get("direct_A") != 4.0 or cutoffs.get("near_A") != 5.0:
        raise ValueError("contact_cutoffs are fixed at direct_A=4.0 and near_A=5.0")

    rules = spec.get("representative_selection_rules")
    if not isinstance(rules, dict):
        raise ValueError("representative_selection_rules must be an object")
    allowed_rules = {
        "max_receptor_clash_pairs",
        "max_chain_B_clash_pairs",
        "min_direct_target_contacts",
        "required_any_direct_residues",
        "max_per_conformer",
        "diversity_key",
    }
    unsupported_rules = sorted(set(rules) - allowed_rules)
    if unsupported_rules:
        raise ValueError(
            f"unsupported representative selection rules: {unsupported_rules}"
        )
    expected_rule_values = {
        "max_receptor_clash_pairs": 0,
        "max_chain_B_clash_pairs": 0,
        "diversity_key": "exact_direct_contact_pattern",
    }
    for field, expected in expected_rule_values.items():
        if rules.get(field) != expected:
            raise ValueError(f"representative_selection_rules.{field} must be {expected!r}")
    for field in ("min_direct_target_contacts", "max_per_conformer"):
        value = rules.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"representative_selection_rules.{field} must be positive")
    if int(rules["max_per_conformer"]) > 3:
        raise ValueError("representative_selection_rules.max_per_conformer must be <= 3")
    required_residues = rules.get("required_any_direct_residues")
    if not isinstance(required_residues, list) or any(
        isinstance(value, bool) or not isinstance(value, int) or value < 1
        for value in required_residues
    ):
        raise ValueError("required_any_direct_residues must be a list of positive integers")

    summary_rules = spec.get("summary_rules")
    threshold = summary_rules.get("recurrent_contact_frequency_min") if isinstance(summary_rules, dict) else None
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not 0 <= threshold <= 1:
        raise ValueError("recurrent_contact_frequency_min must be between 0 and 1")

    input_hashes = spec.get("input_sha256", {})
    if not isinstance(input_hashes, dict):
        raise ValueError("input_sha256 must be an object when provided")
    allowed_hash_fields = {
        "receptor_monomer_path",
        "receptor_dimer_path",
        "conformers",
    }
    if schema_version == "0.4":
        allowed_hash_fields.add("docking_receptor_pdbqt_path")
    unknown_hash_fields = sorted(set(input_hashes) - allowed_hash_fields)
    if unknown_hash_fields:
        raise ValueError(f"input_sha256 contains unknown fields: {unknown_hash_fields}")
    for field in path_fields:
        if field in input_hashes:
            _require_sha256(input_hashes[field], f"input_sha256.{field}")
    conformer_hashes = input_hashes.get("conformers", {})
    if not isinstance(conformer_hashes, dict):
        raise ValueError("input_sha256.conformers must be an object")
    unknown_conformers = sorted(set(conformer_hashes) - set(names))
    if unknown_conformers:
        raise ValueError(
            f"input_sha256.conformers contains unknown conformers: {unknown_conformers}"
        )
    for conformer, hashes in conformer_hashes.items():
        if not isinstance(hashes, dict):
            raise ValueError(f"input_sha256.conformers.{conformer} must be an object")
        allowed_conformer_hashes = {"docking_output_path", "pose_metrics_path"}
        if schema_version == "0.4":
            allowed_conformer_hashes.update(
                {
                    "starting_sdf_path",
                    "ligand_pdbqt_path",
                    "vina_config_path",
                    "vina_log_path",
                }
            )
        unknown_fields = sorted(set(hashes) - allowed_conformer_hashes)
        if unknown_fields:
            raise ValueError(
                f"input_sha256.conformers.{conformer} contains unknown fields: "
                f"{unknown_fields}"
            )
        for field, value in hashes.items():
            _require_sha256(
                value, f"input_sha256.conformers.{conformer}.{field}"
            )


def _validate_v04_provenance(spec: dict[str, Any]) -> None:
    for field in ("template_id", "receptor_id", "comparison_protocol_id"):
        if not isinstance(spec.get(field), str) or not spec[field].strip():
            raise ValueError(f"{field} must be a non-empty string for schema 0.4")
    _require_sha256(
        spec.get("receptor_registry_sha256"), "receptor_registry_sha256"
    )
    _require_sha256(
        spec.get("chemistry_contract_sha256"), "chemistry_contract_sha256"
    )
    sequence = spec["sequence"]
    if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", sequence) is None:
        raise ValueError("sequence must contain standard uppercase amino-acid codes")
    peptide_length = spec.get("peptide_length")
    if isinstance(peptide_length, bool) or not isinstance(peptide_length, int) or peptide_length < 1:
        raise ValueError("peptide_length must be a positive integer")
    if peptide_length != len(sequence):
        raise ValueError(
            f"peptide_length {peptide_length} does not equal sequence length {len(sequence)}"
        )
    if spec.get("screening_protocol_status") not in {
        "formal_standardized",
        "legacy_regression_only",
    }:
        raise ValueError(
            "screening_protocol_status must be formal_standardized or "
            "legacy_regression_only"
        )
    if spec.get("vina_score_primary_source") != "docking_output_pdbqt_remark":
        raise ValueError(
            "vina_score_primary_source must be docking_output_pdbqt_remark"
        )
    tolerance = spec.get("score_crosscheck_tolerance_kcal_mol")
    if isinstance(tolerance, bool) or not isinstance(tolerance, (int, float)) or float(tolerance) != 0.001:
        raise ValueError("score_crosscheck_tolerance_kcal_mol must be 0.001")

    ligand = spec.get("ligand_preparation")
    if not isinstance(ligand, dict):
        raise ValueError("ligand_preparation must be an object")
    if not isinstance(ligand.get("protocol_id"), str) or not ligand["protocol_id"].strip():
        raise ValueError("ligand_preparation.protocol_id must be a non-empty string")
    formal_charge = ligand.get("formal_charge")
    if isinstance(formal_charge, bool) or not isinstance(formal_charge, int):
        raise ValueError("ligand_preparation.formal_charge must be an integer")

    docking = spec.get("docking_protocol")
    if not isinstance(docking, dict):
        raise ValueError("docking_protocol must be an object")
    required_text = {
        "protocol_id": None,
        "engine": "AutoDock Vina",
        "version": None,
        "scoring_function": "vina",
        "score_units": "kcal/mol",
        "score_direction": "lower_is_better",
    }
    for field, expected in required_text.items():
        value = docking.get(field)
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"docking_protocol.{field} must be a non-empty string")
        if expected is not None and value != expected:
            raise ValueError(f"docking_protocol.{field} must be {expected!r}")
    for field in ("exhaustiveness", "num_modes", "energy_range", "cpu"):
        value = docking.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"docking_protocol.{field} must be a positive integer")
    seed = docking.get("seed")
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("docking_protocol.seed must be an integer")
    for field in ("box_center_A", "box_size_A"):
        values = docking.get(field)
        if (
            not isinstance(values, list)
            or len(values) != 3
            or any(isinstance(value, bool) or not isinstance(value, (int, float)) for value in values)
        ):
            raise ValueError(f"docking_protocol.{field} must contain three numbers")
    if any(float(value) <= 0 for value in docking["box_size_A"]):
        raise ValueError("docking_protocol.box_size_A values must be positive")


def validate_legacy_contact_fraction(
    metrics: dict[str, Any], peptide_length: int, model: int, tolerance: float = 1e-9
) -> None:
    """Verify the legacy 5 A fraction against its explicit peptide denominator."""
    residues = metrics.get("contacted_peptide_residues")
    fraction = metrics.get("contact_fraction")
    if not isinstance(residues, list) or not isinstance(fraction, (int, float)):
        raise ValueError(f"MODEL {model} lacks legacy contact_fraction provenance")
    if any(
        isinstance(value, bool) or not isinstance(value, int) or value < 1 or value > peptide_length
        for value in residues
    ):
        raise ValueError(
            f"MODEL {model} contacted peptide residue is outside 1..{peptide_length}"
        )
    if len(residues) != len(set(residues)):
        raise ValueError(f"MODEL {model} contacted peptide residues contain duplicates")
    expected = len(residues) / peptide_length
    if abs(float(fraction) - expected) > tolerance:
        raise ValueError(
            f"MODEL {model} legacy contact_fraction {fraction} does not match "
            f"{len(residues)}/{peptide_length}={expected:.9f}"
        )


def _require_relative_path(value: Any, field: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field} must be a non-empty relative path")
    if Path(value).is_absolute():
        raise ValueError(f"{field} must be relative to --data-root")


def _require_sha256(value: Any, field: str) -> None:
    if not isinstance(value, str) or re.fullmatch(r"[0-9a-fA-F]{64}", value) is None:
        raise ValueError(f"{field} must be exactly 64 hexadecimal characters")


def validate_model_inventory(
    conformer_name: str,
    pdbqt_labels: set[int],
    metrics_labels: set[int],
    expected_count: int,
) -> None:
    actual_count = len(pdbqt_labels)
    if actual_count != expected_count:
        raise ValueError(
            f"{conformer_name} MODEL count mismatch: expected {expected_count}, "
            f"actual {actual_count}"
        )
    if pdbqt_labels != metrics_labels:
        only_pdbqt = sorted(pdbqt_labels - metrics_labels)
        only_metrics = sorted(metrics_labels - pdbqt_labels)
        raise ValueError(
            f"{conformer_name} MODEL labels differ: only in PDBQT {only_pdbqt}; "
            f"only in metrics {only_metrics}"
        )


def read_pose_metrics_table(path: Path) -> tuple[list[str], dict[int, dict[str, Any]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames
        if not header or "model" not in header:
            raise ValueError(f"Pose metrics header has no model field: {path}")
        raw_rows = list(reader)
    fields = set(header)
    result: dict[int, dict[str, Any]] = {}
    for raw_row in raw_rows:
        model = parse_metric(raw_row, fields, "model", int)
        if not isinstance(model, int):
            raise ValueError(f"Pose metrics contains invalid model value: {raw_row.get('model')!r}")
        if model in result:
            raise ValueError(f"Duplicate MODEL {model} in pose metrics: {path}")
        parsed = {
            output_name: parse_metric(raw_row, fields, source_name, converter)
            for output_name, (source_name, converter) in METRIC_FIELDS.items()
        }
        missing = [field for field in REQUIRED_IMPORTED_METRICS if parsed[field] == "NA"]
        if missing:
            raise ValueError(f"MODEL {model} is missing required pose metrics: {missing}")
        parsed["raw_model_row"] = {name: raw_row.get(name, "") for name in header}
        result[model] = parsed
    if not result:
        raise ValueError(f"Pose metrics contains no rows: {path}")
    return list(header), result


def annotate_filter_provenance(
    record: dict[str, Any], rules: dict[str, Any]
) -> dict[str, Any]:
    result = dict(record)
    direct_contacts = set(int(value) for value in record["direct_target_contacts"])
    receptor_clash_free = int(record["receptor_clash_pairs"]) <= int(
        rules["max_receptor_clash_pairs"]
    )
    chain_b_clash_free = int(record["chain_B_clash_pairs"]) <= int(
        rules["max_chain_B_clash_pairs"]
    )
    enough_direct = len(direct_contacts) >= int(rules["min_direct_target_contacts"])
    required_any = set(int(value) for value in rules["required_any_direct_residues"])
    required_hit = not required_any or bool(direct_contacts & required_any)
    interface_eligible = enough_direct and required_hit
    representative_eligible = (
        receptor_clash_free
        and chain_b_clash_free
        and interface_eligible
    )
    reasons: list[str] = []
    if not receptor_clash_free:
        reasons.append("receptor_clash_pairs>0")
    if not chain_b_clash_free:
        reasons.append("chain_B_clash_pairs>0")
    if not enough_direct:
        reasons.append(
            f"direct_target_contact_count<{int(rules['min_direct_target_contacts'])}"
        )
    if not required_hit:
        reasons.append("required_any_direct_residues_not_hit")

    if not receptor_clash_free or not chain_b_clash_free:
        status = "filtered_clash"
    elif not interface_eligible:
        status = "filtered_low_interface"
    else:
        status = "kept"
    result.update(
        {
            "receptor_clash_free": receptor_clash_free,
            "chain_B_clash_free": chain_b_clash_free,
            "interface_eligible": interface_eligible,
            "representative_eligible": representative_eligible,
            "filter_reasons": reasons,
            "selection_status": status,
        }
    )
    return result


def contact_pattern(values: list[int]) -> tuple[int, ...]:
    return tuple(sorted(set(int(value) for value in values)))


def select_contact_pattern_representatives(
    records: list[dict[str, Any]], max_per_conformer: int
) -> list[dict[str, Any]]:
    result = [dict(record) for record in records]
    conformers = list(dict.fromkeys(str(record["conformer"]) for record in result))
    for conformer in conformers:
        conformer_records = [record for record in result if record["conformer"] == conformer]
        eligible = [record for record in conformer_records if record["representative_eligible"]]
        all_counts = Counter(contact_pattern(record["direct_target_contacts"]) for record in conformer_records)
        eligible_counts = Counter(contact_pattern(record["direct_target_contacts"]) for record in eligible)
        for record in conformer_records:
            pattern = contact_pattern(record["direct_target_contacts"])
            record["direct_contact_pattern"] = ";".join(map(str, pattern))
            record["all_pose_pattern_count"] = all_counts[pattern]
            record["all_pose_pattern_prevalence"] = all_counts[pattern] / len(conformer_records)
            record["eligible_pattern_count"] = eligible_counts[pattern]
            record["eligible_pattern_prevalence"] = (
                eligible_counts[pattern] / len(eligible) if eligible else 0.0
            )
            record["representative_rank_within_conformer"] = None

        winners: list[dict[str, Any]] = []
        for pattern in eligible_counts:
            candidates = [
                record
                for record in eligible
                if contact_pattern(record["direct_target_contacts"]) == pattern
            ]
            candidates.sort(
                key=lambda record: (
                    -float(record["legacy_5A_receptor_contact_fraction"]),
                    float(record["vina_score"]),
                    -float(record["minimum_receptor_distance_A"]),
                    int(record["model"]),
                )
            )
            winners.append(candidates[0])
        winners.sort(
            key=lambda record: (
                -len(contact_pattern(record["direct_target_contacts"])),
                -int(record["eligible_pattern_count"]),
                -float(record["legacy_5A_receptor_contact_fraction"]),
                float(record["vina_score"]),
                -float(record["minimum_receptor_distance_A"]),
                int(record["model"]),
            )
        )
        for rank, winner in enumerate(winners[:max_per_conformer], start=1):
            winner["selection_status"] = "representative"
            winner["representative_rank_within_conformer"] = rank
    return result


def build_contact_frequency_rows(
    records: list[dict[str, Any]],
    target_residues: list[int],
    conformer_names: list[str],
    recurrent_threshold: float,
) -> list[dict[str, Any]]:
    by_conformer = {
        name: [record for record in records if record["conformer"] == name]
        for name in conformer_names
    }
    if any(not rows for rows in by_conformer.values()):
        missing = [name for name, rows in by_conformer.items() if not rows]
        raise ValueError(f"No pose records for conformers: {missing}")

    residue_names: dict[int, str] = {}
    for record in records:
        for detail in record.get("target_contact_details", []):
            residue_names[int(detail["target_residue_number"])] = str(
                detail["target_residue"]
            )

    output: list[dict[str, Any]] = []
    for residue in target_residues:
        conformer_direct_frequencies: list[float] = []
        conformer_near_frequencies: list[float] = []
        conformer_presence_count = 0
        recurrent_conformer_support_count = 0
        for name in conformer_names:
            rows = by_conformer[name]
            direct_count = sum(residue in record["direct_target_contacts"] for record in rows)
            near_count = sum(residue in record["near_target_contacts"] for record in rows)
            direct_frequency = direct_count / len(rows)
            near_frequency = near_count / len(rows)
            conformer_direct_frequencies.append(direct_frequency)
            conformer_near_frequencies.append(near_frequency)
            conformer_presence_count += int(direct_count > 0)
            recurrent_conformer_support_count += int(
                direct_frequency >= recurrent_threshold
            )
            output.append(
                _frequency_row(
                    "conformer",
                    name,
                    residue,
                    residue_names.get(residue, str(residue)),
                    len(rows),
                    direct_count,
                    near_count,
                )
            )

        direct_count = sum(residue in record["direct_target_contacts"] for record in records)
        near_count = sum(residue in record["near_target_contacts"] for record in records)
        peptide_row = _frequency_row(
            "peptide",
            "ALL",
            residue,
            residue_names.get(residue, str(residue)),
            len(records),
            direct_count,
            near_count,
        )
        peptide_row.update(
            {
                "conformer_balanced_direct_frequency": sum(conformer_direct_frequencies)
                / len(conformer_names),
                "conformer_balanced_near_frequency": sum(conformer_near_frequencies)
                / len(conformer_names),
                "conformer_presence_count": conformer_presence_count,
                "recurrent_conformer_support_count": (
                    recurrent_conformer_support_count
                ),
            }
        )
        output.insert(len(output) - len(conformer_names), peptide_row)
    return output


def _frequency_row(
    scope: str,
    conformer: str,
    residue: int,
    residue_name: str,
    pose_count: int,
    direct_count: int,
    near_count: int,
) -> dict[str, Any]:
    return {
        "scope": scope,
        "conformer": conformer,
        "target_residue_number": residue,
        "target_residue": residue_name,
        "pose_count": pose_count,
        "direct_count": direct_count,
        "direct_frequency": direct_count / pose_count,
        "near_count": near_count,
        "near_frequency": near_count / pose_count,
        "direct_or_near_count": direct_count + near_count,
        "direct_or_near_frequency": (direct_count + near_count) / pose_count,
        "conformer_balanced_direct_frequency": None,
        "conformer_balanced_near_frequency": None,
        "conformer_presence_count": None,
        "recurrent_conformer_support_count": None,
    }


def shared_recurrent_direct_core(
    peptide_frequency_rows: list[dict[str, Any]], number_of_conformers: int
) -> list[int]:
    return sorted(
        int(row["target_residue_number"])
        for row in peptide_frequency_rows
        if int(row["recurrent_conformer_support_count"]) == number_of_conformers
    )
