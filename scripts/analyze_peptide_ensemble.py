"""CLI for V0.3 peptide-level multi-conformer pose ensemble analysis."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from statistics import mean, median
from typing import Any

from peptide_ensemble_core import (
    annotate_filter_provenance,
    build_contact_frequency_rows,
    read_pose_metrics_table,
    select_contact_pattern_representatives,
    shared_recurrent_direct_core,
    validate_ensemble_spec,
    validate_model_inventory,
)
from pose_analysis_core import (
    compare_legacy_contacts,
    compute_target_residue_contacts,
    parse_pdbqt_atoms,
    parse_receptor_pdb_atoms,
    parse_vina_model_blocks,
    sha256,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
TEST_OUTPUT_ROOT = (REPOSITORY_ROOT / "test_output").resolve()
SCRIPT_PATH = Path(__file__).resolve()
ENSEMBLE_CORE_PATH = SCRIPT_PATH.with_name("peptide_ensemble_core.py")
POSE_CORE_PATH = SCRIPT_PATH.with_name("pose_analysis_core.py")
OUTPUT_NAMES = (
    "pose_table.tsv",
    "pose_target_contacts.tsv",
    "target_contact_frequency.tsv",
    "conformer_summary.tsv",
    "representative_poses.tsv",
    "peptide_summary.md",
    "peptide_summary.json",
    "run_manifest.txt",
)
POSE_COLUMNS = (
    "pose_id",
    "conformer",
    "model",
    "vina_score",
    "legacy_5A_receptor_contact_fraction",
    "direct_target_contacts",
    "direct_target_contact_count",
    "near_target_contacts",
    "near_target_contact_count",
    "legacy_5A_target_contacts",
    "legacy_consistency_status",
    "receptor_clash_pairs",
    "minimum_receptor_distance_A",
    "chain_B_clash_pairs",
    "minimum_chain_B_distance_A",
    "legacy_basic_geometry_pass",
    "legacy_basic_geometry_pass_source",
    "receptor_clash_free",
    "chain_B_clash_free",
    "interface_eligible",
    "representative_eligible",
    "filter_reasons",
    "direct_contact_pattern",
    "all_pose_pattern_count",
    "all_pose_pattern_prevalence",
    "eligible_pattern_count",
    "eligible_pattern_prevalence",
    "selection_status",
    "representative_rank_within_conformer",
)
CONTACT_COLUMNS = (
    "pose_id",
    "conformer",
    "model",
    "target_residue_number",
    "target_residue",
    "classification",
    "min_distance_A",
    "ligand_atom_id",
    "ligand_atom_name",
    "ligand_atom_type",
    "receptor_atom",
    "direct_cutoff_A",
    "near_cutoff_A",
    "provenance",
)
FREQUENCY_COLUMNS = (
    "scope",
    "conformer",
    "target_residue_number",
    "target_residue",
    "pose_count",
    "direct_count",
    "direct_frequency",
    "near_count",
    "near_frequency",
    "direct_or_near_count",
    "direct_or_near_frequency",
    "conformer_balanced_direct_frequency",
    "conformer_balanced_near_frequency",
    "conformer_presence_count",
    "recurrent_conformer_support_count",
)
CONFORMER_COLUMNS = (
    "conformer",
    "total_poses",
    "receptor_clash_free_poses",
    "chain_B_clash_free_poses",
    "interface_eligible_poses",
    "representative_eligible_poses",
    "kept_poses",
    "filtered_clash_poses",
    "filtered_low_interface_poses",
    "representative_poses",
    "best_vina_score",
    "median_vina_score",
    "worst_vina_score",
    "median_legacy_5A_receptor_contact_fraction",
    "mean_direct_target_contact_count",
    "mean_near_target_contact_count",
    "recurrent_direct_target_residues",
    "recurrent_near_target_residues",
    "dominant_direct_contact_patterns",
    "representative_pose_ids",
    "legacy_warning_mismatch_count",
)


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def require_safe_output_directory(path: Path) -> Path:
    resolved = path.resolve()
    if not resolved.is_relative_to(TEST_OUTPUT_ROOT):
        raise ValueError(f"--output-dir must be inside {TEST_OUTPUT_ROOT}: {resolved}")
    if resolved.exists():
        raise FileExistsError(f"Output directory already exists: {resolved}")
    return resolved


def resolve_input(
    data_root: Path,
    raw_path: str,
    label: str,
    expected_hash: str | None,
) -> tuple[Path, dict[str, Any]]:
    root = data_root.resolve()
    candidate = Path(raw_path)
    if candidate.is_absolute():
        raise ValueError(f"{label} must be relative to --data-root")
    resolved = (root / candidate).resolve()
    if not resolved.is_relative_to(root):
        raise ValueError(f"{label} resolves outside --data-root: {raw_path}")
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} is missing: {resolved}")
    actual_hash = sha256(resolved)
    if expected_hash is not None and actual_hash.lower() != expected_hash.lower():
        raise ValueError(
            f"SHA-256 mismatch for {label}: expected {expected_hash}, actual {actual_hash}"
        )
    return resolved, {
        "spec_path": raw_path,
        "resolved_path": str(resolved),
        "sha256": actual_hash,
        "expected_sha256": expected_hash,
        "hash_verification": "MATCH" if expected_hash else "NOT_PROVIDED",
    }


def analyze_peptide_ensemble(
    spec_path: Path, data_root: Path, output_dir: Path
) -> Path:
    resolved_spec = spec_path.resolve()
    spec = load_json(resolved_spec)
    validate_ensemble_spec(spec)
    safe_output = require_safe_output_directory(output_dir)
    expected_hashes = spec.get("input_sha256", {})

    receptor_path, receptor_detail = resolve_input(
        data_root,
        spec["receptor_monomer_path"],
        "receptor_monomer_path",
        expected_hashes.get("receptor_monomer_path"),
    )
    dimer_path, dimer_detail = resolve_input(
        data_root,
        spec["receptor_dimer_path"],
        "receptor_dimer_path",
        expected_hashes.get("receptor_dimer_path"),
    )
    receptor_atoms = parse_receptor_pdb_atoms(receptor_path, spec["receptor_chain"])

    input_details: dict[str, Any] = {
        "receptor_monomer_path": receptor_detail,
        "receptor_dimer_path": {
            **dimer_detail,
            "role": "validation_hash_provenance_only_not_recalculated",
        },
        "conformers": {},
    }
    records: list[dict[str, Any]] = []
    contact_rows: list[dict[str, Any]] = []
    conformer_names = [item["conformer_name"] for item in spec["conformers"]]
    cutoff = spec["contact_cutoffs"]

    for conformer_spec in spec["conformers"]:
        conformer = conformer_spec["conformer_name"]
        conformer_hashes = expected_hashes.get("conformers", {}).get(conformer, {})
        docking_path, docking_detail = resolve_input(
            data_root,
            conformer_spec["docking_output_path"],
            f"{conformer}.docking_output_path",
            conformer_hashes.get("docking_output_path"),
        )
        metrics_path, metrics_detail = resolve_input(
            data_root,
            conformer_spec["pose_metrics_path"],
            f"{conformer}.pose_metrics_path",
            conformer_hashes.get("pose_metrics_path"),
        )
        blocks = parse_vina_model_blocks(docking_path)
        metrics_header, metrics_by_model = read_pose_metrics_table(metrics_path)
        validate_model_inventory(
            conformer,
            set(blocks),
            set(metrics_by_model),
            int(conformer_spec["expected_model_count"]),
        )
        input_details["conformers"][conformer] = {
            "docking_output_path": docking_detail,
            "pose_metrics_path": metrics_detail,
            "metrics_header": metrics_header,
            "model_labels": list(blocks),
        }

        for model, block in blocks.items():
            pose_id = f"{spec['candidate']}_{conformer}_m{model}"
            ligand_atoms = parse_pdbqt_atoms(block)
            details = compute_target_residue_contacts(
                ligand_atoms,
                receptor_atoms,
                spec["target_residues"],
                float(cutoff["direct_A"]),
                float(cutoff["near_A"]),
            )
            direct_contacts = [
                row["target_residue_number"]
                for row in details
                if row["classification"] == "direct"
            ]
            near_contacts = [
                row["target_residue_number"]
                for row in details
                if row["classification"] == "near"
            ]
            metrics = metrics_by_model[model]
            legacy_consistency = compare_legacy_contacts(
                direct_contacts, near_contacts, metrics["target_contacts"]
            )
            record = {
                "pose_id": pose_id,
                "conformer": conformer,
                "model": model,
                "vina_score": metrics["score"],
                "legacy_5A_receptor_contact_fraction": metrics["contact_fraction"],
                "direct_target_contacts": direct_contacts,
                "direct_target_contact_count": len(direct_contacts),
                "near_target_contacts": near_contacts,
                "near_target_contact_count": len(near_contacts),
                "legacy_5A_target_contacts": metrics["target_contacts"],
                "legacy_consistency_status": legacy_consistency["status"],
                "legacy_consistency": legacy_consistency,
                "receptor_clash_pairs": metrics["receptor_clash_pairs"],
                "minimum_receptor_distance_A": metrics["min_receptor_distance"],
                "chain_B_clash_pairs": metrics["chain_b_clash_pairs"],
                "minimum_chain_B_distance_A": metrics["min_chain_b_distance"],
                "legacy_basic_geometry_pass": metrics["basic_geometry_pass"],
                "legacy_basic_geometry_pass_source": "legacy_pose_metrics",
                "target_contact_details": details,
                "provenance": {
                    "direct_near_contacts": "V0.2 heavy-atom coordinate primitive",
                    "score_legacy_5A_receptor_contact_fraction_clashes_legacy_basic": (
                        "legacy_pose_metrics"
                    ),
                    "chain_B_metrics": "existing pose_metrics.tsv; dimer not recalculated",
                },
            }
            records.append(
                annotate_filter_provenance(
                    record, spec["representative_selection_rules"]
                )
            )
            contact_rows.extend(
                {"pose_id": pose_id, "conformer": conformer, "model": model, **row}
                for row in details
            )

    if len(records) != int(spec["expected_total_pose_count"]):
        raise ValueError(
            f"Total pose count mismatch: expected {spec['expected_total_pose_count']}, "
            f"actual {len(records)}"
        )
    records = select_contact_pattern_representatives(
        records, int(spec["representative_selection_rules"]["max_per_conformer"])
    )
    frequency_rows = build_contact_frequency_rows(
        records,
        spec["target_residues"],
        conformer_names,
        recurrent_threshold=float(
            spec["summary_rules"]["recurrent_contact_frequency_min"]
        ),
    )
    conformer_rows = build_conformer_summary(
        records,
        frequency_rows,
        conformer_names,
        float(spec["summary_rules"]["recurrent_contact_frequency_min"]),
    )
    representatives = sorted(
        (record for record in records if record["selection_status"] == "representative"),
        key=lambda record: (
            conformer_names.index(record["conformer"]),
            record["representative_rank_within_conformer"],
        ),
    )
    warnings = build_warnings(records, conformer_names, representatives)
    peptide_frequency_rows = [row for row in frequency_rows if row["scope"] == "peptide"]
    summary = build_summary_json(
        spec,
        records,
        peptide_frequency_rows,
        conformer_rows,
        representatives,
        input_details,
        warnings,
    )

    staging = safe_output.with_name(
        f".{safe_output.name}.incomplete_"
        + hashlib.sha256(str(safe_output).encode()).hexdigest()[:12]
    )
    if staging.exists():
        raise FileExistsError(f"Staging directory already exists: {staging}")
    staging.parent.mkdir(parents=True, exist_ok=True)
    staging.mkdir(exist_ok=False)
    write_tsv(staging / "pose_table.tsv", POSE_COLUMNS, records)
    write_tsv(staging / "pose_target_contacts.tsv", CONTACT_COLUMNS, contact_rows)
    write_tsv(
        staging / "target_contact_frequency.tsv", FREQUENCY_COLUMNS, frequency_rows
    )
    write_tsv(staging / "conformer_summary.tsv", CONFORMER_COLUMNS, conformer_rows)
    write_tsv(staging / "representative_poses.tsv", POSE_COLUMNS, representatives)
    (staging / "peptide_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (staging / "peptide_summary.md").write_text(
        build_summary_markdown(summary), encoding="utf-8"
    )
    (staging / "run_manifest.txt").write_text(
        build_manifest(resolved_spec, spec, data_root, safe_output, input_details, summary),
        encoding="utf-8",
    )
    staging.replace(safe_output)
    return safe_output


def build_conformer_summary(
    records: list[dict[str, Any]],
    frequency_rows: list[dict[str, Any]],
    conformer_names: list[str],
    recurrent_threshold: float,
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for conformer in conformer_names:
        rows = [record for record in records if record["conformer"] == conformer]
        frequencies = [
            row
            for row in frequency_rows
            if row["scope"] == "conformer" and row["conformer"] == conformer
        ]
        status_counts = Counter(record["selection_status"] for record in rows)
        pattern_counts = Counter(record["direct_contact_pattern"] for record in rows)
        dominant_patterns = [
            f"{pattern}:{count}"
            for pattern, count in sorted(
                pattern_counts.items(), key=lambda item: (-item[1], item[0])
            )[:3]
        ]
        output.append(
            {
                "conformer": conformer,
                "total_poses": len(rows),
                "receptor_clash_free_poses": sum(record["receptor_clash_free"] for record in rows),
                "chain_B_clash_free_poses": sum(record["chain_B_clash_free"] for record in rows),
                "interface_eligible_poses": sum(record["interface_eligible"] for record in rows),
                "representative_eligible_poses": sum(record["representative_eligible"] for record in rows),
                "kept_poses": status_counts["kept"],
                "filtered_clash_poses": status_counts["filtered_clash"],
                "filtered_low_interface_poses": status_counts["filtered_low_interface"],
                "representative_poses": status_counts["representative"],
                "best_vina_score": min(float(record["vina_score"]) for record in rows),
                "median_vina_score": median(float(record["vina_score"]) for record in rows),
                "worst_vina_score": max(float(record["vina_score"]) for record in rows),
                "median_legacy_5A_receptor_contact_fraction": median(
                    float(record["legacy_5A_receptor_contact_fraction"])
                    for record in rows
                ),
                "mean_direct_target_contact_count": mean(len(record["direct_target_contacts"]) for record in rows),
                "mean_near_target_contact_count": mean(len(record["near_target_contacts"]) for record in rows),
                "recurrent_direct_target_residues": [
                    row["target_residue_number"]
                    for row in frequencies
                    if row["direct_frequency"] >= recurrent_threshold
                ],
                "recurrent_near_target_residues": [
                    row["target_residue_number"]
                    for row in frequencies
                    if row["near_frequency"] >= recurrent_threshold
                ],
                "dominant_direct_contact_patterns": dominant_patterns,
                "representative_pose_ids": [
                    record["pose_id"]
                    for record in rows
                    if record["selection_status"] == "representative"
                ],
                "legacy_warning_mismatch_count": sum(
                    record["legacy_consistency_status"] == "WARNING_MISMATCH"
                    for record in rows
                ),
            }
        )
    return output


def build_warnings(
    records: list[dict[str, Any]],
    conformer_names: list[str],
    representatives: list[dict[str, Any]],
) -> list[str]:
    warnings: list[str] = []
    mismatches = [
        record["pose_id"]
        for record in records
        if record["legacy_consistency_status"] == "WARNING_MISMATCH"
    ]
    if mismatches:
        warnings.append("WARNING_MISMATCH poses: " + ", ".join(mismatches))
    for conformer in conformer_names:
        if not any(record["conformer"] == conformer for record in representatives):
            warnings.append(f"NO_REPRESENTATIVE_FOR_CONFORMER: {conformer}")
    return warnings


def build_summary_json(
    spec: dict[str, Any],
    records: list[dict[str, Any]],
    peptide_frequency_rows: list[dict[str, Any]],
    conformer_rows: list[dict[str, Any]],
    representatives: list[dict[str, Any]],
    input_details: dict[str, Any],
    warnings: list[str],
) -> dict[str, Any]:
    threshold = float(spec["summary_rules"]["recurrent_contact_frequency_min"])
    recurrent_direct = [
        row["target_residue_number"]
        for row in peptide_frequency_rows
        if row["conformer_balanced_direct_frequency"] >= threshold
    ]
    recurrent_near = [
        row["target_residue_number"]
        for row in peptide_frequency_rows
        if row["conformer_balanced_near_frequency"] >= threshold
    ]
    recurrent_by_conformer = {
        row["conformer"]: set(row["recurrent_direct_target_residues"])
        for row in conformer_rows
    }
    shared = shared_recurrent_direct_core(
        peptide_frequency_rows, len(conformer_rows)
    )
    pairwise: dict[str, float] = {}
    names = list(recurrent_by_conformer)
    for index, left in enumerate(names):
        for right in names[index + 1 :]:
            union = recurrent_by_conformer[left] | recurrent_by_conformer[right]
            pairwise[f"{left}__{right}"] = (
                len(recurrent_by_conformer[left] & recurrent_by_conformer[right])
                / len(union)
                if union
                else 1.0
            )
    return {
        "schema_version": "0.3",
        "analysis_id": spec["analysis_id"],
        "identity": {field: spec[field] for field in ("candidate", "sequence", "group")},
        "pose_counts": {
            "total": len(records),
            "by_conformer": {
                row["conformer"]: row["total_poses"] for row in conformer_rows
            },
        },
        "contact_definition": {
            "direct": "ligand/receptor heavy-atom minimum distance <= 4.0 A",
            "near": "4.0 A < minimum heavy-atom distance <= 5.0 A",
            "receptor": "monomer PDB, explicit spec chain",
            "ligand": "selected Vina MODEL PDBQT coordinates",
            "rounding_before_classification": False,
        },
        "contact_fraction_audit": {
            "metric_name": "legacy_5A_receptor_contact_fraction",
            "definition": "fraction of peptide residues with any heavy atom <= 5.0 A from any chain-A receptor heavy atom",
            "comparison": "distance <= 5.0 A using unrounded coordinates",
            "scope": "whole receptor interface, not target residues",
            "source": "legacy_pose_metrics",
            "generating_logic": "existing pilot analyze_results.py",
            "relationship_to_legacy_contacts": "same heavy-atom <=5.0 A receptor loop as legacy target contacts",
            "v03_role": "secondary_ranking_only",
            "v03_direct_hard_filter": False,
            "legacy_basic_geometry_pass_contact_fraction_threshold": 0.5,
        },
        "legacy_basic_geometry_pass": {
            "source": "legacy_pose_metrics",
            "role": "provenance_only_not_v03_eligibility",
            "definition": (
                "legacy_5A_receptor_contact_fraction >= 0.5 AND "
                "legacy_5A_target_contact_count >= 2 AND receptor_clash_pairs == 0 "
                "AND chain_B_clash_pairs == 0"
            ),
        },
        "selection_rules": spec["representative_selection_rules"],
        "selection_status_counts": dict(Counter(record["selection_status"] for record in records)),
        "recurrent_contacts": {
            "frequency_threshold": threshold,
            "direct_by_conformer_balanced_frequency": recurrent_direct,
            "near_by_conformer_balanced_frequency": recurrent_near,
        },
        "conformer_consistency": {
            "recurrent_direct_shared_core": shared,
            "recurrent_direct_by_conformer": {
                name: sorted(values) for name, values in recurrent_by_conformer.items()
            },
            "pairwise_recurrent_direct_jaccard": pairwise,
            "note": "contact-pattern comparison only; no RMSD or structural clustering",
        },
        "target_contact_frequency": peptide_frequency_rows,
        "conformer_summary": conformer_rows,
        "representative_poses": [
            {column: record.get(column) for column in POSE_COLUMNS}
            for record in representatives
        ],
        "input_files": input_details,
        "warnings": warnings,
        "provenance": {
            "contacts": "V0.2 heavy-atom contact primitive",
            "imported_metrics": "existing pose_metrics.tsv",
            "receptor_dimer": "validated and hashed only; chain-B metrics not recalculated",
        },
    }


def build_summary_markdown(summary: dict[str, Any]) -> str:
    counts = summary["pose_counts"]
    recurrent = summary["recurrent_contacts"]
    lines = [
        f"# {summary['analysis_id']} peptide ensemble analysis",
        "",
        f"- Total poses: `{counts['total']}`",
        "- Poses by conformer: `"
        + ", ".join(f"{key}={value}" for key, value in counts["by_conformer"].items())
        + "`",
        "- Recurrent direct target residues: `"
        + ";".join(map(str, recurrent["direct_by_conformer_balanced_frequency"]))
        + "`",
        "- Recurrent near target residues: `"
        + ";".join(map(str, recurrent["near_by_conformer_balanced_frequency"]))
        + "`",
        "",
        "## Legacy 5 A receptor contact-fraction provenance",
        "",
        summary["contact_fraction_audit"]["definition"] + ".",
        "`legacy_5A_receptor_contact_fraction` is retained as a secondary ranking "
        "value and is not a V0.3 hard filter.",
        "`legacy_basic_geometry_pass` is retained for provenance only and does not "
        "control V0.3 eligibility.",
        "",
        "## Contact-pattern representatives",
        "",
        "V0.3 selects representatives by eligible-pose contact-pattern prevalence; "
        "it does not perform RMSD clustering or claim structural diversity.",
        "",
    ]
    for representative in summary["representative_poses"]:
        lines.append(
            f"- `{representative['pose_id']}`: direct `"
            + ";".join(map(str, representative["direct_target_contacts"]))
            + f"`, score `{representative['vina_score']}`"
        )
    if summary["warnings"]:
        lines.extend(["", "## Warnings", ""])
        lines.extend(f"- {warning}" for warning in summary["warnings"])
    return "\n".join(lines) + "\n"


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: format_tsv_value(row.get(column)) for column in columns})


def format_tsv_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, list):
        return ";".join(map(str, value))
    if isinstance(value, float):
        return f"{value:.9f}"
    return value


def build_manifest(
    spec_path: Path,
    spec: dict[str, Any],
    data_root: Path,
    output_dir: Path,
    input_details: dict[str, Any],
    summary: dict[str, Any],
) -> str:
    lines = [
        f"run_datetime={datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"analysis_id={spec['analysis_id']}",
        f"data_root={data_root.resolve()}",
        f"spec={spec_path}",
        f"spec_sha256={sha256(spec_path)}",
        f"script_sha256={sha256(SCRIPT_PATH)}",
        f"ensemble_core_sha256={sha256(ENSEMBLE_CORE_PATH)}",
        f"pose_core_sha256={sha256(POSE_CORE_PATH)}",
        f"python_version={platform.python_version()}",
        f"python_executable={Path(sys.executable).resolve()}",
        f"pose_count.total={summary['pose_counts']['total']}",
        "contact_definition=heavy atoms; direct d<=4.0 A; near 4.0<d<=5.0 A; no pre-classification rounding",
        "legacy_5A_receptor_contact_fraction.source=legacy_pose_metrics",
        "legacy_5A_receptor_contact_fraction.role=secondary_ranking_only_not_hard_filter",
        "legacy_basic_geometry_pass.source=legacy_pose_metrics",
        "legacy_basic_geometry_pass.role=provenance_only_not_v03_eligibility",
        "pattern_prevalence.denominator=representative-eligible poses within each conformer",
        "representative_diversity=contact-pattern diversity only; no RMSD clustering",
        "chain_B_metrics=imported from pose_metrics.tsv; receptor dimer not recalculated",
    ]
    for field in ("receptor_monomer_path", "receptor_dimer_path"):
        detail = input_details[field]
        for key in ("spec_path", "resolved_path", "sha256", "expected_sha256", "hash_verification"):
            lines.append(f"input.{field}.{key}={detail.get(key) or 'NOT_PROVIDED'}")
    for conformer, conformer_details in input_details["conformers"].items():
        for field in ("docking_output_path", "pose_metrics_path"):
            detail = conformer_details[field]
            for key in ("spec_path", "resolved_path", "sha256", "expected_sha256", "hash_verification"):
                lines.append(
                    f"input.conformers.{conformer}.{field}.{key}="
                    f"{detail.get(key) or 'NOT_PROVIDED'}"
                )
        lines.append(
            f"input.conformers.{conformer}.model_labels="
            + ";".join(map(str, conformer_details["model_labels"]))
        )
    for warning in summary["warnings"]:
        lines.append(f"warning={warning}")
    for name in OUTPUT_NAMES:
        lines.append(f"output={output_dir / name}")
    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Analyze a peptide docking ensemble.")
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        output = analyze_peptide_ensemble(args.spec, args.data_root, args.output_dir)
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print(f"Analysis completed: {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
