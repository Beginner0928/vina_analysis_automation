"""V0.4 Vina-score parsing, cross-checking, and peptide aggregation."""

from __future__ import annotations

import math
import re
from statistics import mean, median
from typing import Any


VINA_RESULT_PATTERN = re.compile(
    r"^REMARK VINA RESULT:\s+([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\b"
)


def parse_vina_scores(model_blocks: dict[int, str]) -> dict[int, float]:
    """Extract exactly one primary Vina score from each MODEL block."""
    scores: dict[int, float] = {}
    for model, block in model_blocks.items():
        matches = [
            match
            for line in block.splitlines()
            if (match := VINA_RESULT_PATTERN.match(line)) is not None
        ]
        if len(matches) != 1:
            raise ValueError(
                f"MODEL {model} must contain exactly one REMARK VINA RESULT; "
                f"found {len(matches)}"
            )
        score = float(matches[0].group(1))
        if not math.isfinite(score):
            raise ValueError(f"MODEL {model} Vina score must be finite")
        scores[int(model)] = score
    return scores


def crosscheck_vina_scores(
    primary_scores: dict[int, float],
    pose_metrics: dict[int, dict[str, Any]],
    tolerance_kcal_mol: float,
) -> dict[int, dict[str, Any]]:
    """Compare output-remark scores with pose-metrics scores by MODEL label."""
    if tolerance_kcal_mol < 0:
        raise ValueError("Vina score tolerance must be non-negative")
    if set(primary_scores) != set(pose_metrics):
        raise ValueError("Vina score sources contain different MODEL labels")

    checks: dict[int, dict[str, Any]] = {}
    for model in sorted(primary_scores):
        primary = float(primary_scores[model])
        secondary = float(pose_metrics[model]["score"])
        difference = primary - secondary
        if abs(difference) > tolerance_kcal_mol:
            raise ValueError(
                "WARNING_SCORE_MISMATCH: "
                f"MODEL {model} docking-output score {primary} differs from "
                f"pose-metrics score {secondary} by {difference:+.6f} kcal/mol "
                f"(tolerance {tolerance_kcal_mol:.6f})"
            )
        checks[model] = {
            "model": model,
            "docking_output_pdbqt_remark_score": primary,
            "pose_metrics_score": secondary,
            "difference_kcal_mol": difference,
            "tolerance_kcal_mol": tolerance_kcal_mol,
            "status": "MATCH",
        }
    return checks


def build_vina_score_rows(
    records: list[dict[str, Any]], conformer_names: list[str]
) -> list[dict[str, Any]]:
    """Build one peptide row and arbitrary-count per-conformer score rows."""
    if not records:
        raise ValueError("Cannot summarize an empty pose ensemble")
    if not conformer_names or len(conformer_names) != len(set(conformer_names)):
        raise ValueError("conformer_names must be a non-empty unique list")

    by_conformer = {
        name: [record for record in records if record["conformer"] == name]
        for name in conformer_names
    }
    missing = [name for name, rows in by_conformer.items() if not rows]
    if missing:
        raise ValueError(f"No pose records for conformers: {missing}")
    unknown = sorted(
        set(str(record["conformer"]) for record in records) - set(conformer_names)
    )
    if unknown:
        raise ValueError(f"Pose records contain unknown conformers: {unknown}")

    conformer_rows = [
        _score_row("conformer", name, by_conformer[name]) for name in conformer_names
    ]
    overall = _score_row("peptide", "ALL", records)
    eligible_conformer_rows = [
        row for row in conformer_rows if row["eligible_pose_count"] > 0
    ]
    if eligible_conformer_rows:
        overall["mean_conformer_median_eligible_vina_score"] = mean(
            float(row["median_eligible_vina_score"])
            for row in eligible_conformer_rows
        )
        conformer_bests = [
            float(row["best_eligible_vina_score"])
            for row in eligible_conformer_rows
        ]
        overall["best_eligible_vina_score_range_across_conformers"] = (
            max(conformer_bests) - min(conformer_bests)
        )
    return [overall, *conformer_rows]


def _score_row(
    scope: str, conformer: str, records: list[dict[str, Any]]
) -> dict[str, Any]:
    raw = sorted(
        ((float(record["vina_score"]), str(record["pose_id"])) for record in records),
        key=lambda item: (item[0], item[1]),
    )
    if any(not math.isfinite(score) for score, _ in raw):
        raise ValueError("Vina scores must be finite")
    eligible = sorted(
        (
            (float(record["vina_score"]), str(record["pose_id"]))
            for record in records
            if bool(record["representative_eligible"])
        ),
        key=lambda item: (item[0], item[1]),
    )
    raw_scores = [score for score, _ in raw]
    eligible_scores = [score for score, _ in eligible]
    row = {
        "scope": scope,
        "conformer": conformer,
        "raw_pose_count": len(raw),
        "eligible_pose_count": len(eligible),
        "eligible_pose_fraction": len(eligible) / len(raw),
        "best_raw_vina_score": raw[0][0],
        "best_raw_pose_id": raw[0][1],
        "median_raw_vina_score": median(raw_scores),
        "best_eligible_vina_score": eligible[0][0] if eligible else None,
        "best_eligible_pose_id": eligible[0][1] if eligible else None,
        "median_eligible_vina_score": median(eligible_scores) if eligible else None,
        "mean_eligible_vina_score": mean(eligible_scores) if eligible else None,
        "q25_eligible_vina_score": _type7_quantile(eligible_scores, 0.25),
        "q75_eligible_vina_score": _type7_quantile(eligible_scores, 0.75),
        "mean_conformer_median_eligible_vina_score": None,
        "best_eligible_vina_score_range_across_conformers": None,
    }
    return row


def _type7_quantile(values: list[float], probability: float) -> float | None:
    """Return Hyndman-Fan type-7 quantile using explicit linear interpolation."""
    if not values:
        return None
    ordered = sorted(float(value) for value in values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] + fraction * (ordered[upper] - ordered[lower])
