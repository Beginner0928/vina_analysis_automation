"""Produce within-group summaries for a completed PARTIAL35 exploratory screen."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from collections import Counter
from pathlib import Path
from typing import Any

from batch_execution_v05 import file_sha256
from partial35_contract_v05 import EXCLUDED_CANDIDATES, PARTIAL35_CANDIDATE_IDS


PROJECT_ROOT = Path(__file__).resolve().parents[1]

GROUP_SUMMARY_FIELDS = (
    "group",
    "canonical_group_candidate_count",
    "screened_candidate_count",
    "technically_complete_candidate_count",
    "technical_failure_candidate_count",
    "candidate_statuses",
    "not_docked_candidate_ids",
    "not_docked_reason",
    "best_eligible_vina_score",
    "best_eligible_candidate_id",
    "median_of_candidate_median_eligible_vina_scores",
    "mean_eligible_pose_fraction",
    "conformer_score_distribution",
    "candidate_median_score_distribution",
    "recurrent_direct_contact_candidate_frequency",
    "representative_pose_ids",
    "warning_count",
    "analysis_status",
    "score_semantics",
    "interpretation_scope",
)


def _list_value(value: str) -> list[Any]:
    if not value:
        return []
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError:
        return value.split(";")
    return parsed if isinstance(parsed, list) else [parsed]


def build_group_summaries(
    candidate_rows: list[dict[str, Any]],
    *,
    candidate_statuses: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    actual_ids = [str(row.get("candidate")) for row in candidate_rows]
    if candidate_statuses is None:
        if actual_ids != list(PARTIAL35_CANDIDATE_IDS):
            raise ValueError("Candidate summary rows must be the exact approved 35 in canonical order")
        candidate_statuses = [
            {"candidate_id": candidate_id, "state": "COMPLETE", "error_code": None}
            for candidate_id in PARTIAL35_CANDIDATE_IDS
        ]
    status_ids = [str(row.get("candidate_id")) for row in candidate_statuses]
    if status_ids != list(PARTIAL35_CANDIDATE_IDS):
        raise ValueError("Candidate statuses must be the exact approved 35 in canonical order")
    complete_ids = [
        str(row["candidate_id"])
        for row in candidate_statuses
        if row.get("state") == "COMPLETE"
    ]
    if actual_ids != complete_ids:
        raise ValueError("Candidate summary rows must match COMPLETE candidate statuses in canonical order")
    excluded_by_group: dict[str, list[str]] = {group: [] for group in "ABCD"}
    for candidate_id, _reason in EXCLUDED_CANDIDATES:
        excluded_by_group[candidate_id[0]].append(candidate_id)

    results: list[dict[str, Any]] = []
    for group in "ABCD":
        rows = [row for row in candidate_rows if row.get("group") == group]
        statuses = [row for row in candidate_statuses if str(row["candidate_id"]).startswith(group)]
        failed_statuses = [row for row in statuses if row.get("state") != "COMPLETE"]
        score_pairs = sorted(
            (
                (float(row["best_eligible_vina_score"]), str(row["candidate"]))
                for row in rows
            ),
            key=lambda item: (item[0], item[1]),
        )
        medians = [float(row["median_eligible_vina_score"]) for row in rows]
        fractions = [float(row["eligible_pose_fraction"]) for row in rows]
        contact_support: Counter[int] = Counter()
        representatives: list[str] = []
        conformer_distribution: dict[str, list[dict[str, Any]]] = {
            conformer: [] for conformer in ("conf01", "conf02", "conf03")
        }
        for row in rows:
            contacts = _list_value(str(row.get("shared_recurrent_direct_core", "")))
            contact_support.update(int(value) for value in contacts)
            representatives.extend(
                str(value)
                for value in _list_value(str(row.get("representative_pose_ids", "")))
            )
            conformer_rows = row.get("_conformer_score_rows", [])
            if not isinstance(conformer_rows, list) or {
                str(value.get("conformer")) for value in conformer_rows
            } != set(conformer_distribution):
                raise ValueError(
                    f"Candidate {row['candidate']} lacks the exact three conformer score rows"
                )
            for value in conformer_rows:
                conformer_distribution[str(value["conformer"])].append(
                    {
                        "candidate_id": row["candidate"],
                        "best_eligible_vina_score": float(
                            value["best_eligible_vina_score"]
                        ),
                        "median_eligible_vina_score": float(
                            value["median_eligible_vina_score"]
                        ),
                    }
                )
        warnings = sum(int(row.get("warning_count", 0)) for row in rows)
        results.append(
            {
                "group": group,
                "canonical_group_candidate_count": 10,
                "screened_candidate_count": len(statuses),
                "technically_complete_candidate_count": len(rows),
                "technical_failure_candidate_count": len(failed_statuses),
                "candidate_statuses": statuses,
                "not_docked_candidate_ids": excluded_by_group[group],
                "not_docked_reason": (
                    "CONFORMER_GENERATION_QC_FAIL_NOT_DOCKED"
                    if excluded_by_group[group]
                    else None
                ),
                "best_eligible_vina_score": score_pairs[0][0] if score_pairs else None,
                "best_eligible_candidate_id": score_pairs[0][1] if score_pairs else None,
                "median_of_candidate_median_eligible_vina_scores": (
                    statistics.median(medians) if medians else None
                ),
                "mean_eligible_pose_fraction": (
                    statistics.mean(fractions) if fractions else None
                ),
                "conformer_score_distribution": conformer_distribution,
                "candidate_median_score_distribution": [
                    {
                        "candidate_id": row["candidate"],
                        "median_eligible_vina_score": float(row["median_eligible_vina_score"]),
                    }
                    for row in rows
                ],
                "recurrent_direct_contact_candidate_frequency": {
                    str(residue): count / len(rows)
                    for residue, count in sorted(contact_support.items())
                } if rows else {},
                "representative_pose_ids": representatives,
                "warning_count": warnings,
                "analysis_status": (
                    "PASS_WITH_TECHNICAL_FAILURES"
                    if failed_statuses
                    else ("PASS_WITH_WARNINGS" if warnings else "PASS")
                ),
                "score_semantics": "AutoDock Vina predicted docking score",
                "interpretation_scope": "WITHIN_GROUP_ONLY",
            }
        )
    return results


def _load_complete_candidate_row(batch_root: Path, candidate_id: str) -> dict[str, Any]:
    attempts = batch_root / "candidates" / candidate_id / "attempts"
    completion_paths = sorted(attempts.glob("*/completion.json")) if attempts.is_dir() else []
    valid: list[tuple[Path, dict[str, Any]]] = []
    for path in completion_paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("status") == "COMPLETE" and payload.get("integrity_status") == "PASS":
            valid.append((path.parent, payload))
    if len(valid) != 1:
        raise ValueError(f"Expected exactly one valid COMPLETE attempt for {candidate_id}; found {len(valid)}")
    attempt, completion = valid[0]
    artifact = completion.get("artifacts", {}).get("candidate_summary")
    if not isinstance(artifact, dict):
        raise ValueError(f"Completion marker lacks candidate_summary for {candidate_id}")
    summary = (attempt / artifact["path"]).resolve()
    if not summary.is_relative_to(attempt.resolve()) or file_sha256(summary) != artifact["sha256"]:
        raise ValueError(f"candidate_summary integrity mismatch for {candidate_id}")
    with summary.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1 or rows[0].get("candidate") != candidate_id:
        raise ValueError(f"candidate_summary identity mismatch for {candidate_id}")
    score_artifact = completion.get("artifacts", {}).get(
        "ensemble.vina_score_summary.tsv"
    )
    if not isinstance(score_artifact, dict):
        raise ValueError(f"Completion marker lacks vina_score_summary for {candidate_id}")
    score_path = (attempt / score_artifact["path"]).resolve()
    if (
        not score_path.is_relative_to(attempt.resolve())
        or file_sha256(score_path) != score_artifact["sha256"]
    ):
        raise ValueError(f"vina_score_summary integrity mismatch for {candidate_id}")
    with score_path.open("r", encoding="utf-8-sig", newline="") as handle:
        score_rows = list(csv.DictReader(handle, delimiter="\t"))
    conformer_rows = [row for row in score_rows if row.get("scope") == "conformer"]
    if [row.get("conformer") for row in conformer_rows] != [
        "conf01",
        "conf02",
        "conf03",
    ]:
        raise ValueError(f"vina_score_summary conformer inventory mismatch for {candidate_id}")
    rows[0]["_conformer_score_rows"] = conformer_rows
    return rows[0]


def summarize_partial35(batch_root: Path, output_root: Path) -> Path:
    if output_root.resolve().is_relative_to(PROJECT_ROOT):
        raise ValueError("PARTIAL35 scientific summary output must remain outside the Git repository")
    if output_root.exists():
        raise FileExistsError(f"Refusing to overwrite existing PARTIAL35 summary: {output_root}")
    batch_runs = sorted((batch_root / "batch_runs").glob("[0-9][0-9][0-9][0-9].json"))
    if not batch_runs:
        raise FileNotFoundError("PARTIAL35 batch has no batch run record")
    batch_run = json.loads(batch_runs[-1].read_text(encoding="utf-8"))
    candidate_statuses = batch_run.get("candidates")
    if not isinstance(candidate_statuses, list):
        raise ValueError("PARTIAL35 batch run lacks candidate status records")
    complete_ids = [
        str(row["candidate_id"])
        for row in candidate_statuses
        if row.get("state") == "COMPLETE"
    ]
    rows = [_load_complete_candidate_row(batch_root, candidate_id) for candidate_id in complete_ids]
    summaries = build_group_summaries(rows, candidate_statuses=candidate_statuses)
    output_root.mkdir(parents=True)
    json_path = output_root / "partial35_group_summary.json"
    with json_path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(
            {
                "schema_version": "0.5-partial35-summary",
                "package_classification": "PARTIAL_EXPLORATORY_SCREENING",
                "completeness_statement": "NOT_A_COMPLETE_TOP40_SCREEN",
                "screened_candidate_count": 35,
                "technically_complete_candidate_count": len(complete_ids),
                "technical_failure_candidate_count": 35 - len(complete_ids),
                "excluded_candidates": [
                    {
                        "candidate_id": candidate_id,
                        "status": "CONFORMER_GENERATION_QC_FAIL_NOT_DOCKED",
                        "reason": reason,
                    }
                    for candidate_id, reason in EXCLUDED_CANDIDATES
                ],
                "cross_group_global_ranking": "NOT_PERFORMED",
                "groups": summaries,
                "status": "PASS",
            },
            handle,
            indent=2,
            ensure_ascii=False,
        )
        handle.write("\n")
    with (output_root / "partial35_group_summary.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=GROUP_SUMMARY_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(
            {
                field: (
                    json.dumps(row[field], separators=(",", ":"))
                    if isinstance(row[field], (list, dict))
                    else row[field]
                )
                for field in GROUP_SUMMARY_FIELDS
            }
            for row in summaries
        )
    with (output_root / "partial35_summary.md").open(
        "x", encoding="utf-8", newline="\n"
    ) as handle:
        handle.write("# PARTIAL35 exploratory screening summary\n\n")
        handle.write("This is not a complete Top40 screen. Comparisons are within group only.\n\n")
        handle.write("The reported values are AutoDock Vina predicted docking scores, not experimental binding free energies.\n\n")
        for row in summaries:
            handle.write(
                f"## Group {row['group']}\n\n"
                f"- Screened: {row['screened_candidate_count']}/10\n"
                f"- Not docked: {', '.join(row['not_docked_candidate_ids']) or 'none'}\n"
                f"- Best eligible predicted score: {row['best_eligible_vina_score']} ({row['best_eligible_candidate_id']})\n"
                f"- Warning count: {row['warning_count']}\n\n"
            )
    return output_root


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        output = summarize_partial35(args.batch_root.resolve(), args.output_root.resolve())
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print(f"PARTIAL35 group summary PASS: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
