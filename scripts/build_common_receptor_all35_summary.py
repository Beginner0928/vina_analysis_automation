#!/usr/bin/env python3
"""Merge the frozen Top12 and remaining23 common-receptor summaries."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


BLOCKED_CANDIDATES = {"A8", "B1", "D3", "D7", "D10"}
CRITICAL_SCORE_FIELDS = (
    "MeanBestScore3",
    "SDBestScore3",
    "GlobalBestScore",
    "Conf01BestScore",
    "Conf02BestScore",
    "Conf03BestScore",
)
MODEL_COUNT_FIELDS = (
    "Conf01ActualModelCount",
    "Conf02ActualModelCount",
    "Conf03ActualModelCount",
)


def read_summary(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        return list(reader.fieldnames), list(reader)


def finite_score(row: dict[str, str], field: str) -> float:
    raw = row.get(field, "")
    try:
        value = float(raw)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{row.get('Candidate', '<unknown>')}: invalid {field}={raw!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"{row.get('Candidate', '<unknown>')}: non-finite {field}={raw!r}")
    return value


def validate_rows(rows: list[dict[str, str]]) -> None:
    candidates = [row.get("Candidate", "").strip() for row in rows]
    if any(not candidate for candidate in candidates):
        raise ValueError("Every row must have a Candidate")
    duplicates = sorted({candidate for candidate in candidates if candidates.count(candidate) > 1})
    if duplicates:
        raise ValueError(f"Duplicate Candidate values: {duplicates}")
    blocked = sorted(BLOCKED_CANDIDATES.intersection(candidates))
    if blocked:
        raise ValueError(f"Blocked candidates must not enter the all35 summary: {blocked}")

    for row in rows:
        candidate = row["Candidate"]
        if not row.get("Sequence", "").strip():
            raise ValueError(f"{candidate}: Sequence is missing")
        if row.get("ConformerCount") != "3":
            raise ValueError(f"{candidate}: ConformerCount must equal 3")
        if row.get("EnsembleAnalysisStatus") != "PASS":
            raise ValueError(f"{candidate}: EnsembleAnalysisStatus must equal PASS")
        for field in CRITICAL_SCORE_FIELDS:
            finite_score(row, field)
        for field in MODEL_COUNT_FIELDS:
            raw = row.get(field, "")
            try:
                count = int(raw)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{candidate}: invalid {field}={raw!r}") from exc
            if not 1 <= count <= 20:
                raise ValueError(f"{candidate}: {field} must be within 1..20, got {count}")


def assign_ranks(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    ranked = sorted(rows, key=lambda row: (finite_score(row, "MeanBestScore3"), row["Candidate"]))
    for rank, row in enumerate(ranked, start=1):
        row["CommonReceptorOverallRank"] = str(rank)

    by_group: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in ranked:
        by_group[row["OriginalSelectionGroup"]].append(row)
    for group_rows in by_group.values():
        group_rows.sort(key=lambda row: (finite_score(row, "MeanBestScore3"), row["Candidate"]))
        for rank, row in enumerate(group_rows, start=1):
            row["WithinOriginalSelectionGroupRank"] = str(rank)
    return ranked


def build_summary(top12: Path, remaining23: Path, output: Path) -> None:
    top_fields, top_rows = read_summary(top12)
    remaining_fields, remaining_rows = read_summary(remaining23)
    if top_fields != remaining_fields:
        raise ValueError("Top12 and remaining23 CSV headers differ")
    if len(top_rows) != 12 or len(remaining_rows) != 23:
        raise ValueError(
            f"Expected 12 Top12 rows and 23 remaining rows, got {len(top_rows)} and {len(remaining_rows)}"
        )

    rows = top_rows + remaining_rows
    if len(rows) != 35:
        raise ValueError(f"Expected 35 combined rows, got {len(rows)}")
    validate_rows(rows)
    ranked = assign_ranks(rows)

    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {output}")
    with output.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=top_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(ranked)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top12", required=True, type=Path)
    parser.add_argument("--remaining23", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    build_summary(args.top12, args.remaining23, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
