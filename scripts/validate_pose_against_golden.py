"""Validate a V0.2 generic pose-analysis package against the V0.1 Golden case."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any


REQUIRED_OUTPUTS = (
    "metrics.json",
    "contacts.tsv",
    "measurements.tsv",
    "summary.md",
    "interface.pse",
    "interface.png",
    "run_manifest.txt",
)


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Not a valid tab-separated table: {path}")
        return list(reader)


def parse_float(value: Any, label: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} is not numeric: {value!r}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{label} is not finite: {value!r}")
    return parsed


def is_true(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


class Report:
    def __init__(self) -> None:
        self.lines: list[str] = []
        self.failed = False

    def pass_item(self, message: str) -> None:
        self.lines.append(f"[PASS] {message}")

    def fail_item(self, message: str) -> None:
        self.failed = True
        self.lines.append(f"[FAIL] {message}")

    def compare_exact(self, label: str, expected: Any, actual: Any) -> None:
        if actual == expected:
            self.pass_item(f"{label} = {actual}")
        else:
            self.fail_item(f"{label}: expected {expected!r}, actual {actual!r}")

    def compare_numeric(
        self, label: str, expected: float, actual: Any, tolerance: float
    ) -> None:
        try:
            actual_value = parse_float(actual, label)
        except ValueError as exc:
            self.fail_item(str(exc))
            return
        delta = actual_value - expected
        status = "PASS" if abs(delta) <= tolerance else "FAIL"
        self.lines.extend(
            [
                f"[{status}] {label}",
                f"expected {expected:.6f}",
                f"actual {actual_value:.6f}",
                f"delta {delta:+.6f}",
            ]
        )
        if status == "FAIL":
            self.failed = True

    def text(self) -> str:
        result = "FAIL" if self.failed else "PASS"
        return "\n".join(self.lines + ["", "FINAL RESULT:", result, ""])


def validate(golden_path: Path, results_dir: Path) -> Report:
    report = Report()
    golden = load_json(golden_path)
    missing = [name for name in REQUIRED_OUTPUTS if not (results_dir / name).is_file()]
    if missing:
        for name in missing:
            report.fail_item(f"required V0.2 output is missing: {name}")
        return report

    metrics = load_json(results_dir / "metrics.json")
    contacts = load_tsv(results_dir / "contacts.tsv")
    measurements = load_tsv(results_dir / "measurements.tsv")
    recomputed = metrics.get("v02_recomputed", {})
    imported = metrics.get("pose_metrics_imported", {}).get("values", {})

    report.compare_exact(
        "state_count", golden["expected_state_count"], recomputed.get("state_count")
    )
    report.compare_exact(
        "selected MODEL label", golden["model"], recomputed.get("selected_model_label")
    )
    actual_contacts = {
        row.get("target_residue", "")
        for row in contacts
        if row.get("classification", "") == "direct"
    }
    expected_contacts = set(golden["expected_target_contacts"])
    if actual_contacts == expected_contacts:
        report.pass_item("target contacts = " + ", ".join(sorted(actual_contacts)))
    else:
        report.fail_item(
            "target contacts: expected "
            + ", ".join(sorted(expected_contacts))
            + "; actual "
            + ", ".join(sorted(actual_contacts))
        )

    tolerances = golden["tolerances"]
    expected_metrics = golden["expected_metrics"]
    report.compare_numeric(
        "Vina score",
        float(expected_metrics["score"]),
        imported.get("score"),
        float(tolerances["score"]),
    )
    report.compare_numeric(
        "contact_fraction",
        float(expected_metrics["contact_fraction"]),
        imported.get("contact_fraction"),
        float(tolerances["contact_fraction"]),
    )
    for field in ("target_contact_count", "receptor_clash_pairs", "chain_b_clash_pairs", "basic_geometry_pass"):
        report.compare_exact(field, expected_metrics[field], imported.get(field))
    for field in ("min_receptor_distance", "min_chain_b_distance"):
        report.compare_numeric(
            field,
            float(expected_metrics[field]),
            imported.get(field),
            1e-9,
        )

    rows_by_interaction = {row.get("interaction", ""): row for row in measurements}
    report.compare_exact(
        "measurement row count", len(golden["measurements"]), len(measurements)
    )
    for expected in golden["measurements"]:
        name = expected["interaction"]
        actual = rows_by_interaction.get(name)
        if actual is None:
            report.fail_item(f"measurement is missing: {name}")
            continue
        report.compare_numeric(
            f"{name} distance",
            float(expected["distance_A"]),
            actual.get("distance_A"),
            float(tolerances["distance_A"]),
        )
        expected_angle = expected.get("angle_deg")
        actual_angle = actual.get("angle_deg", "").strip()
        if expected_angle is None:
            if actual_angle == "":
                report.pass_item(f"{name} angle is blank")
            else:
                report.fail_item(f"{name} angle: expected blank, actual {actual_angle!r}")
        else:
            report.compare_numeric(
                f"{name} angle",
                float(expected_angle),
                actual_angle,
                float(tolerances["angle_deg"]),
            )
        report.compare_exact(
            f"{name} interpretation_source",
            golden["interpretation_source"],
            actual.get("interpretation_source"),
        )
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare a V0.2 pose package with the immutable V0.1 Golden case."
    )
    parser.add_argument("--golden", required=True, type=Path)
    parser.add_argument("--results-dir", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        report = validate(args.golden.resolve(), args.results_dir.resolve())
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        report = Report()
        report.fail_item(f"validation error: {exc}")
    print(report.text(), end="")
    return 1 if report.failed else 0


if __name__ == "__main__":
    sys.exit(main())
