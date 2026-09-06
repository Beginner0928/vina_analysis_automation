"""Regression validation for the V0.1 A4_conf02 MODEL 1 analysis bundle."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
GOLDEN_PATH = REPOSITORY_ROOT / "templates" / "A4_conf02_m1_golden.json"

REQUIRED_OUTPUTS = (
    "A4_conf02_m1_measurements.tsv",
    "A4_conf02_m1_contacts.tsv",
    "A4_conf02_m1_metrics.json",
    "A4_conf02_m1_summary.md",
    "run_manifest.txt",
    "A4_conf02_m1_interface.pse",
    "A4_conf02_m1_interface.png",
)


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return data


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) < 2:
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
        final = "FAIL" if self.failed else "PASS"
        return "\n".join(self.lines + ["", "FINAL RESULT:", final, ""])


def validate() -> tuple[Report, Path | None]:
    report = Report()
    if not GOLDEN_PATH.is_file():
        report.fail_item(f"golden JSON is missing: {GOLDEN_PATH}")
        return report, None

    golden = load_json(GOLDEN_PATH)
    output_dir = Path(golden["output_directory"])
    if not output_dir.is_dir():
        report.fail_item(f"analysis output directory is missing: {output_dir}")
        return report, None

    missing = [name for name in REQUIRED_OUTPUTS if not (output_dir / name).is_file()]
    if missing:
        for name in missing:
            report.fail_item(f"required analysis output is missing: {name}")
        return report, output_dir

    try:
        metrics = load_json(output_dir / "A4_conf02_m1_metrics.json")
        measurement_rows = load_tsv(output_dir / "A4_conf02_m1_measurements.tsv")
        contact_rows = load_tsv(output_dir / "A4_conf02_m1_contacts.tsv")
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        report.fail_item(f"cannot parse analysis results: {exc}")
        return report, output_dir

    report.compare_exact(
        "state_count", golden["expected_state_count"], metrics.get("state_count")
    )

    expected_contacts = set(golden["expected_target_contacts"])
    actual_contacts = {
        row.get("target_residue", "")
        for row in contact_rows
        if is_true(row.get("within_4A", ""))
    }
    if actual_contacts == expected_contacts:
        report.pass_item("target contacts = " + ", ".join(sorted(actual_contacts)))
    else:
        report.fail_item(
            "target contacts: expected "
            + ", ".join(sorted(expected_contacts))
            + "; actual "
            + ", ".join(sorted(actual_contacts))
        )

    pose_metrics = metrics.get("pose_metrics", {})
    tolerances = golden["tolerances"]
    report.compare_numeric(
        "Vina score",
        float(golden["expected_metrics"]["score"]),
        pose_metrics.get("score"),
        float(tolerances["score"]),
    )
    report.compare_numeric(
        "contact_fraction",
        float(golden["expected_metrics"]["contact_fraction"]),
        pose_metrics.get("contact_fraction"),
        float(tolerances["contact_fraction"]),
    )

    rows_by_interaction = {
        row.get("interaction", ""): row for row in measurement_rows
    }
    report.compare_exact(
        "measurement row count", len(golden["measurements"]), len(measurement_rows)
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
                report.fail_item(
                    f"{name} angle: expected blank, actual {actual_angle!r}"
                )
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

    return report, output_dir


def main() -> int:
    try:
        report, output_dir = validate()
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        report = Report()
        report.fail_item(f"validation error: {exc}")
        output_dir = None

    text = report.text()
    print(text, end="")
    if output_dir is not None:
        report_path = output_dir / "validation_report.txt"
        report_path.write_text(text, encoding="utf-8", newline="\n")
    return 1 if report.failed else 0


if __name__ == "__main__":
    sys.exit(main())
