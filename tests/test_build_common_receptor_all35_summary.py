from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_common_receptor_all35_summary.py"

FIELDS = [
    "Candidate",
    "OriginalSelectionGroup",
    "Sequence",
    "CommonReceptorOverallRank",
    "WithinOriginalSelectionGroupRank",
    "MeanBestScore3",
    "SDBestScore3",
    "GlobalBestScore",
    "Conf01ActualModelCount",
    "Conf02ActualModelCount",
    "Conf03ActualModelCount",
    "Conf01BestScore",
    "Conf02BestScore",
    "Conf03BestScore",
    "ConformerCount",
    "EnsembleAnalysisStatus",
]


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def make_row(candidate: str, score: float, group: str) -> dict[str, str]:
    return {
        "Candidate": candidate,
        "OriginalSelectionGroup": group,
        "Sequence": "ACDEFG",
        "CommonReceptorOverallRank": "999",
        "WithinOriginalSelectionGroupRank": "999",
        "MeanBestScore3": str(score),
        "SDBestScore3": "0.1",
        "GlobalBestScore": str(score - 0.2),
        "Conf01ActualModelCount": "19" if candidate == "A01" else "20",
        "Conf02ActualModelCount": "20",
        "Conf03ActualModelCount": "20",
        "Conf01BestScore": str(score),
        "Conf02BestScore": str(score),
        "Conf03BestScore": str(score),
        "ConformerCount": "3",
        "EnsembleAnalysisStatus": "PASS",
    }


class BuildCommonReceptorAll35SummaryTests(unittest.TestCase):
    def test_cli_merges_exactly_35_rows_and_recomputes_deterministic_ranks(self) -> None:
        self.assertTrue(SCRIPT.exists(), "the all35 summary builder must exist")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            top12 = root / "top12.csv"
            remaining23 = root / "remaining23.csv"
            output = root / "all35.csv"

            top_rows = [make_row(f"A{i:02d}", -7.0 + i / 10, "A") for i in range(1, 13)]
            remaining_rows = [
                make_row(f"B{i:02d}", -5.8 + i / 10, "B") for i in range(1, 24)
            ]
            # A12 and B01 tie; Candidate is the deterministic secondary key.
            remaining_rows[0]["MeanBestScore3"] = top_rows[-1]["MeanBestScore3"]
            write_rows(top12, top_rows)
            write_rows(remaining23, remaining_rows)

            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--top12",
                    str(top12),
                    "--remaining23",
                    str(remaining23),
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

            with output.open(encoding="utf-8-sig", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 35)
            self.assertEqual([int(row["CommonReceptorOverallRank"]) for row in rows], list(range(1, 36)))
            tied = [row["Candidate"] for row in rows if row["MeanBestScore3"] == "-5.8"]
            self.assertEqual(tied, ["A12", "B01"])
            self.assertEqual(len({row["Candidate"] for row in rows}), 35)


if __name__ == "__main__":
    unittest.main()
