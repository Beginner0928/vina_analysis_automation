from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from partial35_contract_v05 import PARTIAL35_CANDIDATE_IDS  # noqa: E402
from summarize_partial35_screen import build_group_summaries  # noqa: E402


class Partial35GroupSummaryTests(unittest.TestCase):
    def rows(self) -> list[dict]:
        rows = []
        for ordinal, candidate_id in enumerate(PARTIAL35_CANDIDATE_IDS):
            group = candidate_id[0]
            row = {
                    "candidate": candidate_id,
                    "group": group,
                    "best_eligible_vina_score": str(-5.0 - ordinal / 10),
                    "median_eligible_vina_score": str(-4.0 - ordinal / 20),
                    "eligible_pose_fraction": "0.75",
                    "shared_recurrent_direct_core": "154;159",
                    "representative_pose_ids": f"{candidate_id}_conf01_m1",
                    "warning_count": "0",
                    "analysis_status": "PASS",
                }
            row["_conformer_score_rows"] = [
                {
                    "conformer": conformer,
                    "best_eligible_vina_score": str(-5.0 - ordinal / 10 - index / 100),
                    "median_eligible_vina_score": str(-4.0 - ordinal / 20 - index / 100),
                }
                for index, conformer in enumerate(("conf01", "conf02", "conf03"))
            ]
            rows.append(row)
        return rows

    def test_group_summaries_preserve_partial_counts_and_no_cross_group_rank(self) -> None:
        result = build_group_summaries(self.rows())
        self.assertEqual(
            {"A": 9, "B": 9, "C": 10, "D": 7},
            {row["group"]: row["screened_candidate_count"] for row in result},
        )
        self.assertTrue(all(row["canonical_group_candidate_count"] == 10 for row in result))
        self.assertEqual(["A8"], result[0]["not_docked_candidate_ids"])
        self.assertEqual(["B1"], result[1]["not_docked_candidate_ids"])
        self.assertEqual([], result[2]["not_docked_candidate_ids"])
        self.assertEqual(["D3", "D7", "D10"], result[3]["not_docked_candidate_ids"])
        self.assertTrue(all(row["interpretation_scope"] == "WITHIN_GROUP_ONLY" for row in result))
        self.assertTrue(all("global_rank" not in row for row in result))
        self.assertEqual("AutoDock Vina predicted docking score", result[0]["score_semantics"])
        self.assertEqual(
            {"conf01", "conf02", "conf03"},
            set(result[0]["conformer_score_distribution"]),
        )
        self.assertEqual(9, len(result[0]["conformer_score_distribution"]["conf01"]))

    def test_missing_or_extra_candidate_fails_closed(self) -> None:
        rows = self.rows()
        with self.assertRaisesRegex(ValueError, "exact approved 35"):
            build_group_summaries(rows[:-1])
        changed = self.rows()
        changed.append(dict(changed[0], candidate="A8"))
        with self.assertRaisesRegex(ValueError, "exact approved 35"):
            build_group_summaries(changed)

    def test_candidate_technical_failure_is_reported_without_cross_group_ranking(self) -> None:
        rows = self.rows()[1:]
        statuses = [
            {
                "candidate_id": candidate_id,
                "state": "FAILED_DOCKING" if candidate_id == "A1" else "COMPLETE",
                "error_code": "DOCKING_OR_OUTPUT_AUDIT_FAILED" if candidate_id == "A1" else None,
            }
            for candidate_id in PARTIAL35_CANDIDATE_IDS
        ]
        result = build_group_summaries(rows, candidate_statuses=statuses)
        group_a = result[0]
        self.assertEqual(8, group_a["technically_complete_candidate_count"])
        self.assertEqual(1, group_a["technical_failure_candidate_count"])
        self.assertEqual("FAILED_DOCKING", group_a["candidate_statuses"][0]["state"])
        self.assertNotIn("global_rank", group_a)


if __name__ == "__main__":
    unittest.main()
