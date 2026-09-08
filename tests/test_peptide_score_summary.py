from __future__ import annotations

import sys
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "scripts"))


class VinaRemarkScoreTests(unittest.TestCase):
    def test_scores_are_parsed_from_each_model_remark(self) -> None:
        from peptide_score_summary import parse_vina_scores

        blocks = {
            7: "MODEL 7\nREMARK VINA RESULT:    -7.125      0.000      0.000\nENDMDL\n",
            12: "MODEL 12\nREMARK VINA RESULT: -6.500 1.2 2.3\nENDMDL\n",
        }

        self.assertEqual({7: -7.125, 12: -6.5}, parse_vina_scores(blocks))

    def test_missing_or_duplicate_score_remark_is_rejected(self) -> None:
        from peptide_score_summary import parse_vina_scores

        with self.assertRaisesRegex(ValueError, r"MODEL 1.*exactly one"):
            parse_vina_scores({1: "MODEL 1\nENDMDL\n"})
        with self.assertRaisesRegex(ValueError, r"MODEL 2.*exactly one"):
            parse_vina_scores(
                {
                    2: (
                        "MODEL 2\n"
                        "REMARK VINA RESULT: -7.0 0 0\n"
                        "REMARK VINA RESULT: -6.0 0 0\n"
                        "ENDMDL\n"
                    )
                }
            )

    def test_pose_metric_scores_are_cross_checked_by_model_label(self) -> None:
        from peptide_score_summary import crosscheck_vina_scores

        metrics = {1: {"score": -7.0}, 3: {"score": -6.0}}
        checks = crosscheck_vina_scores({1: -7.0005, 3: -6.0}, metrics, 0.001)

        self.assertEqual("MATCH", checks[1]["status"])
        self.assertAlmostEqual(-0.0005, checks[1]["difference_kcal_mol"])
        self.assertEqual("MATCH", checks[3]["status"])

    def test_score_mismatch_is_fatal_and_explicit(self) -> None:
        from peptide_score_summary import crosscheck_vina_scores

        with self.assertRaisesRegex(
            ValueError, r"WARNING_SCORE_MISMATCH.*MODEL 2.*-7\.0.*-6\.9"
        ):
            crosscheck_vina_scores({2: -7.0}, {2: {"score": -6.9}}, 0.001)


class VinaScoreAggregationTests(unittest.TestCase):
    def test_raw_and_eligible_scores_are_kept_separate(self) -> None:
        from peptide_score_summary import build_vina_score_rows

        records = [
            self.record("c1_m1", "c1", -9.0, False),
            self.record("c1_m2", "c1", -8.0, True),
            self.record("c1_m3", "c1", -6.0, True),
            self.record("c2_m1", "c2", -7.0, True),
            self.record("c2_m2", "c2", -5.0, True),
        ]

        rows = build_vina_score_rows(records, ["c1", "c2"])
        overall = rows[0]

        self.assertEqual("peptide", overall["scope"])
        self.assertEqual(-9.0, overall["best_raw_vina_score"])
        self.assertEqual("c1_m1", overall["best_raw_pose_id"])
        self.assertEqual(-8.0, overall["best_eligible_vina_score"])
        self.assertEqual("c1_m2", overall["best_eligible_pose_id"])
        self.assertEqual(-7.0, overall["median_raw_vina_score"])
        self.assertEqual(-6.5, overall["median_eligible_vina_score"])
        self.assertEqual(-6.5, overall["mean_eligible_vina_score"])
        self.assertEqual(-7.25, overall["q25_eligible_vina_score"])
        self.assertEqual(-5.75, overall["q75_eligible_vina_score"])
        self.assertEqual(-6.5, overall["mean_conformer_median_eligible_vina_score"])
        self.assertEqual(1.0, overall["best_eligible_vina_score_range_across_conformers"])

        by_conformer = {row["conformer"]: row for row in rows[1:]}
        self.assertEqual(-8.0, by_conformer["c1"]["best_eligible_vina_score"])
        self.assertEqual(-7.0, by_conformer["c1"]["median_eligible_vina_score"])
        self.assertEqual(-7.0, by_conformer["c2"]["best_eligible_vina_score"])
        self.assertEqual(-6.0, by_conformer["c2"]["median_eligible_vina_score"])

    def test_filtered_poses_never_enter_eligible_statistics(self) -> None:
        from peptide_score_summary import build_vina_score_rows

        rows = build_vina_score_rows(
            [
                self.record("filtered", "only", -20.0, False),
                self.record("eligible", "only", -5.0, True),
            ],
            ["only"],
        )

        overall = rows[0]
        self.assertEqual(-20.0, overall["best_raw_vina_score"])
        self.assertEqual(-5.0, overall["best_eligible_vina_score"])
        self.assertEqual(1, overall["eligible_pose_count"])

    def test_no_eligible_poses_use_null_statistics(self) -> None:
        from peptide_score_summary import build_vina_score_rows

        rows = build_vina_score_rows(
            [self.record("filtered", "only", -9.0, False)], ["only"]
        )

        overall = rows[0]
        self.assertEqual(0, overall["eligible_pose_count"])
        self.assertIsNone(overall["best_eligible_vina_score"])
        self.assertIsNone(overall["best_eligible_pose_id"])
        self.assertIsNone(overall["mean_conformer_median_eligible_vina_score"])
        self.assertIsNone(
            overall["best_eligible_vina_score_range_across_conformers"]
        )

    def test_arbitrary_conformer_count_is_supported(self) -> None:
        from peptide_score_summary import build_vina_score_rows

        records = [
            self.record(f"{name}_m1", name, -float(index), True)
            for index, name in enumerate(("a", "b", "c", "d"), start=1)
        ]

        rows = build_vina_score_rows(records, ["a", "b", "c", "d"])

        self.assertEqual(5, len(rows))
        self.assertEqual(["ALL", "a", "b", "c", "d"], [row["conformer"] for row in rows])

    @staticmethod
    def record(
        pose_id: str, conformer: str, score: float, eligible: bool
    ) -> dict[str, object]:
        return {
            "pose_id": pose_id,
            "conformer": conformer,
            "vina_score": score,
            "representative_eligible": eligible,
        }


if __name__ == "__main__":
    unittest.main()
