from __future__ import annotations

import csv
import json
import subprocess
import sys
import unittest
import uuid
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
COMPETITION_ROOT = REPOSITORY_ROOT.parent
TEST_OUTPUT_ROOT = REPOSITORY_ROOT / "test_output"
ANALYZER = REPOSITORY_ROOT / "scripts" / "analyze_peptide_ensemble.py"
SPEC = REPOSITORY_ROOT / "specs" / "A4_ensemble_v03.json"
REQUIRED_OUTPUTS = {
    "pose_table.tsv",
    "pose_target_contacts.tsv",
    "target_contact_frequency.tsv",
    "conformer_summary.tsv",
    "representative_poses.tsv",
    "peptide_summary.md",
    "peptide_summary.json",
    "run_manifest.txt",
}


class A4EnsembleIntegrationTests(unittest.TestCase):
    def test_a4_enumerates_59_poses_and_expected_conformer_counts(self) -> None:
        output_dir, completed = self.run_analysis("complete")

        self.assertEqual(0, completed.returncode, completed.stderr)
        self.assertEqual(REQUIRED_OUTPUTS, {path.name for path in output_dir.iterdir()})
        with (output_dir / "pose_table.tsv").open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(59, len(rows))
        self.assertEqual(59, len({row["pose_id"] for row in rows}))
        counts = {
            conformer: sum(row["conformer"] == conformer for row in rows)
            for conformer in ("conf01", "conf02", "conf03")
        }
        self.assertEqual({"conf01": 20, "conf02": 19, "conf03": 20}, counts)
        by_pose = {row["pose_id"]: row for row in rows}
        for pose_id in (
            "A4_conf01_m5",
            "A4_conf01_m12",
            "A4_conf01_m13",
            "A4_conf03_m5",
        ):
            self.assertEqual("0", by_pose[pose_id]["legacy_basic_geometry_pass"])
            self.assertEqual("legacy_pose_metrics", by_pose[pose_id]["legacy_basic_geometry_pass_source"])
            self.assertEqual("filtered_low_interface", by_pose[pose_id]["selection_status"])
            self.assertNotIn("geometry_eligible", by_pose[pose_id])

        with (output_dir / "pose_target_contacts.tsv").open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            contact_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(59 * 9, len(contact_rows))
        glu163 = next(
            row
            for row in contact_rows
            if row["pose_id"] == "A4_conf03_m1"
            and row["target_residue"] == "GLU163"
        )
        self.assertEqual("near", glu163["classification"])
        self.assertAlmostEqual(4.905039, float(glu163["min_distance_A"]), places=6)

        summary = json.loads(
            (output_dir / "peptide_summary.json").read_text(encoding="utf-8")
        )
        self.assertEqual(59, summary["pose_counts"]["total"])
        self.assertEqual(
            {"conf01": 20, "conf02": 19, "conf03": 20},
            summary["pose_counts"]["by_conformer"],
        )
        self.assertEqual(
            "secondary_ranking_only",
            summary["contact_fraction_audit"]["v03_role"],
        )
        self.assertIs(
            False,
            summary["contact_fraction_audit"].get("v03_direct_hard_filter"),
        )
        self.assertEqual(
            0.5,
            summary["contact_fraction_audit"].get(
                "legacy_basic_geometry_pass_contact_fraction_threshold"
            ),
        )
        self.assertEqual(
            {"kept": 35, "filtered_low_interface": 15, "representative": 9},
            summary["selection_status_counts"],
        )
        frequency_by_residue = {
            row["target_residue_number"]: row
            for row in summary["target_contact_frequency"]
        }
        self.assertEqual(3, frequency_by_residue[150]["conformer_presence_count"])
        self.assertEqual(1, frequency_by_residue[150]["recurrent_conformer_support_count"])
        self.assertEqual(
            [154, 159, 160, 161],
            summary["conformer_consistency"]["recurrent_direct_shared_core"],
        )

    def test_missing_conformer_input_fails_without_final_package(self) -> None:
        spec = json.loads(SPEC.read_text(encoding="utf-8"))
        spec["conformers"][0]["docking_output_path"] = "missing/conf01.pdbqt"
        output_dir, completed = self.run_analysis("missing", spec)

        self.assertNotEqual(0, completed.returncode)
        self.assertIn("conf01", completed.stderr)
        self.assertIn("docking_output_path", completed.stderr)
        self.assertFalse(output_dir.exists())

    def test_wrong_expected_model_count_fails_explicitly(self) -> None:
        spec = json.loads(SPEC.read_text(encoding="utf-8"))
        spec["conformers"][0]["expected_model_count"] = 21
        spec["expected_total_pose_count"] = 60
        output_dir, completed = self.run_analysis("wrong_count", spec)

        self.assertNotEqual(0, completed.returncode)
        self.assertRegex(completed.stderr, r"conf01.*expected 21.*actual 20")
        self.assertFalse(output_dir.exists())

    def run_analysis(
        self, label: str, spec: dict[str, object] | None = None
    ) -> tuple[Path, subprocess.CompletedProcess[str]]:
        run_id = uuid.uuid4().hex
        output_dir = TEST_OUTPUT_ROOT / f"A4_v03_{label}_{run_id}" / "A4"
        spec_path = SPEC
        if spec is not None:
            spec_path = TEST_OUTPUT_ROOT / f"A4_v03_{label}_{run_id}.json"
            spec_path.parent.mkdir(parents=True, exist_ok=True)
            spec_path.write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8")
        completed = subprocess.run(
            [
                sys.executable,
                str(ANALYZER),
                "--spec",
                str(spec_path),
                "--data-root",
                str(COMPETITION_ROOT),
                "--output-dir",
                str(output_dir),
            ],
            cwd=REPOSITORY_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        return output_dir, completed


if __name__ == "__main__":
    unittest.main()
