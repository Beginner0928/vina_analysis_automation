from __future__ import annotations

import json
import csv
import subprocess
import sys
import unittest
import uuid
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
COMPETITION_ROOT = REPOSITORY_ROOT.parent
TEST_OUTPUT_ROOT = REPOSITORY_ROOT / "test_output"
ANALYZER = REPOSITORY_ROOT / "scripts" / "analyze_pose.py"
GOLDEN_VALIDATOR = REPOSITORY_ROOT / "scripts" / "validate_pose_against_golden.py"
REQUIRED_OUTPUTS = (
    "metrics.json",
    "contacts.tsv",
    "measurements.tsv",
    "summary.md",
    "interface.pse",
    "interface.png",
    "run_manifest.txt",
)


def a4_conf03_spec(model: int) -> dict[str, object]:
    return {
        "schema_version": "0.2",
        "analysis_id": f"A4_conf03_m{model}",
        "candidate": "A4",
        "sequence": "SPSNFITMYDW",
        "group": "A",
        "conformer": "conf03",
        "model": model,
        "docking_output_path": "docking_test/vina_pilot_20260906/outputs/A4_conf03_out.pdbqt",
        "pose_metrics_path": "docking_test/vina_pilot_20260906/analysis/A4_conf03_pose_metrics.tsv",
        "receptor_monomer_path": "GPT_vina_branch_20260906/receptor_pdb/h26d3_af3_s2718_m2_TfR1_A.pdb",
        "receptor_dimer_path": "GPT_vina_branch_20260906/receptor_dimer/h26d3_af3_s2718_m2_TfR1_AB.pdb",
        "receptor_chain": "A",
        "target_residues": [150, 151, 154, 158, 159, 160, 161, 163, 385],
        "direct_contact_cutoff_A": 4.0,
        "near_contact_cutoff_A": 5.0,
        "manual_interactions": [],
    }


class AnalyzePoseCliNegativeTests(unittest.TestCase):
    def test_missing_model_999_fails_without_analysis_package(self) -> None:
        run_id = uuid.uuid4().hex
        spec_path = TEST_OUTPUT_ROOT / f"negative_model_999_{run_id}.json"
        output_dir = TEST_OUTPUT_ROOT / f"negative_model_999_{run_id}"
        TEST_OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
        spec_path.write_text(
            json.dumps(a4_conf03_spec(999), indent=2) + "\n", encoding="utf-8"
        )

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

        self.assertNotEqual(0, completed.returncode)
        self.assertIn("MODEL 999 is absent from Vina output", completed.stderr)
        self.assertFalse(output_dir.exists())


class AnalyzePoseCliSmokeTests(unittest.TestCase):
    def test_a4_conf03_model_1_generates_provenance_separated_package(self) -> None:
        run_id = uuid.uuid4().hex
        output_dir = TEST_OUTPUT_ROOT / f"A4_conf03_m1_v02_{run_id}"
        completed = subprocess.run(
            [
                sys.executable,
                str(ANALYZER),
                "--spec",
                str(REPOSITORY_ROOT / "specs" / "A4_conf03_m1_v02.json"),
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

        self.assertEqual(0, completed.returncode, completed.stderr)
        self.assertNotIn("[WARNING]", completed.stderr)
        for name in REQUIRED_OUTPUTS:
            path = output_dir / name
            self.assertTrue(path.is_file(), f"missing output: {name}")
            self.assertGreater(path.stat().st_size, 0, f"empty output: {name}")

        metrics = json.loads((output_dir / "metrics.json").read_text(encoding="utf-8"))
        recomputed = metrics["v02_recomputed"]
        imported = metrics["pose_metrics_imported"]
        self.assertEqual(1, recomputed["selected_model_label"])
        self.assertEqual(20, recomputed["state_count"])
        self.assertEqual(20, recomputed["parsed_model_block_count"])
        self.assertEqual(20, recomputed["pymol_state_count"])
        self.assertEqual(
            [154, 158, 159, 160, 161],
            recomputed["direct_target_contacts"],
        )
        self.assertEqual([163], recomputed["near_target_contacts"])
        self.assertEqual("v02_recomputed", recomputed["provenance"])
        self.assertEqual(
            "raw_PDB_and_PDBQT_coordinate_heavy_atom_recomputation",
            recomputed["field_provenance"]["direct_and_near_target_contacts"],
        )
        self.assertAlmostEqual(-6.202, imported["values"]["score"], places=6)
        self.assertEqual("existing_pose_metrics_tsv", imported["provenance"])
        self.assertEqual(
            "existing_pose_metrics_tsv", imported["field_provenance"]["chain_b_clash_pairs"]
        )
        contact_check = metrics["consistency_checks"]["legacy_5A_target_contacts"]
        self.assertEqual("EXPECTED_DEFINITION_DIFFERENCE", contact_check["status"])
        self.assertEqual([163], contact_check["legacy_only_explained_by_near"])
        self.assertEqual(
            [154, 158, 159, 160, 161, 163],
            metrics["pose_metrics_imported"]["legacy_5A_target_contacts"],
        )

        with (output_dir / "contacts.tsv").open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            contact_rows = list(csv.DictReader(handle, delimiter="\t"))
        glu163 = next(row for row in contact_rows if row["target_residue"] == "GLU163")
        self.assertEqual("near", glu163["classification"])
        self.assertAlmostEqual(4.905039, float(glu163["min_distance_A"]), places=6)
        self.assertEqual("60", glu163["ligand_atom_id"])
        self.assertEqual("O", glu163["ligand_atom_name"])
        self.assertEqual("OA", glu163["ligand_atom_type"])
        self.assertEqual("OE2", glu163["receptor_atom"])

        measurement_lines = (output_dir / "measurements.tsv").read_text(
            encoding="utf-8"
        ).splitlines()
        self.assertEqual(1, len(measurement_lines))
        summary = (output_dir / "summary.md").read_text(encoding="utf-8")
        self.assertIn("not recalculated from the AB dimer", summary)
        self.assertIn("Legacy-only residues explained by near contacts: `163`", summary)
        manifest = (output_dir / "run_manifest.txt").read_text(encoding="utf-8")
        self.assertIn(
            "result.legacy_contact_consistency=EXPECTED_DEFINITION_DIFFERENCE",
            manifest,
        )
        self.assertNotIn("WARNING_MISMATCH", manifest)

    def test_hashless_spec_still_records_actual_input_hashes(self) -> None:
        run_id = uuid.uuid4().hex
        source_spec = json.loads(
            (REPOSITORY_ROOT / "specs" / "A4_conf03_m1_v02.json").read_text(
                encoding="utf-8"
            )
        )
        source_spec.pop("input_sha256")
        spec_path = TEST_OUTPUT_ROOT / f"A4_conf03_m1_hashless_{run_id}.json"
        output_dir = TEST_OUTPUT_ROOT / f"A4_conf03_m1_hashless_{run_id}"
        spec_path.write_text(
            json.dumps(source_spec, indent=2) + "\n", encoding="utf-8"
        )

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

        self.assertEqual(0, completed.returncode, completed.stderr)
        manifest = (output_dir / "run_manifest.txt").read_text(encoding="utf-8")
        self.assertRegex(
            manifest, r"input\.docking_output_path\.sha256=[0-9a-f]{64}"
        )
        self.assertIn(
            "input.docking_output_path.expected_sha256=NOT_PROVIDED", manifest
        )
        self.assertIn(
            "input.docking_output_path.hash_verification=NOT_PROVIDED", manifest
        )


class AnalyzePoseCliGoldenRegressionTests(unittest.TestCase):
    def test_a4_conf02_model_1_passes_immutable_v01_golden(self) -> None:
        run_id = uuid.uuid4().hex
        output_dir = TEST_OUTPUT_ROOT / f"A4_conf02_m1_v02_{run_id}"
        analysis = subprocess.run(
            [
                sys.executable,
                str(ANALYZER),
                "--spec",
                str(REPOSITORY_ROOT / "specs" / "A4_conf02_m1_v02.json"),
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
        self.assertEqual(0, analysis.returncode, analysis.stderr)

        validation = subprocess.run(
            [
                sys.executable,
                str(GOLDEN_VALIDATOR),
                "--golden",
                str(REPOSITORY_ROOT / "templates" / "A4_conf02_m1_golden.json"),
                "--results-dir",
                str(output_dir),
            ],
            cwd=REPOSITORY_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(0, validation.returncode, validation.stdout + validation.stderr)
        self.assertIn("FINAL RESULT:\nPASS", validation.stdout.replace("\r\n", "\n"))


if __name__ == "__main__":
    unittest.main()
