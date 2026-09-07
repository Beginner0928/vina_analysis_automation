from __future__ import annotations

import sys
import unittest
import uuid
import hashlib
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = REPOSITORY_ROOT / "scripts"
TEST_OUTPUT_ROOT = REPOSITORY_ROOT / "test_output"
sys.path.insert(0, str(SCRIPTS_DIR))


class ModelBlockParsingTests(unittest.TestCase):
    def test_model_blocks_are_keyed_by_vina_model_label(self) -> None:
        from pose_analysis_core import parse_vina_model_blocks

        run_dir = TEST_OUTPUT_ROOT / f"v02_unit_{uuid.uuid4().hex}"
        run_dir.mkdir(parents=True, exist_ok=False)
        pdbqt = run_dir / "nonsequential_models.pdbqt"
        pdbqt.write_text(
            "MODEL 7\nATOM      1  C   LIG L   1       0.000   0.000   0.000\nENDMDL\n"
            "MODEL 42\nATOM      1  C   LIG L   1       1.000   0.000   0.000\nENDMDL\n",
            encoding="utf-8",
        )

        blocks = parse_vina_model_blocks(pdbqt)

        self.assertEqual([7, 42], list(blocks))
        self.assertTrue(blocks[7].startswith("MODEL 7\n"))
        self.assertTrue(blocks[42].startswith("MODEL 42\n"))

    def test_missing_model_label_fails_with_available_labels(self) -> None:
        from pose_analysis_core import require_model_block

        blocks = {7: "MODEL 7\nENDMDL\n", 42: "MODEL 42\nENDMDL\n"}

        with self.assertRaisesRegex(
            ValueError, r"MODEL 999.*available MODEL labels: 7, 42"
        ):
            require_model_block(blocks, 999)

    def test_duplicate_model_labels_are_rejected(self) -> None:
        from pose_analysis_core import parse_vina_model_blocks

        run_dir = TEST_OUTPUT_ROOT / f"v02_unit_{uuid.uuid4().hex}"
        run_dir.mkdir(parents=True, exist_ok=False)
        pdbqt = run_dir / "duplicate_models.pdbqt"
        pdbqt.write_text(
            "MODEL 7\nATOM      1  C   LIG L   1       0.000   0.000   0.000\nENDMDL\n"
            "MODEL 7\nATOM      1  C   LIG L   1       1.000   0.000   0.000\nENDMDL\n",
            encoding="utf-8",
        )

        with self.assertRaisesRegex(ValueError, r"Duplicate MODEL label 7"):
            parse_vina_model_blocks(pdbqt)


class InputResolutionTests(unittest.TestCase):
    def test_relative_input_paths_are_resolved_against_data_root_and_hashed(self) -> None:
        from pose_analysis_core import resolve_inputs

        run_dir = TEST_OUTPUT_ROOT / f"v02_unit_{uuid.uuid4().hex}"
        data_root = run_dir / "competition"
        input_path = data_root / "source" / "pose.pdbqt"
        input_path.parent.mkdir(parents=True, exist_ok=False)
        input_path.write_bytes(b"MODEL 1\nENDMDL\n")
        expected_hash = hashlib.sha256(input_path.read_bytes()).hexdigest()
        spec = {
            "docking_output_path": "source/pose.pdbqt",
            "pose_metrics_path": "source/pose.pdbqt",
            "receptor_monomer_path": "source/pose.pdbqt",
            "receptor_dimer_path": "source/pose.pdbqt",
            "input_sha256": {"docking_output_path": expected_hash},
        }

        inputs, actual_hashes = resolve_inputs(spec, data_root)

        self.assertEqual(input_path.resolve(), inputs["docking_output_path"])
        self.assertEqual(expected_hash, actual_hashes["docking_output_path"])

    def test_provided_input_hash_mismatch_is_fatal(self) -> None:
        from pose_analysis_core import resolve_inputs

        run_dir = TEST_OUTPUT_ROOT / f"v02_unit_{uuid.uuid4().hex}"
        data_root = run_dir / "competition"
        input_path = data_root / "source" / "pose.pdbqt"
        input_path.parent.mkdir(parents=True, exist_ok=False)
        input_path.write_bytes(b"MODEL 1\nENDMDL\n")
        spec = {
            "docking_output_path": "source/pose.pdbqt",
            "pose_metrics_path": "source/pose.pdbqt",
            "receptor_monomer_path": "source/pose.pdbqt",
            "receptor_dimer_path": "source/pose.pdbqt",
            "input_sha256": {"docking_output_path": "0" * 64},
        }

        with self.assertRaisesRegex(
            ValueError, r"SHA-256 mismatch.*docking_output_path"
        ):
            resolve_inputs(spec, data_root)


class ContactConsistencyTests(unittest.TestCase):
    def test_target_contact_mismatch_is_a_warning_with_differences(self) -> None:
        from pose_analysis_core import compare_target_contacts

        result = compare_target_contacts([150, 154, 160], [150, 159, 160])

        self.assertEqual("WARNING_MISMATCH", result["status"])
        self.assertEqual([154], result["only_in_v02_recomputed"])
        self.assertEqual([159], result["only_in_pose_metrics"])

    def test_absent_pose_metric_contacts_are_reported_as_not_available(self) -> None:
        from pose_analysis_core import compare_target_contacts

        result = compare_target_contacts([150, 154], "NA")

        self.assertEqual("NOT_AVAILABLE", result["status"])


class ContactClassificationTests(unittest.TestCase):
    def test_direct_near_and_none_boundaries_use_unrounded_distance(self) -> None:
        from pose_analysis_core import classify_contact

        cases = (
            (3.999999, "direct"),
            (4.0, "direct"),
            (4.000001, "near"),
            (5.0, "near"),
            (5.000001, "none"),
        )
        for distance_A, expected in cases:
            with self.subTest(distance_A=distance_A):
                self.assertEqual(
                    expected,
                    classify_contact(
                        distance_A, direct_cutoff_A=4.0, near_cutoff_A=5.0
                    ),
                )

    def test_legacy_only_near_residue_is_expected_definition_difference(self) -> None:
        from pose_analysis_core import compare_legacy_contacts

        result = compare_legacy_contacts(
            direct_contacts=[154, 161],
            near_contacts=[163],
            legacy_5A_contacts=[154, 161, 163],
        )

        self.assertEqual("EXPECTED_DEFINITION_DIFFERENCE", result["status"])
        self.assertEqual([163], result["legacy_only_explained_by_near"])
        self.assertEqual([], result["unexplained_residues"])

    def test_unexplained_legacy_difference_remains_warning(self) -> None:
        from pose_analysis_core import compare_legacy_contacts

        result = compare_legacy_contacts(
            direct_contacts=[154], near_contacts=[], legacy_5A_contacts=[154, 999]
        )

        self.assertEqual("WARNING_MISMATCH", result["status"])
        self.assertEqual([999], result["unexplained_residues"])

    def test_hydrogen_is_excluded_from_minimum_contact_distance(self) -> None:
        from pose_analysis_core import CoordinateAtom, compute_target_residue_contacts

        ligand_atoms = [
            CoordinateAtom(1, "H", "HD", (1.0, 0.0, 0.0), None, "UNL", ""),
            CoordinateAtom(2, "O", "OA", (4.5, 0.0, 0.0), None, "UNL", ""),
        ]
        receptor_atoms = [
            CoordinateAtom(10, "CA", "C", (0.0, 0.0, 0.0), 163, "GLU", "A")
        ]

        rows = compute_target_residue_contacts(
            ligand_atoms,
            receptor_atoms,
            target_residues=[163],
            direct_cutoff_A=4.0,
            near_cutoff_A=5.0,
        )

        self.assertEqual("near", rows[0]["classification"])
        self.assertEqual(4.5, rows[0]["min_distance_A"])
        self.assertEqual(2, rows[0]["ligand_atom_id"])


class SpecValidationTests(unittest.TestCase):
    def _valid_spec(self) -> dict[str, object]:
        return {
            "schema_version": "0.2",
            "analysis_id": "A4_conf03_m1",
            "candidate": "A4",
            "sequence": "SPSNFITMYDW",
            "group": "A",
            "conformer": "conf03",
            "model": 1,
            "docking_output_path": "docking/pose.pdbqt",
            "pose_metrics_path": "docking/metrics.tsv",
            "receptor_monomer_path": "receptors/A.pdb",
            "receptor_dimer_path": "receptors/AB.pdb",
            "receptor_chain": "A",
            "target_residues": [150, 154],
            "direct_contact_cutoff_A": 4.0,
            "near_contact_cutoff_A": 5.0,
            "manual_interactions": [],
        }

    def test_version_controlled_spec_rejects_absolute_input_paths(self) -> None:
        from pose_analysis_core import validate_spec

        spec = self._valid_spec()
        spec["docking_output_path"] = "E:/machine-specific/pose.pdbqt"

        with self.assertRaisesRegex(ValueError, r"must be relative to --data-root"):
            validate_spec(spec)

    def test_model_must_be_positive_integer_label(self) -> None:
        from pose_analysis_core import validate_spec

        spec = self._valid_spec()
        spec["model"] = 0

        with self.assertRaisesRegex(ValueError, r"model.*positive integer"):
            validate_spec(spec)

    def test_malformed_optional_hash_is_rejected(self) -> None:
        from pose_analysis_core import validate_spec

        spec = self._valid_spec()
        spec["input_sha256"] = {"docking_output_path": "not-a-hash"}

        with self.assertRaisesRegex(ValueError, r"input_sha256.*64 hexadecimal"):
            validate_spec(spec)

    def test_formal_contact_cutoffs_are_fixed(self) -> None:
        from pose_analysis_core import validate_spec

        spec = self._valid_spec()
        spec["direct_contact_cutoff_A"] = 4.1

        with self.assertRaisesRegex(ValueError, r"contact cutoffs are fixed"):
            validate_spec(spec)

    def test_pymol_selection_tokens_reject_expression_injection(self) -> None:
        from pose_analysis_core import validate_spec

        spec = self._valid_spec()
        spec["receptor_chain"] = "A or all"

        with self.assertRaisesRegex(ValueError, r"receptor_chain.*safe PyMOL token"):
            validate_spec(spec)


class PoseMetricsTests(unittest.TestCase):
    def test_pose_metrics_row_is_selected_by_model_label(self) -> None:
        from pose_analysis_core import read_model_metrics

        run_dir = TEST_OUTPUT_ROOT / f"v02_unit_{uuid.uuid4().hex}"
        run_dir.mkdir(parents=True, exist_ok=False)
        metrics_path = run_dir / "metrics.tsv"
        metrics_path.write_text(
            "model\tscore\ttarget_contacts\tchain_b_clash_pairs\n"
            "7\t-5.1\t150;154\t0\n"
            "42\t-6.2\t159;160\t1\n",
            encoding="utf-8",
        )

        header, metrics = read_model_metrics(metrics_path, 42)

        self.assertEqual(["model", "score", "target_contacts", "chain_b_clash_pairs"], header)
        self.assertEqual(42, metrics["model"])
        self.assertEqual(-6.2, metrics["score"])
        self.assertEqual([159, 160], metrics["target_contacts"])
        self.assertEqual(1, metrics["chain_b_clash_pairs"])


if __name__ == "__main__":
    unittest.main()
