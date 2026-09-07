from __future__ import annotations

import sys
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "scripts"))


def base_spec() -> dict[str, object]:
    return {
        "schema_version": "0.3",
        "analysis_id": "synthetic_ensemble",
        "candidate": "X1",
        "sequence": "AAAA",
        "group": "X",
        "receptor_monomer_path": "receptors/A.pdb",
        "receptor_dimer_path": "receptors/AB.pdb",
        "receptor_chain": "A",
        "target_residues": [10, 20],
        "expected_total_pose_count": 3,
        "conformers": [
            {
                "conformer_name": "c1",
                "docking_output_path": "outputs/c1.pdbqt",
                "pose_metrics_path": "metrics/c1.tsv",
                "expected_model_count": 2,
            },
            {
                "conformer_name": "c2",
                "docking_output_path": "outputs/c2.pdbqt",
                "pose_metrics_path": "metrics/c2.tsv",
                "expected_model_count": 1,
            },
        ],
        "contact_cutoffs": {"direct_A": 4.0, "near_A": 5.0},
        "representative_selection_rules": {
            "max_receptor_clash_pairs": 0,
            "max_chain_B_clash_pairs": 0,
            "min_direct_target_contacts": 2,
            "required_any_direct_residues": [],
            "max_per_conformer": 3,
            "diversity_key": "exact_direct_contact_pattern",
        },
        "summary_rules": {"recurrent_contact_frequency_min": 0.25},
    }


class EnsembleSpecTests(unittest.TestCase):
    def test_generic_spec_accepts_non_59_total_equal_to_conformer_sum(self) -> None:
        from peptide_ensemble_core import validate_ensemble_spec

        validate_ensemble_spec(base_spec())

    def test_expected_total_must_equal_conformer_sum(self) -> None:
        from peptide_ensemble_core import validate_ensemble_spec

        spec = base_spec()
        spec["expected_total_pose_count"] = 59

        with self.assertRaisesRegex(
            ValueError, r"expected_total_pose_count 59.*conformer.*3"
        ):
            validate_ensemble_spec(spec)

    def test_model_inventory_requires_exact_labels_and_count(self) -> None:
        from peptide_ensemble_core import validate_model_inventory

        with self.assertRaisesRegex(
            ValueError, r"conf01.*expected 3.*actual 2"
        ):
            validate_model_inventory("conf01", {1, 2}, {1, 2}, 3)

        with self.assertRaisesRegex(
            ValueError, r"only in PDBQT.*2.*only in metrics.*3"
        ):
            validate_model_inventory("conf01", {1, 2}, {1, 3}, 2)

    def test_machine_specific_absolute_paths_are_rejected(self) -> None:
        from peptide_ensemble_core import validate_ensemble_spec

        spec = base_spec()
        spec["receptor_monomer_path"] = "E:/local/receptor.pdb"

        with self.assertRaisesRegex(ValueError, r"relative to --data-root"):
            validate_ensemble_spec(spec)

    def test_malformed_optional_input_hash_is_rejected(self) -> None:
        from peptide_ensemble_core import validate_ensemble_spec

        spec = base_spec()
        spec["input_sha256"] = {"receptor_monomer_path": "not-a-hash"}

        with self.assertRaisesRegex(ValueError, r"input_sha256.*64 hexadecimal"):
            validate_ensemble_spec(spec)

    def test_contact_fraction_cannot_be_configured_as_hard_filter(self) -> None:
        from peptide_ensemble_core import validate_ensemble_spec

        spec = base_spec()
        rules = spec["representative_selection_rules"]
        assert isinstance(rules, dict)
        rules["min_contact_fraction"] = 0.6

        with self.assertRaisesRegex(ValueError, r"unsupported.*min_contact_fraction"):
            validate_ensemble_spec(spec)


class ContactFrequencyTests(unittest.TestCase):
    def test_balanced_frequency_equal_weights_conformers(self) -> None:
        from peptide_ensemble_core import build_contact_frequency_rows

        records = [
            {"conformer": "c1", "direct_target_contacts": [10], "near_target_contacts": []},
            {"conformer": "c1", "direct_target_contacts": [], "near_target_contacts": [10]},
            {"conformer": "c2", "direct_target_contacts": [10], "near_target_contacts": [20]},
        ]

        rows = build_contact_frequency_rows(
            records, [10, 20], ["c1", "c2"], recurrent_threshold=0.25
        )
        peptide_10 = next(
            row
            for row in rows
            if row["scope"] == "peptide" and row["target_residue_number"] == 10
        )

        self.assertAlmostEqual(2 / 3, peptide_10["direct_frequency"])
        self.assertAlmostEqual(0.75, peptide_10["conformer_balanced_direct_frequency"])
        self.assertAlmostEqual(0.25, peptide_10["conformer_balanced_near_frequency"])
        self.assertEqual(2, peptide_10["conformer_presence_count"])
        self.assertEqual(2, peptide_10["recurrent_conformer_support_count"])

    def test_presence_and_recurrent_support_are_independent(self) -> None:
        from peptide_ensemble_core import build_contact_frequency_rows

        records = []
        for conformer in ("c1", "c2", "c3"):
            records.extend(
                {
                    "conformer": conformer,
                    "direct_target_contacts": [10] if model == 1 else [],
                    "near_target_contacts": [],
                }
                for model in range(1, 6)
            )

        rows = build_contact_frequency_rows(
            records, [10], ["c1", "c2", "c3"], recurrent_threshold=0.25
        )
        peptide_row = next(row for row in rows if row["scope"] == "peptide")

        self.assertEqual(3, peptide_row["conformer_presence_count"])
        self.assertEqual(0, peptide_row["recurrent_conformer_support_count"])

    def test_shared_recurrent_core_uses_recurrent_support_not_presence(self) -> None:
        from peptide_ensemble_core import shared_recurrent_direct_core

        rows = [
            {
                "target_residue_number": 10,
                "conformer_presence_count": 3,
                "recurrent_conformer_support_count": 0,
            },
            {
                "target_residue_number": 20,
                "conformer_presence_count": 3,
                "recurrent_conformer_support_count": 3,
            },
        ]

        self.assertEqual([20], shared_recurrent_direct_core(rows, 3))


class FilterProvenanceTests(unittest.TestCase):
    def rules(self) -> dict[str, object]:
        return base_spec()["representative_selection_rules"]  # type: ignore[return-value]

    def test_all_filter_failures_are_retained(self) -> None:
        from peptide_ensemble_core import annotate_filter_provenance

        record = {
            "legacy_basic_geometry_pass": 0,
            "receptor_clash_pairs": 2,
            "chain_B_clash_pairs": 1,
            "direct_target_contacts": [],
            "legacy_5A_receptor_contact_fraction": 0.9,
        }

        result = annotate_filter_provenance(record, self.rules())

        self.assertFalse(result["receptor_clash_free"])
        self.assertFalse(result["chain_B_clash_free"])
        self.assertFalse(result["interface_eligible"])
        self.assertFalse(result["representative_eligible"])
        self.assertEqual("filtered_clash", result["selection_status"])
        self.assertEqual(
            [
                "receptor_clash_pairs>0",
                "chain_B_clash_pairs>0",
                "direct_target_contact_count<2",
            ],
            result["filter_reasons"],
        )

    def test_contact_fraction_is_not_a_hard_filter(self) -> None:
        from peptide_ensemble_core import annotate_filter_provenance

        record = {
            "legacy_basic_geometry_pass": 1,
            "receptor_clash_pairs": 0,
            "chain_B_clash_pairs": 0,
            "direct_target_contacts": [10, 20],
            "legacy_5A_receptor_contact_fraction": 0.1,
        }

        result = annotate_filter_provenance(record, self.rules())

        self.assertTrue(result["interface_eligible"])
        self.assertTrue(result["representative_eligible"])
        self.assertEqual("kept", result["selection_status"])
        self.assertEqual([], result["filter_reasons"])

    def test_legacy_basic_geometry_failure_does_not_control_eligibility(self) -> None:
        from peptide_ensemble_core import annotate_filter_provenance

        record = {
            "legacy_basic_geometry_pass": 0,
            "receptor_clash_pairs": 0,
            "chain_B_clash_pairs": 0,
            "direct_target_contacts": [10, 20],
            "legacy_5A_receptor_contact_fraction": 0.1,
        }

        result = annotate_filter_provenance(record, self.rules())

        self.assertTrue(result["receptor_clash_free"])
        self.assertTrue(result["chain_B_clash_free"])
        self.assertTrue(result["interface_eligible"])
        self.assertTrue(result["representative_eligible"])
        self.assertEqual("kept", result["selection_status"])
        self.assertNotIn("geometry_eligible", result)


class RepresentativeSelectionTests(unittest.TestCase):
    def test_pattern_prevalence_denominator_is_eligible_poses_only(self) -> None:
        from peptide_ensemble_core import select_contact_pattern_representatives

        records = [
            self.record("c1", 1, [10, 20], True, -6.0, 0.8),
            self.record("c1", 2, [10, 20], False, -7.0, 0.9),
            self.record("c1", 3, [10, 30], True, -5.9, 0.7),
            self.record("c1", 4, [10, 30], True, -5.8, 0.6),
        ]

        selected = select_contact_pattern_representatives(records, 2)
        representatives = [
            record for record in selected if record["selection_status"] == "representative"
        ]

        representatives.sort(key=lambda record: record["representative_rank_within_conformer"])
        self.assertEqual([3, 1], [record["model"] for record in representatives])
        by_model = {record["model"]: record for record in selected}
        self.assertAlmostEqual(1 / 3, by_model[1]["eligible_pattern_prevalence"])
        self.assertAlmostEqual(2 / 3, by_model[3]["eligible_pattern_prevalence"])
        self.assertAlmostEqual(2 / 4, by_model[1]["all_pose_pattern_prevalence"])

    @staticmethod
    def record(
        conformer: str,
        model: int,
        contacts: list[int],
        eligible: bool,
        score: float,
        contact_fraction: float,
    ) -> dict[str, object]:
        return {
            "pose_id": f"X_{conformer}_m{model}",
            "conformer": conformer,
            "model": model,
            "direct_target_contacts": contacts,
            "near_target_contacts": [],
            "representative_eligible": eligible,
            "selection_status": "kept" if eligible else "filtered_low_interface",
            "vina_score": score,
            "legacy_5A_receptor_contact_fraction": contact_fraction,
            "minimum_receptor_distance_A": 2.5,
        }


if __name__ == "__main__":
    unittest.main()
