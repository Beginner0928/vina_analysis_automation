from __future__ import annotations

import copy
import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from candidate_manifest_v05 import (  # noqa: E402
    load_json_contract,
    select_candidates,
    validate_candidate_manifest,
    validate_ligand_inventory,
    validate_screening_protocol,
)


class CandidateManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads(
            (ROOT / "specs" / "top40_candidates_v05.json").read_text(encoding="utf-8")
        )

    def test_canonical_inventory_is_exactly_a1_to_d10(self) -> None:
        candidates = validate_candidate_manifest(copy.deepcopy(self.manifest))
        expected = [f"{group}{index}" for group in "ABCD" for index in range(1, 11)]
        self.assertEqual(expected, [row["candidate_id"] for row in candidates])
        self.assertEqual(40, len(candidates))
        self.assertEqual({"A": 10, "B": 10, "C": 10, "D": 10}, self.manifest["expected_group_counts"])

    def test_group_selection_is_canonical_and_deterministic(self) -> None:
        candidates = validate_candidate_manifest(copy.deepcopy(self.manifest))
        ab = select_candidates(candidates, groups=["B", "A"])
        cd = select_candidates(candidates, groups=["D", "C"])
        self.assertEqual([f"A{i}" for i in range(1, 11)] + [f"B{i}" for i in range(1, 11)], [row["candidate_id"] for row in ab])
        self.assertEqual([f"C{i}" for i in range(1, 11)] + [f"D{i}" for i in range(1, 11)], [row["candidate_id"] for row in cd])

    def test_explicit_subset_uses_manifest_order_not_cli_order(self) -> None:
        candidates = validate_candidate_manifest(copy.deepcopy(self.manifest))
        selected = select_candidates(candidates, candidate_ids=["B3", "A2", "A1"])
        self.assertEqual(["A1", "A2", "B3"], [row["candidate_id"] for row in selected])

    def test_selector_conflict_and_missing_selector_are_rejected(self) -> None:
        candidates = validate_candidate_manifest(copy.deepcopy(self.manifest))
        with self.assertRaisesRegex(ValueError, "exactly one"):
            select_candidates(candidates, groups=["A"], candidate_ids=["A1"])
        with self.assertRaisesRegex(ValueError, "exactly one"):
            select_candidates(candidates)

    def test_duplicate_requested_candidate_and_unknown_candidate_are_rejected(self) -> None:
        candidates = validate_candidate_manifest(copy.deepcopy(self.manifest))
        with self.assertRaisesRegex(ValueError, "duplicate"):
            select_candidates(candidates, candidate_ids=["A1", "A1"])
        with self.assertRaisesRegex(ValueError, "unknown candidate"):
            select_candidates(candidates, candidate_ids=["A1", "Z9"])

    def test_duplicate_manifest_candidate_is_rejected(self) -> None:
        manifest = copy.deepcopy(self.manifest)
        manifest["candidates"][1] = copy.deepcopy(manifest["candidates"][0])
        with self.assertRaisesRegex(ValueError, "duplicate candidate_id"):
            validate_candidate_manifest(manifest)

    def test_wrong_group_mapping_and_invalid_group_are_rejected(self) -> None:
        wrong_mapping = copy.deepcopy(self.manifest)
        wrong_mapping["candidates"][0]["group"] = "B"
        with self.assertRaisesRegex(ValueError, "ID-prefix/group"):
            validate_candidate_manifest(wrong_mapping)

        invalid_group = copy.deepcopy(self.manifest)
        invalid_group["candidates"][0]["group"] = "Z"
        with self.assertRaisesRegex(ValueError, "invalid group"):
            validate_candidate_manifest(invalid_group)

    def test_invalid_sequence_and_length_mismatch_are_rejected(self) -> None:
        invalid_sequence = copy.deepcopy(self.manifest)
        invalid_sequence["candidates"][0]["sequence"] = "ACDZ"
        invalid_sequence["candidates"][0]["peptide_length"] = 4
        with self.assertRaisesRegex(ValueError, "standard amino-acid"):
            validate_candidate_manifest(invalid_sequence)

        wrong_length = copy.deepcopy(self.manifest)
        wrong_length["candidates"][0]["peptide_length"] += 1
        with self.assertRaisesRegex(ValueError, "peptide_length"):
            validate_candidate_manifest(wrong_length)

    def test_source_provenance_values_follow_additional_properties_schema(self) -> None:
        malformed = copy.deepcopy(self.manifest)
        malformed["source_provenance"]["sources"]["top40_candidates.tsv"] = "not-a-hash"
        schema = json.loads(
            (ROOT / "specs" / "top40_candidates_v05.schema.json").read_text(
                encoding="utf-8"
            )
        )
        from candidate_manifest_v05 import validate_json_schema

        with self.assertRaisesRegex(ValueError, "does not match"):
            validate_json_schema(malformed, schema)

    def test_missing_or_extra_formal_candidate_is_rejected(self) -> None:
        missing = copy.deepcopy(self.manifest)
        missing["candidates"].pop()
        with self.assertRaisesRegex(ValueError, "exactly 40"):
            validate_candidate_manifest(missing)

        extra = copy.deepcopy(self.manifest)
        extra["candidates"].append(
            {"candidate_id": "A11", "group": "A", "sequence": "AAA", "peptide_length": 3}
        )
        with self.assertRaisesRegex(ValueError, "exactly 40"):
            validate_candidate_manifest(extra)


class LigandInventoryAndProtocolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads((ROOT / "specs" / "top40_candidates_v05.json").read_text(encoding="utf-8"))
        cls.inventory = json.loads((ROOT / "specs" / "top40_ligand_inventory_v05.json").read_text(encoding="utf-8"))
        cls.protocol = json.loads((ROOT / "specs" / "screening_protocol_v05.json").read_text(encoding="utf-8"))

    def test_inventory_has_all_120_slots_without_fabricating_readiness(self) -> None:
        candidates = validate_candidate_manifest(copy.deepcopy(self.manifest))
        protocol = validate_screening_protocol(copy.deepcopy(self.protocol))
        slots = validate_ligand_inventory(copy.deepcopy(self.inventory), candidates, protocol)
        counts: dict[str, int] = {}
        for slot in slots.values():
            counts[slot["readiness_status"]] = counts.get(slot["readiness_status"], 0) + 1
        self.assertEqual(120, len(slots))
        self.assertEqual(
            {"approved_for_formal_screening": 6, "conversion_audit_only": 6, "not_prepared": 108},
            counts,
        )

    def test_required_conformer_count_is_protocol_driven_not_hard_coded(self) -> None:
        candidates = validate_candidate_manifest(copy.deepcopy(self.manifest))
        protocol = copy.deepcopy(self.protocol)
        protocol["formal_conformers"] = {
            "required_conformer_count": 2,
            "required_conformer_names": ["shape_a", "shape_b"],
        }
        protocol = validate_screening_protocol(protocol, enforce_frozen_values=False)
        inventory = {
            "schema_version": "0.5",
            "inventory_id": "synthetic_two_conformer_inventory",
            "expected_slot_count": 80,
            "slots": [
                {
                    "candidate_id": candidate["candidate_id"],
                    "conformer_name": conformer,
                    "readiness_status": "not_prepared",
                }
                for candidate in candidates
                for conformer in ("shape_a", "shape_b")
            ],
        }
        slots = validate_ligand_inventory(inventory, candidates, protocol)
        self.assertEqual(80, len(slots))

    def test_inventory_rejects_missing_slot_and_duplicate_slot(self) -> None:
        candidates = validate_candidate_manifest(copy.deepcopy(self.manifest))
        protocol = validate_screening_protocol(copy.deepcopy(self.protocol))

        missing = copy.deepcopy(self.inventory)
        missing["slots"].pop()
        with self.assertRaisesRegex(ValueError, "slot inventory"):
            validate_ligand_inventory(missing, candidates, protocol)

        duplicate = copy.deepcopy(self.inventory)
        duplicate["slots"][-1] = copy.deepcopy(duplicate["slots"][0])
        with self.assertRaisesRegex(ValueError, "duplicate ligand slot"):
            validate_ligand_inventory(duplicate, candidates, protocol)

    def test_inventory_schema_rejects_unexpected_slot_fields(self) -> None:
        malformed = copy.deepcopy(self.inventory)
        malformed["slots"][0]["fabricated_ready_flag"] = True
        schema = json.loads(
            (ROOT / "specs" / "top40_ligand_inventory_v05.schema.json").read_text(
                encoding="utf-8"
            )
        )
        from candidate_manifest_v05 import validate_json_schema

        with self.assertRaisesRegex(ValueError, "unexpected fields"):
            validate_json_schema(malformed, schema)

    def test_frozen_protocol_values_are_exact(self) -> None:
        protocol = validate_screening_protocol(copy.deepcopy(self.protocol))
        self.assertEqual(3, protocol["formal_conformers"]["required_conformer_count"])
        self.assertEqual(["conf01", "conf02", "conf03"], protocol["formal_conformers"]["required_conformer_names"])
        self.assertEqual(
            {
                "engine": "AutoDock Vina",
                "version": "1.2.7",
                "scoring_function": "vina",
                "exhaustiveness": 32,
                "num_modes": 20,
                "energy_range": 5,
                "seed": 1701,
                "cpu": 8,
            },
            protocol["docking_protocol"],
        )
        self.assertEqual(1, protocol["batch_execution"]["vina_process_concurrency"])

    def test_frozen_ligand_preparation_values_are_exact(self) -> None:
        protocol = copy.deepcopy(self.protocol)
        protocol["ligand_preparation"]["source_coordinates"] = "regenerated_coordinates"
        with self.assertRaisesRegex(ValueError, "ligand preparation"):
            validate_screening_protocol(protocol)


if __name__ == "__main__":
    unittest.main()
