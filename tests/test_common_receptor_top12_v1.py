from __future__ import annotations

import copy
import hashlib
import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT.parent
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from plan_common_receptor_top12_v1 import (  # noqa: E402
    assert_new_output_root,
    load_contract,
    resolve_jobs,
    validate_jobs,
)


PROTOCOL_PATH = ROOT / "specs" / "common_receptor_top12_v1.json"
SELECTION_PATH = ROOT / "specs" / "top12_common_receptor_selection_v1.json"
INVENTORY_PATH = (
    DATA_ROOT
    / "docking_test_v1_top12"
    / "partial35_exploratory_v05"
    / "inputs"
    / "top40_ligand_inventory_partial35_v05.json"
)


class CommonReceptorTop12ContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.protocol = load_contract(PROTOCOL_PATH)
        cls.selection = load_contract(SELECTION_PATH)

    def test_protocol_freezes_one_complete_group_a_receptor_and_fixed_box(self) -> None:
        receptor = self.protocol["receptor"]
        self.assertEqual("A", receptor["registry_group"])
        self.assertEqual("h26d3_af3_s2718_m2", receptor["template_id"])
        self.assertEqual(
            "GPT_vina_branch_20260906/receptor_pdbqt/h26d3_af3_s2718_m2_TfR1_A.pdbqt",
            receptor["receptor_pdbqt"],
        )
        self.assertEqual(
            "f95e4f3eceb2774a711466c3d5690e14e256122f85d42df6dc459ea820f81ca1",
            receptor["sha256"]["receptor_pdbqt"],
        )
        self.assertEqual(
            [-10.827158, 23.888659, -11.674443], self.protocol["box"]["center_A"]
        )
        self.assertEqual(
            [30.591256, 37.457546, 33.838885], self.protocol["box"]["size_A"]
        )
        self.assertEqual(
            {
                "version": "1.2.7",
                "exhaustiveness": 16,
                "num_modes": 20,
                "energy_range": 5,
                "seed": 1701,
                "cpu": 8,
                "concurrency": 1,
            },
            {key: value for key, value in self.protocol["vina"].items() if key != "executable_sha256"},
        )
        self.assertEqual(
            "e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5",
            self.protocol["vina"]["executable_sha256"],
        )

    def test_selection_freezes_exact_sequences_order_and_three_conformers(self) -> None:
        expected = [
            ("A3", "PIPTHSDCPWGIHQ"),
            ("A9", "YPHDMIQEFIRP"),
            ("A4", "SPSNFITMYDW"),
            ("B3", "IPVGNCRYMQ"),
            ("B7", "KRVDSTDNCDFQM"),
            ("B4", "IKECPYFSQYGSM"),
            ("C2", "RDICQFQFHK"),
            ("C7", "KVFVKWPDSRY"),
            ("C8", "SKFGNPPSWGD"),
            ("D8", "ITPQVPWIRY"),
            ("D4", "GSYGYNDSIWCM"),
            ("D6", "PSATRDTYYTDTSRY"),
        ]
        self.assertEqual(
            expected,
            [(row["candidate_id"], row["sequence"]) for row in self.selection["candidates"]],
        )
        self.assertEqual(["conf01", "conf02", "conf03"], self.selection["conformers"])

    def test_checked_in_source_hashes_match_registry_and_selection(self) -> None:
        registry_path = ROOT / "specs" / "receptor_registry_v04.json"
        self.assertEqual(
            self.protocol["receptor_registry"]["sha256"],
            hashlib.sha256(registry_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            self.protocol["selection"]["sha256"],
            hashlib.sha256(SELECTION_PATH.read_bytes()).hexdigest(),
        )


class CommonReceptorTop12JobTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.protocol = load_contract(PROTOCOL_PATH)
        cls.selection = load_contract(SELECTION_PATH)
        cls.inventory = load_contract(INVENTORY_PATH)
        cls.jobs = resolve_jobs(
            protocol=cls.protocol,
            selection=cls.selection,
            inventory=cls.inventory,
            data_root=DATA_ROOT,
            output_root=DATA_ROOT / "docking_test" / "common_receptor_top12_v1",
        )

    def test_jobs_cover_12_by_3_and_have_one_receptor_box_and_vina_identity(self) -> None:
        audit = validate_jobs(self.jobs)
        self.assertEqual(36, audit["job_count"])
        self.assertEqual(12, audit["candidate_count"])
        self.assertEqual(3, audit["conformer_count"])
        self.assertEqual(1, audit["unique_receptor_pdbqt_count"])
        self.assertEqual(1, audit["unique_receptor_pdbqt_sha256_count"])
        self.assertEqual(1, audit["unique_center_count"])
        self.assertEqual(1, audit["unique_size_count"])
        self.assertEqual(1, audit["unique_vina_parameter_set_count"])
        self.assertEqual("PASS", audit["status"])

    def test_resolver_consumes_selection_order_instead_of_hard_coding_candidates(self) -> None:
        altered = copy.deepcopy(self.selection)
        altered["candidates"][0], altered["candidates"][1] = (
            altered["candidates"][1],
            altered["candidates"][0],
        )
        jobs = resolve_jobs(
            protocol=self.protocol,
            selection=altered,
            inventory=self.inventory,
            data_root=DATA_ROOT,
            output_root=DATA_ROOT / "docking_test" / "unused_common_receptor_test",
        )
        self.assertEqual(
            ["A9", "A9", "A9", "A3", "A3", "A3"],
            [row["candidate_id"] for row in jobs[:6]],
        )

    def test_resolver_rejects_ligand_hash_drift_and_duplicate_slot(self) -> None:
        drift = copy.deepcopy(self.inventory)
        target = next(row for row in drift["slots"] if row.get("candidate_id") == "A3")
        target["ligand_pdbqt_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "ligand.*SHA-256"):
            resolve_jobs(
                protocol=self.protocol,
                selection=self.selection,
                inventory=drift,
                data_root=DATA_ROOT,
                output_root=DATA_ROOT / "docking_test" / "unused_common_receptor_test",
            )
        duplicate = copy.deepcopy(self.inventory)
        duplicate["slots"].append(copy.deepcopy(duplicate["slots"][0]))
        with self.assertRaisesRegex(ValueError, "Duplicate ligand inventory slot"):
            resolve_jobs(
                protocol=self.protocol,
                selection=self.selection,
                inventory=duplicate,
                data_root=DATA_ROOT,
                output_root=DATA_ROOT / "docking_test" / "unused_common_receptor_test",
            )

    def test_new_output_root_guard_rejects_any_existing_path(self) -> None:
        with self.assertRaises(FileExistsError):
            assert_new_output_root(ROOT)
        assert_new_output_root(DATA_ROOT / "docking_test" / "definitely_absent_common_receptor_test")


if __name__ == "__main__":
    unittest.main()
