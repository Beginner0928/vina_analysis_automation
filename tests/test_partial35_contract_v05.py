from __future__ import annotations

import copy
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT.parent
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from candidate_manifest_v05 import load_candidate_manifest  # noqa: E402
from competition_conformer_core import load_generation_protocol, sha256_file  # noqa: E402
from partial35_contract_v05 import (  # noqa: E402
    EXCLUDED_CANDIDATES,
    PARTIAL35_CANDIDATE_IDS,
    canonical_candidate_set_sha256,
    load_partial35_contract,
    partial_package_directory,
    validate_partial35_contract,
)


CONTRACT_PATH = ROOT / "specs" / "partial35_exploratory_package_v05.json"
PROTOCOL_PATH = ROOT / "specs" / "competition_conformer_protocol_v1.json"
MANIFEST_PATH = ROOT / "specs" / "top40_candidates_v05.json"


class Partial35ContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest, cls.candidates = load_candidate_manifest(MANIFEST_PATH)
        cls.protocol = load_generation_protocol(PROTOCOL_PATH)
        cls.protocol_sha256 = sha256_file(PROTOCOL_PATH)

    def test_frozen_exact_set_and_canonical_order_are_accepted(self) -> None:
        contract = load_partial35_contract(CONTRACT_PATH)
        selected = validate_partial35_contract(
            contract,
            self.candidates,
            self.protocol,
            self.protocol_sha256,
            canonical_manifest_sha256=sha256_file(MANIFEST_PATH),
        )
        self.assertEqual(list(PARTIAL35_CANDIDATE_IDS), [row["candidate_id"] for row in selected])
        self.assertEqual(35, len(selected))
        self.assertEqual(105, contract["expected_formal_sdf_count"])
        self.assertEqual(
            {
                "A8": "CONFORMER_GENERATION_QC_FAIL",
                "B1": "CONFORMER_GENERATION_QC_FAIL",
                "D3": "CONFORMER_GENERATION_QC_FAIL",
                "D7": "CONFORMER_GENERATION_QC_FAIL",
                "D10": "CONFORMER_GENERATION_QC_FAIL",
            },
            dict(EXCLUDED_CANDIDATES),
        )
        self.assertEqual("NOT_A_COMPLETE_TOP40_SCREEN", contract["completeness_statement"])

    def test_candidate_set_digest_covers_canonical_identity_rows(self) -> None:
        contract = load_partial35_contract(CONTRACT_PATH)
        selected = [
            row for row in self.candidates if row["candidate_id"] in set(PARTIAL35_CANDIDATE_IDS)
        ]
        digest = canonical_candidate_set_sha256(selected)
        self.assertEqual(contract["included_candidate_set_sha256"], digest)
        changed = copy.deepcopy(selected)
        changed[0]["sequence"] += "A"
        self.assertNotEqual(digest, canonical_candidate_set_sha256(changed))

    def test_unknown_duplicate_or_wrong_set_is_rejected(self) -> None:
        base = load_partial35_contract(CONTRACT_PATH)
        cases = []
        unknown = copy.deepcopy(base)
        unknown["included_candidate_ids"][-1] = "Z9"
        cases.append(unknown)
        duplicate = copy.deepcopy(base)
        duplicate["included_candidate_ids"][-1] = duplicate["included_candidate_ids"][0]
        cases.append(duplicate)
        incomplete = copy.deepcopy(base)
        incomplete["included_candidate_ids"].pop()
        incomplete["included_candidate_count"] = 34
        cases.append(incomplete)
        for contract in cases:
            with self.subTest(ids=contract["included_candidate_ids"][-2:]), self.assertRaises(ValueError):
                validate_partial35_contract(
                    contract,
                    self.candidates,
                    self.protocol,
                    self.protocol_sha256,
                    canonical_manifest_sha256=sha256_file(MANIFEST_PATH),
                )

    def test_partial_package_cannot_claim_top40_or_use_reserved_directory(self) -> None:
        contract = load_partial35_contract(CONTRACT_PATH)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            expected = (
                root
                / "standardized_ligands"
                / "partial_packages"
                / "competition_v1_partial35_mmff_qc_pass"
            )
            self.assertEqual(expected, partial_package_directory(root, contract, self.protocol))
            self.assertNotEqual(
                root / "standardized_ligands" / self.protocol["protocol_id"],
                expected,
            )
        changed = copy.deepcopy(contract)
        changed["completeness_statement"] = "COMPLETE_TOP40_SCREEN"
        with self.assertRaisesRegex(ValueError, "complete"):
            validate_partial35_contract(
                changed,
                self.candidates,
                self.protocol,
                self.protocol_sha256,
                canonical_manifest_sha256=sha256_file(MANIFEST_PATH),
            )

    def test_parent_protocol_science_and_hash_are_frozen(self) -> None:
        contract = load_partial35_contract(CONTRACT_PATH)
        self.assertEqual(
            "f79031dda94e148b851774b33082117ffde70f12b4cbf70b7f90a3266226a4a7",
            self.protocol_sha256,
        )
        self.assertEqual(self.protocol_sha256, contract["parent_generation_protocol"]["sha256"])
        self.assertEqual(30, self.protocol["embedding"]["requested_source_conformer_count"])
        self.assertEqual("MMFF94s", self.protocol["optimization"]["force_field"])
        self.assertEqual([2.0, 1.5, 1.0, 0.5], self.protocol["selection"]["rmsd_threshold_ladder_A"])

    def test_partial35_provenance_text_is_lf_and_portable(self) -> None:
        paths = (
            "specs/partial35_exploratory_package_v05.json",
            "specs/partial35_exploratory_package_v05.schema.json",
            "scripts/partial35_contract_v05.py",
            "scripts/generate_partial35_conformers.py",
            "scripts/audit_partial35_conformers.py",
        )
        for relative in paths:
            with self.subTest(path=relative):
                content = (ROOT / relative).read_bytes()
                self.assertNotIn(b"\r\n", content)
                self.assertNotRegex(content.decode("utf-8"), r"[A-Za-z]:\\")
                attributes = subprocess.check_output(
                    ["git", "check-attr", "text", "eol", "--", relative],
                    cwd=ROOT,
                    text=True,
                )
                self.assertIn(": text: set", attributes)
                self.assertIn(": eol: lf", attributes)

    def test_historical_a4_b3_c2_d8_sdf_hashes_still_match_frozen_inventory(self) -> None:
        inventory = json.loads(
            (ROOT / "specs" / "top40_ligand_inventory_v05.json").read_text(
                encoding="utf-8"
            )
        )
        selected = [
            row
            for row in inventory["slots"]
            if row["candidate_id"] in {"A4", "B3", "C2", "D8"}
            and row["readiness_status"] != "not_prepared"
        ]
        self.assertEqual(12, len(selected))
        for row in selected:
            path = DATA_ROOT / row["source_sdf_path"]
            self.assertEqual(
                row["source_sdf_sha256"],
                hashlib.sha256(path.read_bytes()).hexdigest(),
                f"historical SDF drift: {row['candidate_id']} {row['conformer_name']}",
            )


if __name__ == "__main__":
    unittest.main()
