from __future__ import annotations

import copy
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from competition_conformer_core import (  # noqa: E402
    PROTOCOL_ID,
    load_generation_protocol,
    protocol_identity_sha256,
    validate_generation_protocol,
)


PROTOCOL_PATH = ROOT / "specs/competition_conformer_protocol_v1.json"


class CompetitionConformerContractTests(unittest.TestCase):
    def test_checked_in_contract_freezes_competition_v1(self) -> None:
        protocol = load_generation_protocol(PROTOCOL_PATH)

        self.assertEqual(PROTOCOL_ID, protocol["protocol_id"])
        self.assertEqual("COMPETITION_EXPLORATORY_SCREENING", protocol["classification"])
        self.assertEqual("2026.03.6", protocol["software"]["rdkit_version"])
        self.assertEqual(
            {
                "method": "ETKDGv3",
                "requested_source_conformer_count": 30,
                "use_random_coords": True,
                "prune_rms_threshold_A": 0.25,
                "enforce_chirality": True,
                "num_threads": 1,
            },
            protocol["embedding"],
        )
        self.assertEqual("MMFF94s", protocol["optimization"]["force_field"])
        self.assertEqual(1000, protocol["optimization"]["max_iterations"])
        self.assertEqual(0, protocol["optimization"]["eligible_convergence_status"])
        self.assertEqual("none", protocol["optimization"]["fallback"])
        self.assertEqual([2.0, 1.5, 1.0, 0.5], protocol["selection"]["rmsd_threshold_ladder_A"])
        self.assertEqual(3, protocol["selection"]["final_conformer_count"])
        self.assertEqual("rdMolAlign.GetAlignmentTransform", protocol["rmsd"]["rdkit_api"])

    def test_parameter_drift_is_rejected(self) -> None:
        protocol = load_generation_protocol(PROTOCOL_PATH)
        changed = copy.deepcopy(protocol)
        changed["embedding"]["num_threads"] = 0

        with self.assertRaisesRegex(ValueError, "embedding"):
            validate_generation_protocol(changed)

    def test_identity_seed_selection_fail_closed_and_output_drift_are_rejected(self) -> None:
        protocol = load_generation_protocol(PROTOCOL_PATH)
        mutations = [
            ("candidate authority", ("candidate_authority", "path"), "other.json"),
            ("chemistry contract", ("chemistry_contract", "contract_id"), "other"),
            ("seed", ("seed", "algorithm"), "automatic_random"),
            ("selection", ("selection", "algorithm"), "other"),
            ("fail_closed", ("fail_closed", "fewer_than_three_converged"), False),
            ("output", ("output", "sdf_filename"), "{candidate_id}.sdf"),
        ]

        for label, (section, field), value in mutations:
            with self.subTest(label=label):
                changed = copy.deepcopy(protocol)
                changed[section][field] = value
                with self.assertRaises(ValueError):
                    validate_generation_protocol(changed)

    def test_removing_required_sdf_provenance_field_is_rejected(self) -> None:
        protocol = load_generation_protocol(PROTOCOL_PATH)
        changed = copy.deepcopy(protocol)
        changed["sdf_provenance_fields"].remove("SourceConformerID")

        with self.assertRaisesRegex(ValueError, "provenance"):
            validate_generation_protocol(changed)

    def test_contract_paths_are_portable(self) -> None:
        protocol = load_generation_protocol(PROTOCOL_PATH)

        for contract in (protocol["candidate_authority"], protocol["chemistry_contract"]):
            self.assertFalse(Path(contract["path"]).is_absolute())
            self.assertNotIn("E:\\", contract["path"])

    def test_protocol_byte_change_changes_identity_hash(self) -> None:
        original = PROTOCOL_PATH.read_bytes()
        changed = original.replace(b'"MMFF94s"', b'"MMFF94"', 1)

        self.assertNotEqual(
            protocol_identity_sha256(original),
            protocol_identity_sha256(changed),
        )

    def test_new_protocol_artifacts_are_lf_and_covered_by_repository_policy(self) -> None:
        paths = [
            "specs/competition_conformer_protocol_v1.json",
            "specs/competition_conformer_protocol_v1.schema.json",
            "scripts/competition_conformer_core.py",
            "docs/v0.5_competition_conformers.md",
        ]
        for relative in paths:
            with self.subTest(path=relative):
                self.assertNotIn(b"\r\n", (ROOT / relative).read_bytes())
                attributes = subprocess.check_output(
                    ["git", "check-attr", "text", "eol", "--", relative],
                    cwd=ROOT,
                    text=True,
                )
                self.assertIn(": text: set", attributes)
                self.assertIn(": eol: lf", attributes)

    def test_git_tracks_no_raw_science_or_generated_output(self) -> None:
        tracked = subprocess.check_output(
            ["git", "ls-files"], cwd=ROOT, text=True
        ).splitlines()
        forbidden_suffixes = (".pdb", ".pdbqt", ".sdf", ".png", ".pse")
        forbidden_directories = ("docking_test/", "test_output/", "GPT_vina_branch_20260906/")

        self.assertEqual(
            [],
            [
                path
                for path in tracked
                if path.lower().endswith(forbidden_suffixes)
                or any(path.startswith(directory) for directory in forbidden_directories)
            ],
        )


if __name__ == "__main__":
    unittest.main()
