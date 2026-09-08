from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "scripts"))

from peptide_ensemble_core import validate_ensemble_spec  # noqa: E402


class V04SchemaContractTests(unittest.TestCase):
    def test_repository_lf_policy_protects_provenance_critical_text(self) -> None:
        attributes_path = REPOSITORY_ROOT / ".gitattributes"
        if not attributes_path.is_file():
            self.fail("repository .gitattributes is missing")

        required_rules = {
            ".gitattributes text eol=lf",
            ".gitignore text eol=lf",
            "*.json text eol=lf",
            "*.py text eol=lf",
            "*.md text eol=lf",
            "*.tsv text eol=lf",
            "*.csv text eol=lf",
            "*.txt text eol=lf",
            "*.yml text eol=lf",
            "*.yaml text eol=lf",
            "*.toml text eol=lf",
        }
        actual_rules = {
            line.strip()
            for line in attributes_path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        }
        self.assertTrue(required_rules.issubset(actual_rules))

        registry_relative = "specs/receptor_registry_v04.json"
        attributes = subprocess.run(
            ["git", "check-attr", "text", "eol", "--", registry_relative],
            cwd=REPOSITORY_ROOT,
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        self.assertIn("text: set", attributes)
        self.assertIn("eol: lf", attributes)

        registry_bytes = (REPOSITORY_ROOT / registry_relative).read_bytes()
        self.assertNotIn(b"\r\n", registry_bytes)
        locked_specs = (
            "A4_ensemble_v04.json",
            "A4_formal_ensemble_v04.json",
            "A4_standardized_docking_v04.json",
            "B3_ensemble_v04.json",
            "B3_standardized_docking_v04.json",
        )
        expected_hashes = {
            json.loads((REPOSITORY_ROOT / "specs" / name).read_text(encoding="utf-8"))[
                "receptor_registry_sha256"
            ].lower()
            for name in locked_specs
        }
        self.assertEqual(1, len(expected_hashes))
        self.assertEqual(
            expected_hashes.pop(),
            hashlib.sha256(registry_bytes).hexdigest(),
        )

    def test_versioned_ensemble_specs_match_declared_top_level_contract(self) -> None:
        schema = json.loads(
            (REPOSITORY_ROOT / "specs/peptide_ensemble_v04.schema.json").read_text(
                encoding="utf-8"
            )
        )
        allowed = set(schema["properties"])
        required = set(schema["required"])
        for name in ("A4_ensemble_v04.json", "B3_ensemble_v04.json"):
            with self.subTest(spec=name):
                spec = json.loads(
                    (REPOSITORY_ROOT / "specs" / name).read_text(encoding="utf-8")
                )
                self.assertFalse(set(spec) - allowed)
                self.assertFalse(required - set(spec))
                validate_ensemble_spec(spec)

    def test_registry_is_the_only_standardized_receptor_mapping_source(self) -> None:
        spec = json.loads(
            (
                REPOSITORY_ROOT / "specs/B3_standardized_docking_v04.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(
            "vina_analysis_automation/specs/receptor_registry_v04.json",
            spec["receptor_registry_path"],
        )
        self.assertEqual(64, len(spec["receptor_registry_sha256"]))
        self.assertNotIn("receptor_sources", spec)
        self.assertNotIn("box_registry_path", spec)
        self.assertNotIn("box_registry_sha256", spec)


if __name__ == "__main__":
    unittest.main()
