from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "scripts"))

from peptide_ensemble_core import validate_ensemble_spec  # noqa: E402


class V04SchemaContractTests(unittest.TestCase):
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
