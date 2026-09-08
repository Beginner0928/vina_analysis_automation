from __future__ import annotations

import hashlib
import json
import sys
import tomllib
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from environment_preflight_v05 import (  # noqa: E402
    audit_runtime_environment,
    parse_requirements_lock,
)


class EnvironmentContractTests(unittest.TestCase):
    def test_lock_contains_verified_windows_python312_runtime(self) -> None:
        lock_path = ROOT / "requirements-v05-lock.txt"
        locked = parse_requirements_lock(lock_path.read_text(encoding="utf-8"))
        self.assertEqual("2026.3.6", locked["rdkit"])
        self.assertEqual("0.8.0", locked["meeko"])
        self.assertEqual("2.6.1", locked["prody"])
        self.assertIn("numpy", locked)
        self.assertIn("scipy", locked)

    def test_environment_audit_rejects_wrong_python_or_package_version(self) -> None:
        lock = {"rdkit": "2026.3.6", "meeko": "0.8.0", "prody": "2.6.1"}
        installed = dict(lock)
        with self.assertRaisesRegex(ValueError, "Python 3.12"):
            audit_runtime_environment(lock, python_version=(3, 11, 9), installed_versions=installed)
        installed["rdkit"] = "2025.9.1"
        with self.assertRaisesRegex(ValueError, "rdkit"):
            audit_runtime_environment(lock, python_version=(3, 12, 5), installed_versions=installed)

    def test_environment_audit_passes_current_verified_versions(self) -> None:
        lock_path = ROOT / "requirements-v05-lock.txt"
        locked = parse_requirements_lock(lock_path.read_text(encoding="utf-8"))
        result = audit_runtime_environment(locked)
        self.assertEqual("PASS", result["status"])
        self.assertEqual("3.12", result["required_python_minor"])


class ProvenanceTextContractTests(unittest.TestCase):
    PROTECTED = [
        ROOT / "specs" / "top40_candidates_v05.json",
        ROOT / "specs" / "top40_candidates_v05.schema.json",
        ROOT / "specs" / "top40_ligand_inventory_v05.json",
        ROOT / "specs" / "top40_ligand_inventory_v05.schema.json",
        ROOT / "specs" / "screening_protocol_v05.json",
        ROOT / "specs" / "screening_protocol_v05.schema.json",
        ROOT / "requirements-v05-lock.txt",
        ROOT / "pyproject.toml",
    ]

    def test_new_contracts_are_lf_and_have_repository_eol_policy(self) -> None:
        attributes = (ROOT / ".gitattributes").read_text(encoding="utf-8")
        for rule in ("*.json text eol=lf", "*.txt text eol=lf", "*.toml text eol=lf", "*.py text eol=lf"):
            self.assertIn(rule, attributes)
        for path in self.PROTECTED:
            raw = path.read_bytes()
            self.assertNotIn(b"\r\n", raw, path.name)
            self.assertTrue(raw.endswith(b"\n"), path.name)
            self.assertEqual(64, len(hashlib.sha256(raw).hexdigest()))

    def test_contracts_and_schemas_are_valid_json_objects(self) -> None:
        for path in self.PROTECTED:
            if path.suffix != ".json":
                continue
            value = json.loads(path.read_text(encoding="utf-8"))
            self.assertIsInstance(value, dict, path.name)
        for schema in ROOT.glob("specs/*_v05.schema.json"):
            value = json.loads(schema.read_text(encoding="utf-8"))
            self.assertEqual("https://json-schema.org/draft/2020-12/schema", value["$schema"])

    def test_pyproject_keeps_pymol_out_of_production_dependencies(self) -> None:
        value = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        self.assertEqual(">=3.12,<3.13", value["project"]["requires-python"])
        dependencies = "\n".join(value["project"]["dependencies"]).lower()
        self.assertNotIn("pymol", dependencies)


if __name__ == "__main__":
    unittest.main()
