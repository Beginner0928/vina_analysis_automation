from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT.parent
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from candidate_manifest_v05 import load_candidate_manifest  # noqa: E402
from competition_conformer_core import load_generation_protocol, sha256_file  # noqa: E402
from generate_partial35_conformers import (  # noqa: E402
    build_parser as build_generation_parser,
    generate_partial35_package,
)
from partial35_contract_v05 import (  # noqa: E402
    PARTIAL35_CANDIDATE_IDS,
    load_partial35_contract,
)


CONTRACT_PATH = ROOT / "specs" / "partial35_exploratory_package_v05.json"
PROTOCOL_PATH = ROOT / "specs" / "competition_conformer_protocol_v1.json"
MANIFEST_PATH = ROOT / "specs" / "top40_candidates_v05.json"


def fake_result(candidate: dict) -> dict:
    return {
        "candidate": candidate,
        "generation_seed": 1,
        "generated_source_conformer_count": 30,
        "mmff_converged_conformer_count": 9,
        "selection": {
            "selected_source_conformer_ids": [1, 2, 3],
            "selected_energies": [1.0, 2.0, 3.0],
            "rmsd_threshold_used_A": 2.0,
            "pairwise_rmsd_A": {
                "conf01-conf02": 2.1,
                "conf01-conf03": 2.2,
                "conf02-conf03": 2.3,
            },
        },
    }


class Partial35ConformerCliTests(unittest.TestCase):
    def test_cli_requires_contract_protocol_and_data_root_only(self) -> None:
        parser = build_generation_parser()
        with self.assertRaises(SystemExit):
            parser.parse_args([])
        parsed = parser.parse_args(
            ["--contract", str(CONTRACT_PATH), "--protocol", str(PROTOCOL_PATH), "--data-root", str(DATA_ROOT)]
        )
        self.assertEqual(CONTRACT_PATH, parsed.contract)
        self.assertFalse(hasattr(parsed, "force"))
        self.assertFalse(hasattr(parsed, "all"))

    def test_generation_is_canonical_sequential_and_records_partial_provenance(self) -> None:
        contract = load_partial35_contract(CONTRACT_PATH)
        protocol = load_generation_protocol(PROTOCOL_PATH)
        _, candidates = load_candidate_manifest(MANIFEST_PATH)
        calls: list[str] = []

        def generate(candidate, _protocol):
            calls.append(candidate["candidate_id"])
            return fake_result(candidate)

        def write(results, _protocol, _protocol_sha, _chemistry, package, **_kwargs):
            package.mkdir(parents=True)
            manifest = package / "preparation_manifest.json"
            manifest.write_text(
                json.dumps({"candidate_count": len(results), "formal_sdf_count": len(results) * 3}) + "\n",
                encoding="utf-8",
                newline="\n",
            )
            return manifest

        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory) / "partial"
            manifest = generate_partial35_package(
                contract=contract,
                protocol=protocol,
                protocol_sha256=sha256_file(PROTOCOL_PATH),
                canonical_candidates=candidates,
                canonical_manifest_sha256=sha256_file(MANIFEST_PATH),
                chemistry_contract={"contract_id": "v04_formal_standardized_peptide_chemistry"},
                chemistry_contract_sha256=contract["chemistry_contract"]["sha256"],
                package_directory=package,
                generator_git_commit="a" * 40,
                generate_fn=generate,
                write_fn=write,
            )
            self.assertEqual(list(PARTIAL35_CANDIDATE_IDS), calls)
            payload = json.loads(manifest.read_text(encoding="utf-8"))
            self.assertEqual(35, payload["included_candidate_count"])
            self.assertEqual(105, payload["expected_formal_sdf_count"])
            self.assertEqual("NOT_A_COMPLETE_TOP40_SCREEN", payload["completeness_statement"])
            self.assertEqual("a" * 40, payload["generator_git_commit"])
            self.assertEqual(list(PARTIAL35_CANDIDATE_IDS), payload["included_candidate_ids"])
            self.assertEqual(5, len(payload["excluded_candidates"]))
            self.assertRegex(payload["core_preparation_manifest_sha256"], r"^[0-9a-f]{64}$")

    def test_existing_package_is_never_overwritten(self) -> None:
        contract = load_partial35_contract(CONTRACT_PATH)
        protocol = load_generation_protocol(PROTOCOL_PATH)
        _, candidates = load_candidate_manifest(MANIFEST_PATH)
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory) / "existing"
            package.mkdir()
            sentinel = package / "sentinel.txt"
            sentinel.write_text("keep", encoding="utf-8")
            with self.assertRaisesRegex(FileExistsError, "overwrite"):
                generate_partial35_package(
                    contract=contract,
                    protocol=protocol,
                    protocol_sha256=sha256_file(PROTOCOL_PATH),
                    canonical_candidates=candidates,
                    canonical_manifest_sha256=sha256_file(MANIFEST_PATH),
                    chemistry_contract={"contract_id": "v04_formal_standardized_peptide_chemistry"},
                    chemistry_contract_sha256=contract["chemistry_contract"]["sha256"],
                    package_directory=package,
                    generator_git_commit="b" * 40,
                )
            self.assertEqual("keep", sentinel.read_text(encoding="utf-8"))

    def test_generation_subsystem_has_no_vina_or_process_invocation(self) -> None:
        paths = [
            SCRIPTS / "partial35_contract_v05.py",
            SCRIPTS / "generate_partial35_conformers.py",
            SCRIPTS / "audit_partial35_conformers.py",
        ]
        combined = "\n".join(path.read_text(encoding="utf-8") for path in paths)
        self.assertNotIn("subprocess", combined)
        self.assertNotIn("run_standardized_vina", combined)
        self.assertNotIn("AutoDock", combined)


if __name__ == "__main__":
    unittest.main()
