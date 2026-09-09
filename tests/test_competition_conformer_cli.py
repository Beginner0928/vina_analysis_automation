from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from audit_competition_conformers import build_parser as build_audit_parser  # noqa: E402
from competition_conformer_core import (  # noqa: E402
    PROTOCOL_ID,
    audit_generation_package,
    audit_standardized_peptide,
    load_generation_protocol,
    sha256_file,
    write_generation_package,
)
from generate_competition_conformers import build_parser as build_production_parser  # noqa: E402
from smoke_competition_conformers import (  # noqa: E402
    SMOKE_CANDIDATES,
    build_parser as build_smoke_parser,
    main as smoke_main,
)


PROTOCOL_PATH = ROOT / "specs/competition_conformer_protocol_v1.json"
REQUIRED_SDF_FIELDS = {
    "CandidateID",
    "Group",
    "Sequence",
    "PeptideLength",
    "FormalCharge",
    "AtomCount",
    "HeavyAtomCount",
    "GenerationProtocolID",
    "GenerationProtocolSHA256",
    "RDKitVersion",
    "PythonVersion",
    "Platform",
    "GenerationSeed",
    "RequestedSourceConformerCount",
    "GeneratedSourceConformerCount",
    "MMFFConvergedConformerCount",
    "MMFFVariant",
    "MMFFMaxIterations",
    "MMFFConvergenceStatus",
    "SourceConformerID",
    "FormalConformerOrdinal",
    "MMFF94sEnergy",
    "SelectionAlgorithm",
    "SelectionRMSDThresholdUsed",
    "PairwiseRMSD_conf01_conf02_A",
    "PairwiseRMSD_conf01_conf03_A",
    "PairwiseRMSD_conf02_conf03_A",
    "ChemistryContractID",
    "ChemistryContractSHA256",
    "BondGraphHash",
    "StereochemistryHash",
}


def synthetic_generation_result() -> dict:
    candidate = {"candidate_id": "T1", "group": "T", "sequence": "AAA", "peptide_length": 3}
    molecule = __import__("competition_conformer_core").build_standardized_peptide("AAA")
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 1234
    parameters.numThreads = 1
    conformer_ids = list(AllChem.EmbedMultipleConfs(molecule, numConfs=3, params=parameters))
    if len(conformer_ids) != 3:
        raise AssertionError("Synthetic fixture did not generate exactly three conformers")
    chemistry_audit = audit_standardized_peptide(molecule, candidate)
    return {
        "molecule": molecule,
        "candidate": candidate,
        "chemistry_audit": chemistry_audit,
        "generation_seed": 1234,
        "requested_source_conformer_count": 30,
        "generated_source_conformer_count": 3,
        "mmff_converged_conformer_count": 3,
        "optimization_rows": [
            {"source_conformer_id": conformer_ids[0], "convergence_status": 0, "energy": -3.0},
            {"source_conformer_id": conformer_ids[1], "convergence_status": 0, "energy": -2.0},
            {"source_conformer_id": conformer_ids[2], "convergence_status": 0, "energy": -1.0},
        ],
        "selection": {
            "selected_source_conformer_ids": conformer_ids,
            "selected_energies": [-3.0, -2.0, -1.0],
            "rmsd_threshold_used_A": 1.0,
            "pairwise_rmsd_A": {
                "conf01-conf02": 1.25,
                "conf01-conf03": 1.50,
                "conf02-conf03": 1.75,
            },
            "eligible_conformer_count": 3,
        },
    }


class CompetitionConformerOutputTests(unittest.TestCase):
    def setUp(self) -> None:
        self.protocol = load_generation_protocol(PROTOCOL_PATH)
        self.chemistry_contract = {
            "contract_id": "v04_formal_standardized_peptide_chemistry"
        }

    def test_package_writes_required_sdf_provenance_and_external_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory) / PROTOCOL_ID
            manifest = write_generation_package(
                [synthetic_generation_result()],
                self.protocol,
                sha256_file(PROTOCOL_PATH),
                self.chemistry_contract,
                package,
                candidate_manifest_sha256="a" * 64,
                chemistry_contract_sha256="b" * 64,
            )

            self.assertEqual(package / "preparation_manifest.json", manifest)
            self.assertTrue((package / "generation_summary.tsv").is_file())
            self.assertTrue((package / "SHA256SUMS.tsv").is_file())
            saved = json.loads(manifest.read_text(encoding="utf-8"))
            self.assertEqual(PROTOCOL_ID, saved["protocol_id"])
            self.assertEqual(1, saved["candidate_count"])
            self.assertEqual(3, saved["formal_sdf_count"])
            hashes = {}
            for line in (package / "SHA256SUMS.tsv").read_text(encoding="utf-8").splitlines()[1:]:
                relative_path, digest = line.split("\t")
                hashes[relative_path] = digest
            self.assertEqual(3, len(hashes))

            for ordinal in range(1, 4):
                sdf = package / "T1" / f"T1_conf{ordinal:02d}.sdf"
                molecule = next(m for m in Chem.SDMolSupplier(str(sdf), removeHs=False) if m is not None)
                self.assertEqual(REQUIRED_SDF_FIELDS, set(molecule.GetPropNames()))
                self.assertEqual(f"conf{ordinal:02d}", molecule.GetProp("FormalConformerOrdinal"))
                self.assertFalse(molecule.HasProp("SDFSHA256"))
                self.assertEqual(sha256_file(sdf), hashes[sdf.relative_to(package).as_posix()])

    def test_existing_protocol_directory_is_never_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory) / PROTOCOL_ID
            package.mkdir()
            sentinel = package / "sentinel.txt"
            sentinel.write_text("preserve", encoding="utf-8")

            with self.assertRaisesRegex(FileExistsError, "Refusing to overwrite"):
                write_generation_package(
                    [synthetic_generation_result()],
                    self.protocol,
                    sha256_file(PROTOCOL_PATH),
                    self.chemistry_contract,
                    package,
                    candidate_manifest_sha256="a" * 64,
                    chemistry_contract_sha256="b" * 64,
                )

            self.assertEqual("preserve", sentinel.read_text(encoding="utf-8"))

    def test_audit_cli_accepts_intact_package_and_rejects_hash_drift(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory) / PROTOCOL_ID
            write_generation_package(
                [synthetic_generation_result()],
                self.protocol,
                sha256_file(PROTOCOL_PATH),
                self.chemistry_contract,
                package,
                candidate_manifest_sha256="a" * 64,
                chemistry_contract_sha256="b" * 64,
            )
            self.assertEqual("PASS", audit_generation_package(package)["status"])
            sdf = package / "T1/T1_conf01.sdf"
            with sdf.open("ab") as handle:
                handle.write(b"\n")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                audit_generation_package(package)

    def test_audit_rejects_manifest_scientific_provenance_drift(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory) / PROTOCOL_ID
            manifest_path = write_generation_package(
                [synthetic_generation_result()],
                self.protocol,
                sha256_file(PROTOCOL_PATH),
                self.chemistry_contract,
                package,
                candidate_manifest_sha256="a" * 64,
                chemistry_contract_sha256="b" * 64,
            )
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["candidates"][0]["formal_conformers"][0]["mmff94s_energy"] = 999.0
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "manifest/SDF provenance"):
                audit_generation_package(package)

    def test_audit_rejects_generation_summary_drift(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory) / PROTOCOL_ID
            write_generation_package(
                [synthetic_generation_result()],
                self.protocol,
                sha256_file(PROTOCOL_PATH),
                self.chemistry_contract,
                package,
                candidate_manifest_sha256="a" * 64,
                chemistry_contract_sha256="b" * 64,
            )
            summary = package / "generation_summary.tsv"
            summary.write_text(
                summary.read_text(encoding="utf-8").replace("\tPASS\n", "\tFAIL\n"),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "generation_summary"):
                audit_generation_package(package)


class CompetitionConformerCliSafetyTests(unittest.TestCase):
    def test_audit_cli_requires_protocol_data_root_and_package(self) -> None:
        parser = build_audit_parser()

        with self.assertRaises(SystemExit):
            parser.parse_args(["--protocol-directory", "X"])

    def test_production_cli_requires_explicit_all_selection(self) -> None:
        parser = build_production_parser()
        with self.assertRaises(SystemExit):
            parser.parse_args(["--protocol", str(PROTOCOL_PATH), "--data-root", "X"])

    def test_smoke_cli_has_fixed_a4_b3_candidate_set_and_external_root(self) -> None:
        self.assertEqual(("A4", "B3"), SMOKE_CANDIDATES)
        parser = build_smoke_parser()
        args = parser.parse_args(
            [
                "--protocol",
                str(PROTOCOL_PATH),
                "--data-root",
                "X",
                "--smoke-root",
                "Y",
            ]
        )
        self.assertEqual(Path("Y"), args.smoke_root)
        self.assertNotIn("candidates", vars(args))

    def test_smoke_cli_rejects_output_inside_git_repository_before_generation(self) -> None:
        output = ROOT / "test_output/forbidden_competition_smoke"

        exit_code = smoke_main(
            [
                "--protocol",
                str(PROTOCOL_PATH),
                "--data-root",
                str(ROOT.parent),
                "--smoke-root",
                str(output),
            ]
        )

        self.assertEqual(1, exit_code)
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
