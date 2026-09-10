from __future__ import annotations

import copy
import hashlib
import json
import sys
import tempfile
import unittest
import uuid
from pathlib import Path
from unittest.mock import patch

from rdkit import Chem
from rdkit.Chem import AllChem


ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT.parent
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from candidate_manifest_v05 import (  # noqa: E402
    load_candidate_manifest,
    validate_ligand_inventory,
    validate_screening_protocol,
)
from competition_conformer_core import build_standardized_peptide  # noqa: E402
from partial35_contract_v05 import (  # noqa: E402
    EXCLUDED_CANDIDATES,
    PARTIAL35_CANDIDATE_IDS,
)
from prepare_partial35_ligands import (  # noqa: E402
    build_partial35_inventory,
    build_partial35_screening_protocol,
)
from prepare_standardized_vina import prepare_one  # noqa: E402
from standardized_vina_inputs import audit_source_sdf, sha256  # noqa: E402
from batch_preflight_v05 import run_dry_run  # noqa: E402
from tests.test_batch_preflight_v05 import unused_output, write_contract_fixture  # noqa: E402
from audit_partial35_ligands import validate_partial35_inventory_contract  # noqa: E402


MANIFEST_PATH = ROOT / "specs" / "top40_candidates_v05.json"
PROTOCOL_PATH = ROOT / "specs" / "screening_protocol_v05.json"
B3_SPEC = ROOT / "specs" / "B3_standardized_docking_v04.json"
B3_ROOT = DATA_ROOT / "docking_test" / "vina_standardized_v04_20260907" / "B3"
VINA = DATA_ROOT / "tools" / "AutoDockVina" / "vina.exe"


class NTerminalProAuditTests(unittest.TestCase):
    def test_protonated_n_terminal_proline_uses_two_hydrogens(self) -> None:
        molecule = build_standardized_peptide("PAA")
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = 123
        self.assertEqual(0, AllChem.EmbedMolecule(molecule, parameters))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "P1_conf01.sdf"
            molecule.SetProp("_Name", path.stem)
            molecule.SetProp("Sequence", "PAA")
            molecule.SetProp("FormalCharge", "0")
            writer = Chem.SDWriter(str(path))
            writer.write(molecule)
            writer.close()
            audit = audit_source_sdf(path, "P1", "PAA", 0)
        self.assertEqual(2, audit["chemistry_checks"]["n_terminus_hydrogen_count"])
        self.assertEqual("PASS", audit["chemistry_checks"]["n_terminus_rule_status"])


class Partial35InventoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        _, cls.candidates = load_candidate_manifest(MANIFEST_PATH)
        cls.base_protocol = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))

    def artifact_records(self) -> list[dict]:
        records = []
        for candidate_id in PARTIAL35_CANDIDATE_IDS:
            for conformer in ("conf01", "conf02", "conf03"):
                stem = f"{candidate_id}_{conformer}"
                records.append(
                    {
                        "candidate_id": candidate_id,
                        "conformer_name": conformer,
                        "source_sdf_path": f"standardized_ligands/partial_packages/p/{candidate_id}/{stem}.sdf",
                        "source_sdf_sha256": "1" * 64,
                        "audit_path": f"docking_test/partial35/inputs/audit/{stem}.json",
                        "audit_sha256": "2" * 64,
                        "ligand_pdbqt_path": f"docking_test/partial35/inputs/ligands/{candidate_id}/{stem}.pdbqt",
                        "ligand_pdbqt_sha256": "3" * 64,
                        "pdbqt_audit_path": f"docking_test/partial35/inputs/audit/{stem}.json",
                        "pdbqt_audit_sha256": "2" * 64,
                        "generation_provenance": {
                            "protocol_id": "tfr1_top40_sdfgen_competition_v1_b3derived",
                            "protocol_sha256": "4" * 64,
                            "partial_package_id": "competition_v1_partial35_mmff_qc_pass",
                            "partial_package_manifest_sha256": "5" * 64,
                        },
                    }
                )
        return records

    def test_derived_inventory_is_full_canonical_slots_with_exact_partial_statuses(self) -> None:
        inventory = build_partial35_inventory(
            self.candidates,
            self.artifact_records(),
            candidate_manifest_sha256=hashlib.sha256(MANIFEST_PATH.read_bytes()).hexdigest(),
            partial_package_manifest_sha256="5" * 64,
            partial_package_manifest_path="standardized_ligands/partial_packages/p/partial_package_manifest.json",
            ligand_preparation_contract_sha256="6" * 64,
        )
        self.assertEqual(120, inventory["expected_slot_count"])
        self.assertEqual(105, inventory["approved_exploratory_slot_count"])
        self.assertEqual(15, inventory["qc_excluded_slot_count"])
        slots = inventory["slots"]
        self.assertEqual(
            [(row["candidate_id"], row["conformer_name"]) for row in slots],
            [
                (candidate["candidate_id"], conformer)
                for candidate in self.candidates
                for conformer in ("conf01", "conf02", "conf03")
            ],
        )
        approved = [row for row in slots if row["readiness_status"] == "approved_for_exploratory_screening"]
        excluded = [row for row in slots if row["readiness_status"] == "conformer_generation_qc_fail_not_docked"]
        self.assertEqual(105, len(approved))
        self.assertEqual(15, len(excluded))
        self.assertEqual(
            {candidate_id for candidate_id, _reason in EXCLUDED_CANDIDATES},
            {row["candidate_id"] for row in excluded},
        )
        validate_ligand_inventory(inventory, self.candidates, self.base_protocol)

    def test_derived_protocol_preserves_frozen_science_and_locks_inventory(self) -> None:
        inventory_path = Path("top40_ligand_inventory_partial35_v05.json")
        derived = build_partial35_screening_protocol(
            self.base_protocol,
            inventory_path=inventory_path,
            inventory_sha256="a" * 64,
        )
        self.assertEqual("tfr1_partial35_exploratory_screening_v05", derived["protocol_id"])
        self.assertEqual(self.base_protocol["ligand_preparation"], derived["ligand_preparation"])
        self.assertEqual(self.base_protocol["docking_protocol"], derived["docking_protocol"])
        self.assertEqual(self.base_protocol["vina_executable"], derived["vina_executable"])
        self.assertEqual(self.base_protocol["batch_execution"], derived["batch_execution"])
        self.assertEqual(inventory_path.as_posix(), derived["contracts"]["ligand_inventory"]["path"])
        self.assertEqual("a" * 64, derived["contracts"]["ligand_inventory"]["sha256"])
        validate_screening_protocol(derived)

    def test_exploratory_inventory_rejects_missing_pdbqt_provenance(self) -> None:
        inventory = build_partial35_inventory(
            self.candidates,
            self.artifact_records(),
            candidate_manifest_sha256=hashlib.sha256(MANIFEST_PATH.read_bytes()).hexdigest(),
            partial_package_manifest_sha256="5" * 64,
            partial_package_manifest_path="standardized_ligands/partial_packages/p/partial_package_manifest.json",
            ligand_preparation_contract_sha256="6" * 64,
        )
        del inventory["slots"][0]["ligand_pdbqt_sha256"]
        with self.assertRaisesRegex(ValueError, "PDBQT"):
            validate_ligand_inventory(inventory, self.candidates, self.base_protocol)

    def test_gate_c_contract_rejects_wrong_count_or_excluded_set(self) -> None:
        inventory = build_partial35_inventory(
            self.candidates,
            self.artifact_records(),
            candidate_manifest_sha256=hashlib.sha256(MANIFEST_PATH.read_bytes()).hexdigest(),
            partial_package_manifest_sha256="5" * 64,
            partial_package_manifest_path="standardized_ligands/partial_packages/p/partial_package_manifest.json",
            ligand_preparation_contract_sha256="6" * 64,
        )
        result = validate_partial35_inventory_contract(
            inventory,
            self.candidates,
            self.base_protocol,
            expected_partial_package_manifest_sha256="5" * 64,
        )
        self.assertEqual(105, result["approved_pdbqt_count"])
        self.assertEqual(15, result["qc_excluded_slot_count"])
        changed = copy.deepcopy(inventory)
        changed["slots"][0]["readiness_status"] = "conversion_audit_only"
        with self.assertRaisesRegex(ValueError, "105"):
            validate_partial35_inventory_contract(
                changed,
                self.candidates,
                self.base_protocol,
                expected_partial_package_manifest_sha256="5" * 64,
            )


class PrepreparedPdbqtMaterializationTests(unittest.TestCase):
    def test_prepare_one_copies_hash_locked_pdbqt_without_meeko_regeneration(self) -> None:
        spec = json.loads(B3_SPEC.read_text(encoding="utf-8"))
        conf = spec["conformers"][0]
        source_pdbqt = B3_ROOT / "ligands" / "B3_conf01.pdbqt"
        conf["ligand_pdbqt_path"] = source_pdbqt.relative_to(DATA_ROOT).as_posix()
        conf["ligand_pdbqt_sha256"] = sha256(source_pdbqt)
        scratch = ROOT / "test_output" / f"partial35_materialize_{uuid.uuid4().hex}"
        scratch.mkdir(parents=True)
        spec_path = scratch / "spec.json"
        spec_path.write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8", newline="\n")
        output = scratch / "attempt"
        with patch(
            "prepare_standardized_vina.prepare_backbone_rigid_ligand",
            side_effect=AssertionError("Meeko regeneration must not run"),
        ):
            audit_path = prepare_one(spec_path, DATA_ROOT, output, "conf01")
        copied = output / "ligands" / "B3_conf01.pdbqt"
        self.assertEqual(sha256(source_pdbqt), sha256(copied))
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        self.assertEqual("copied_hash_locked_preprepared_pdbqt", audit["generated_pdbqt"]["materialization_mode"])
        self.assertEqual(conf["ligand_pdbqt_sha256"], audit["generated_pdbqt"]["source_prepared_pdbqt_sha256"])


class ExploratoryPdbqtPreflightTests(unittest.TestCase):
    @staticmethod
    def promote_a4(inventory: dict, *, wrong_hash: bool = False) -> None:
        for slot in inventory["slots"]:
            if slot["candidate_id"] != "A4":
                continue
            run_id = f"A4_{slot['conformer_name']}"
            relative = f"docking_test/vina_standardized_v04_20260907/A4/ligands/{run_id}.pdbqt"
            actual_hash = sha256(DATA_ROOT / relative)
            slot.update(
                {
                    "readiness_status": "approved_for_exploratory_screening",
                    "ligand_pdbqt_path": relative,
                    "ligand_pdbqt_sha256": "0" * 64 if wrong_hash else actual_hash,
                    "pdbqt_audit_path": slot["audit_path"],
                    "pdbqt_audit_sha256": slot["audit_sha256"],
                    "generation_provenance": {
                        "protocol_id": "tfr1_top40_sdfgen_competition_v1_b3derived",
                        "protocol_sha256": "4" * 64,
                        "partial_package_id": "competition_v1_partial35_mmff_qc_pass",
                        "partial_package_manifest_sha256": "5" * 64,
                    },
                    "approval_basis": "PARTIAL35_SDF_AND_PDBQT_AUDIT_PASS",
                }
            )

    def test_dry_run_accepts_and_hash_audits_preprepared_pdbqt(self) -> None:
        manifest, protocol = write_contract_fixture(
            "partial35_prepared_a4",
            inventory_change=lambda inventory: self.promote_a4(inventory),
        )
        result = run_dry_run(
            manifest,
            protocol,
            DATA_ROOT,
            unused_output("partial35_prepared_a4_out"),
            VINA,
            candidate_ids=["A4"],
        )
        self.assertEqual("DRY_RUN_PASS", result["outcome"])
        conformers = result["plan"][0]["conformers"]
        self.assertTrue(all(row["pdbqt_artifact_validation_status"] == "PASS" for row in conformers))
        self.assertTrue(all("resolved_ligand_pdbqt" in row for row in conformers))
        self.assertEqual(0, result["vina_docking_process_count"])

    def test_wrong_prepared_pdbqt_hash_is_candidate_failure(self) -> None:
        manifest, protocol = write_contract_fixture(
            "partial35_wrong_pdbqt",
            inventory_change=lambda inventory: self.promote_a4(inventory, wrong_hash=True),
        )
        result = run_dry_run(
            manifest,
            protocol,
            DATA_ROOT,
            unused_output("partial35_wrong_pdbqt_out"),
            VINA,
            candidate_ids=["A4"],
        )
        self.assertEqual("DRY_RUN_FAIL_CANDIDATE", result["outcome"])
        self.assertTrue(any("PDBQT" in row["error_code"] for row in result["failures"]))
        self.assertEqual(0, result["vina_docking_process_count"])


if __name__ == "__main__":
    unittest.main()
