from __future__ import annotations

import sys
import unittest
import json
from pathlib import Path

from rdkit import Chem


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from competition_conformer_core import (  # noqa: E402
    PROTOCOL_ID,
    audit_standardized_peptide,
    bond_graph_hash,
    build_etkdg_parameters,
    build_standardized_peptide,
    derive_candidate_seed,
    expected_formal_charge,
    load_generation_protocol,
    select_formal_conformers,
    stereochemistry_hash,
)


PROTOCOL_PATH = ROOT / "specs/competition_conformer_protocol_v1.json"


def atom_by_site(molecule, residue_number: int, atom_name: str):
    matches = []
    for atom in molecule.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info and info.GetResidueNumber() == residue_number and info.GetName().strip() == atom_name:
            matches.append(atom)
    if len(matches) != 1:
        raise AssertionError(f"Expected one atom at residue {residue_number} {atom_name}; got {len(matches)}")
    return matches[0]


def explicit_hydrogen_neighbors(atom) -> int:
    return sum(neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors())


class CompetitionConformerChemistryTests(unittest.TestCase):
    def test_seed_literals_match_frozen_algorithm(self) -> None:
        self.assertEqual(1420034032, derive_candidate_seed(PROTOCOL_ID, "A4"))
        self.assertEqual(524445188, derive_candidate_seed(PROTOCOL_ID, "B3"))
        self.assertEqual(2119376187, derive_candidate_seed(PROTOCOL_ID, "C2"))
        self.assertEqual(891847413, derive_candidate_seed(PROTOCOL_ID, "D8"))

    def test_expected_formal_charge_uses_frozen_residue_contract(self) -> None:
        self.assertEqual(-1, expected_formal_charge("SPSNFITMYDW"))
        self.assertEqual(1, expected_formal_charge("IPVGNCRYMQ"))
        self.assertEqual(1, expected_formal_charge("RDICQFQFHK"))
        self.assertEqual(1, expected_formal_charge("ITPQVPWIRY"))

    def test_builds_frozen_termini_and_sidechain_protonation(self) -> None:
        sequence = "HCKDERY"
        molecule = build_standardized_peptide(sequence)

        self.assertEqual(expected_formal_charge(sequence), Chem.GetFormalCharge(molecule))
        n_terminal = atom_by_site(molecule, 1, "N")
        c_terminal = atom_by_site(molecule, len(sequence), "OXT")
        self.assertEqual((1, 3), (n_terminal.GetFormalCharge(), explicit_hydrogen_neighbors(n_terminal)))
        self.assertEqual((-1, 0), (c_terminal.GetFormalCharge(), explicit_hydrogen_neighbors(c_terminal)))
        self.assertEqual(1, atom_by_site(molecule, 3, "NZ").GetFormalCharge())
        self.assertEqual(-1, atom_by_site(molecule, 4, "OD2").GetFormalCharge())
        self.assertEqual(-1, atom_by_site(molecule, 5, "OE2").GetFormalCharge())
        self.assertEqual(1, atom_by_site(molecule, 6, "NH1").GetFormalCharge())

    def test_histidine_is_neutral_hie_like(self) -> None:
        molecule = build_standardized_peptide("HA")
        nd1 = atom_by_site(molecule, 1, "ND1")
        ne2 = atom_by_site(molecule, 1, "NE2")

        self.assertEqual((0, 0), (nd1.GetFormalCharge(), explicit_hydrogen_neighbors(nd1)))
        self.assertEqual((0, 1), (ne2.GetFormalCharge(), explicit_hydrogen_neighbors(ne2)))

    def test_cysteine_and_tyrosine_remain_neutral_protonated_sidechains(self) -> None:
        molecule = build_standardized_peptide("CY")
        sulfur = atom_by_site(molecule, 1, "SG")
        tyrosine_oxygen = atom_by_site(molecule, 2, "OH")

        self.assertEqual((0, 1), (sulfur.GetFormalCharge(), explicit_hydrogen_neighbors(sulfur)))
        self.assertEqual((0, 1), (tyrosine_oxygen.GetFormalCharge(), explicit_hydrogen_neighbors(tyrosine_oxygen)))
        self.assertFalse(any(
            bond.GetBeginAtom().GetAtomicNum() == 16 and bond.GetEndAtom().GetAtomicNum() == 16
            for bond in molecule.GetBonds()
        ))

    def test_audit_reports_identity_graph_and_standard_l_stereochemistry(self) -> None:
        row = {"candidate_id": "C2", "group": "C", "sequence": "RDICQFQFHK", "peptide_length": 10}
        molecule = build_standardized_peptide(row["sequence"])

        audit = audit_standardized_peptide(molecule, row)

        self.assertEqual("PASS", audit["status"])
        self.assertEqual(row["sequence"], audit["sequence"])
        self.assertEqual(row["peptide_length"], audit["peptide_length"])
        self.assertEqual(1, audit["formal_charge"])
        self.assertTrue(audit["chemistry_checks"]["standard_L_stereochemistry"])
        self.assertEqual(0, audit["chemistry_checks"]["sulfur_sulfur_bond_count"])
        self.assertEqual(64, len(bond_graph_hash(molecule)))
        self.assertEqual(64, len(stereochemistry_hash(molecule)))

    def test_audit_rejects_sequence_identity_mismatch(self) -> None:
        molecule = build_standardized_peptide("AA")
        wrong = {"candidate_id": "A1", "group": "A", "sequence": "AG", "peptide_length": 2}

        with self.assertRaisesRegex(ValueError, "heavy-atom identity/order"):
            audit_standardized_peptide(molecule, wrong)

    def test_audit_rejects_unintended_disulfide(self) -> None:
        molecule = build_standardized_peptide("CC")
        sulfurs = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 16]
        rw = Chem.RWMol(molecule)
        rw.AddBond(sulfurs[0], sulfurs[1], Chem.BondType.SINGLE)

        with self.assertRaisesRegex(ValueError, "S-S"):
            audit_standardized_peptide(
                rw.GetMol(),
                {"candidate_id": "X1", "group": "X", "sequence": "CC", "peptide_length": 2},
            )

    def test_all_canonical_top40_sequences_pass_chemistry_construction_audit(self) -> None:
        manifest = json.loads((ROOT / "specs/top40_candidates_v05.json").read_text(encoding="utf-8"))

        audits = [
            audit_standardized_peptide(build_standardized_peptide(row["sequence"]), row)
            for row in manifest["candidates"]
        ]

        self.assertEqual(40, len(audits))
        self.assertEqual({"PASS"}, {audit["status"] for audit in audits})


class CompetitionConformerGenerationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.protocol = load_generation_protocol(PROTOCOL_PATH)

    def test_etkdg_parameters_are_frozen_and_serial(self) -> None:
        params = build_etkdg_parameters(self.protocol, 1420034032)

        self.assertEqual(1420034032, params.randomSeed)
        self.assertTrue(params.useRandomCoords)
        self.assertEqual(0.25, params.pruneRmsThresh)
        self.assertTrue(params.enforceChirality)
        self.assertEqual(1, params.numThreads)
        self.assertEqual(2, params.ETversion)

    def test_selection_uses_energy_then_source_id_and_preserves_ids(self) -> None:
        rows = [
            {"source_conformer_id": 7, "convergence_status": 0, "energy": -10.0},
            {"source_conformer_id": 3, "convergence_status": 0, "energy": -10.0},
            {"source_conformer_id": 8, "convergence_status": 0, "energy": -9.0},
            {"source_conformer_id": 9, "convergence_status": 0, "energy": -8.0},
        ]
        distances = {
            frozenset((3, 8)): 2.1,
            frozenset((3, 9)): 2.2,
            frozenset((8, 9)): 2.3,
            frozenset((7, 3)): 0.4,
            frozenset((7, 8)): 2.4,
            frozenset((7, 9)): 2.5,
        }

        result = select_formal_conformers(
            rows,
            lambda left, right: distances[frozenset((left, right))],
            self.protocol,
        )

        self.assertEqual([3, 8, 9], result["selected_source_conformer_ids"])
        self.assertEqual([-10.0, -9.0, -8.0], result["selected_energies"])
        self.assertEqual(2.0, result["rmsd_threshold_used_A"])
        self.assertEqual(
            {"conf01-conf02": 2.1, "conf01-conf03": 2.2, "conf02-conf03": 2.3},
            result["pairwise_rmsd_A"],
        )

    def test_selection_uses_highest_threshold_that_builds_full_set(self) -> None:
        rows = [
            {"source_conformer_id": 1, "convergence_status": 0, "energy": 0.0},
            {"source_conformer_id": 2, "convergence_status": 0, "energy": 1.0},
            {"source_conformer_id": 3, "convergence_status": 0, "energy": 2.0},
        ]
        distances = {
            frozenset((1, 2)): 1.8,
            frozenset((1, 3)): 1.7,
            frozenset((2, 3)): 1.6,
        }

        result = select_formal_conformers(
            rows,
            lambda left, right: distances[frozenset((left, right))],
            self.protocol,
        )

        self.assertEqual(1.5, result["rmsd_threshold_used_A"])

    def test_unconverged_conformers_are_excluded(self) -> None:
        rows = [
            {"source_conformer_id": 1, "convergence_status": 0, "energy": 0.0},
            {"source_conformer_id": 2, "convergence_status": 1, "energy": 0.1},
            {"source_conformer_id": 3, "convergence_status": 0, "energy": 1.0},
            {"source_conformer_id": 4, "convergence_status": 0, "energy": 2.0},
        ]
        distances = {
            frozenset((1, 3)): 2.1,
            frozenset((1, 4)): 2.2,
            frozenset((3, 4)): 2.3,
        }

        result = select_formal_conformers(
            rows,
            lambda left, right: distances[frozenset((left, right))],
            self.protocol,
        )

        self.assertEqual([1, 3, 4], result["selected_source_conformer_ids"])
        self.assertNotIn(2, result["selected_source_conformer_ids"])

    def test_fewer_than_three_converged_conformers_fails_closed(self) -> None:
        rows = [
            {"source_conformer_id": 1, "convergence_status": 0, "energy": 0.0},
            {"source_conformer_id": 2, "convergence_status": 1, "energy": 1.0},
            {"source_conformer_id": 3, "convergence_status": 0, "energy": 2.0},
        ]

        with self.assertRaisesRegex(RuntimeError, "fewer than 3 converged"):
            select_formal_conformers(rows, lambda _left, _right: 9.0, self.protocol)

    def test_no_three_member_set_at_point_five_fails_closed(self) -> None:
        rows = [
            {"source_conformer_id": 1, "convergence_status": 0, "energy": 0.0},
            {"source_conformer_id": 2, "convergence_status": 0, "energy": 1.0},
            {"source_conformer_id": 3, "convergence_status": 0, "energy": 2.0},
        ]

        with self.assertRaisesRegex(RuntimeError, "0.5 A"):
            select_formal_conformers(rows, lambda _left, _right: 0.49, self.protocol)


if __name__ == "__main__":
    unittest.main()
