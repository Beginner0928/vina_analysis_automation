from __future__ import annotations

import json
import sys
import unittest
import uuid
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
COMPETITION_ROOT = REPOSITORY_ROOT.parent
SITE_PACKAGES = COMPETITION_ROOT / ".venv" / "Lib" / "site-packages"
sys.path.insert(0, str(SITE_PACKAGES))
sys.path.insert(0, str(REPOSITORY_ROOT / "scripts"))

from standardized_vina_inputs import (  # noqa: E402
    audit_ligand_pdbqt,
    audit_source_sdf,
    compare_pdbqt_coordinates,
    prepare_backbone_rigid_ligand,
    read_box_registry,
)


class StandardizedVinaInputTests(unittest.TestCase):
    def test_b3_standardized_sdf_prepares_without_coordinate_regeneration(self) -> None:
        source = (
            COMPETITION_ROOT
            / "GPT_vina_branch_20260906"
            / "peptides"
            / "B3"
            / "B3_conf01.sdf"
        )
        audit = audit_source_sdf(source, "B3", "IPVGNCRYMQ", 1)
        self.assertEqual(163, audit["atom_count"])
        self.assertEqual(81, audit["heavy_atom_count"])
        self.assertEqual(1, audit["formal_charge"])
        self.assertEqual(12, audit["source_conformer_id"])
        self.assertEqual(1, audit["chemistry_checks"]["n_terminus_formal_charge"])
        self.assertEqual(3, audit["chemistry_checks"]["n_terminus_hydrogen_count"])
        self.assertEqual(-1, audit["chemistry_checks"]["c_terminus_oxt_formal_charge"])
        self.assertEqual(1, audit["chemistry_checks"]["arginine_nh1_formal_charge"])
        self.assertEqual(0, audit["chemistry_checks"]["sulfur_sulfur_bond_count"])
        self.assertTrue(audit["chemistry_checks"]["heavy_atom_chirality_matches_L_template"])

        output = (
            REPOSITORY_ROOT
            / "test_output"
            / f"standardized_input_{uuid.uuid4().hex}"
            / "B3_conf01.pdbqt"
        )
        source_hash_before = audit["sha256"]
        prep = prepare_backbone_rigid_ligand(source, output)
        source_hash_after = audit_source_sdf(source, "B3", "IPVGNCRYMQ", 1)["sha256"]
        self.assertEqual(source_hash_before, source_hash_after)
        self.assertEqual(21, prep["torsdof"])

        pdbqt = audit_ligand_pdbqt(output, source, "IPVGNCRYMQ")
        self.assertEqual(81, pdbqt["heavy_atom_count"])
        self.assertEqual(81, pdbqt["mapped_heavy_atom_count"])
        self.assertEqual(list(range(1, 11)), pdbqt["mapped_residue_numbers"])
        coordinate_audit = compare_pdbqt_coordinates(source, output)
        self.assertEqual(103, coordinate_audit["mapped_atom_count"])
        self.assertLessEqual(coordinate_audit["maximum_coordinate_delta_A"], 0.001)
        self.assertEqual("MATCH_WITHIN_PDBQT_PRECISION", coordinate_audit["status"])

    def test_group_b_box_is_read_from_registry(self) -> None:
        box = read_box_registry(
            COMPETITION_ROOT
            / "GPT_vina_branch_20260906"
            / "config"
            / "vina_boxes.tsv",
            "h26d3_boltz_msa_s2718",
        )
        self.assertEqual((-15.739, 16.133, 16.584), box["center_A"])
        self.assertEqual((38.023, 26.643, 27.375), box["size_A"])


if __name__ == "__main__":
    unittest.main()
