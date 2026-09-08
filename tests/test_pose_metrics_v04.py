from __future__ import annotations

import sys
import unittest
import uuid
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
COMPETITION_ROOT = REPOSITORY_ROOT.parent
SITE_PACKAGES = COMPETITION_ROOT / ".venv" / "Lib" / "site-packages"
sys.path.insert(0, str(SITE_PACKAGES))
sys.path.insert(0, str(REPOSITORY_ROOT / "scripts"))

from pose_analysis_core import CoordinateAtom  # noqa: E402
from pose_metrics_v04 import (  # noqa: E402
    compute_legacy_pose_metrics,
    ligand_serial_to_residue,
)
from standardized_vina_inputs import prepare_backbone_rigid_ligand  # noqa: E402


class PoseMetricsV04Tests(unittest.TestCase):
    def test_standardized_b3_mapping_uses_sdf_identity_not_neutral_hydrogen_count(self) -> None:
        source = COMPETITION_ROOT / "GPT_vina_branch_20260906/peptides/B3/B3_conf01.sdf"
        pdbqt = (
            REPOSITORY_ROOT
            / "test_output"
            / f"pose_metrics_mapping_{uuid.uuid4().hex}"
            / "B3_conf01.pdbqt"
        )
        prepare_backbone_rigid_ligand(source, pdbqt)
        mapping = ligand_serial_to_residue(source, pdbqt, "IPVGNCRYMQ")
        self.assertEqual(81, len(mapping))
        self.assertEqual(list(range(1, 11)), sorted(set(mapping.values())))

    def test_legacy_metric_semantics_remain_explicit(self) -> None:
        ligand = [
            CoordinateAtom(1, "C1", "C", (0.0, 0.0, 0.0), None, "", ""),
            CoordinateAtom(2, "C2", "C", (10.0, 0.0, 0.0), None, "", ""),
        ]
        receptor = [
            CoordinateAtom(101, "CA", "C", (3.0, 0.0, 0.0), 150, "ASN", "A"),
            CoordinateAtom(102, "CA", "C", (14.5, 0.0, 0.0), 151, "ALA", "A"),
        ]
        chain_b = [
            CoordinateAtom(201, "CA", "C", (30.0, 0.0, 0.0), 1, "GLY", "B")
        ]
        metrics = compute_legacy_pose_metrics(
            model=7,
            score=-5.25,
            ligand_atoms=ligand,
            receptor_atoms=receptor,
            chain_b_atoms=chain_b,
            ligand_residue_by_serial={1: 1, 2: 2},
            peptide_length=2,
            target_residues={150, 151},
        )
        self.assertEqual(1.0, metrics["contact_fraction"])
        self.assertEqual("1;2", metrics["contacted_peptide_residues"])
        self.assertEqual("150;151", metrics["target_contacts"])
        self.assertEqual(0, metrics["receptor_clash_pairs"])
        self.assertEqual(1, metrics["basic_geometry_pass"])
        self.assertEqual("legacy_heavy_atoms_distance_le_5A", metrics["contact_definition"])


if __name__ == "__main__":
    unittest.main()
