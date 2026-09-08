from __future__ import annotations

import json
import sys
import unittest
import uuid
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
COMPETITION_ROOT = REPOSITORY_ROOT.parent
sys.path.insert(0, str(REPOSITORY_ROOT / "scripts"))

from receptor_registry import (  # noqa: E402
    resolve_receptor_bundle,
    validate_spec_against_bundle,
)


REGISTRY = REPOSITORY_ROOT / "specs" / "receptor_registry_v04.json"
EXPECTED_TEMPLATES = {
    "A": "h26d3_af3_s2718_m2",
    "B": "h26d3_boltz_msa_s2718",
    "C": "h26d3_rank009_rb_s1701_dom",
    "D": "h26d3_rank016_mdref_s1701_r1",
}
EXPECTED_BOXES = {
    "A": ((-11.144, 23.985, -11.463), (28.44, 36.605, 33.415)),
    "B": ((-15.739, 16.133, 16.584), (38.023, 26.643, 27.375)),
    "C": ((-3.026, 24.096, -5.173), (30.501, 37.743, 27.644)),
    "D": ((28.128, 3.3, -9.748), (40.255, 29.623, 30.63)),
}
EXPECTED_C_D_PROVENANCE = {
    "C": {
        "normalized_receptor_monomer_pdb_sha256": "ff73d6851c20a4d47c70247a100a5c019372b28f4223afbdd4fed6d6604a99b6",
        "docking_receptor_pdbqt_sha256": "efa6f988bbc2384526137b9a963c195b557322529045d84b460dcf11a519af95",
        "receptor_dimer_pdb_sha256": "a51400e0843812a5bfe8fe1fc8cf43079f293e46f697326ec0101a51ba857935",
    },
    "D": {
        "normalized_receptor_monomer_pdb_sha256": "329a16214e0c0fafafeadc2e7fd9911ebf681f47d171edc0304c9c8e166ada53",
        "docking_receptor_pdbqt_sha256": "561ef56c9691cd339593b2ffda6439d551e6b245e93a0271470ed940b04f38d2",
        "receptor_dimer_pdb_sha256": "45b6e859bef755ea9c2a7898f19f2ebd8cc939ac4aaa62d9ca7ca5221acc9fec",
    },
}


class ReceptorRegistryV04Tests(unittest.TestCase):
    def test_all_groups_resolve_to_locked_template_bundles(self) -> None:
        for group, template in EXPECTED_TEMPLATES.items():
            with self.subTest(group=group):
                bundle = resolve_receptor_bundle(REGISTRY, COMPETITION_ROOT, group)
                self.assertEqual(template, bundle["template_id"])
                self.assertEqual(group, bundle["group"])
                self.assertEqual("PASS", bundle["preflight_status"])
                self.assertEqual(
                    [150, 151, 154, 158, 159, 160, 161, 163, 385],
                    bundle["target_residues"],
                )

    def test_unknown_group_fails_without_fallback(self) -> None:
        with self.assertRaisesRegex(ValueError, "Unknown receptor group"):
            resolve_receptor_bundle(REGISTRY, COMPETITION_ROOT, "Z")

    def test_c_and_d_use_normalized_production_monomers(self) -> None:
        c = resolve_receptor_bundle(REGISTRY, COMPETITION_ROOT, "C")
        d = resolve_receptor_bundle(REGISTRY, COMPETITION_ROOT, "D")
        self.assertTrue(c["receptor_monomer_pdb"].name.endswith("_normalized.pdb"))
        self.assertTrue(d["receptor_monomer_pdb"].name.endswith("_normalized.pdb"))
        self.assertIn("rank009", c["receptor_monomer_pdb"].name)
        self.assertIn("rank016", d["receptor_monomer_pdb"].name)

    def test_c_and_d_require_explicit_approved_normalization_provenance(self) -> None:
        registry = self.changed_registry()
        for group, expected_hashes in EXPECTED_C_D_PROVENANCE.items():
            with self.subTest(group=group):
                entry = registry["groups"][group]
                self.assertTrue(entry["normalized_receptor_required"])
                self.assertEqual("PASS", entry["normalization_status"])
                self.assertEqual("approved", entry["production_status"])
                self.assertEqual("RESOLVED", entry["receptor_preparation_blocker"])
                self.assertEqual(expected_hashes, entry["production_provenance_sha256"])

                bundle = resolve_receptor_bundle(REGISTRY, COMPETITION_ROOT, group)
                self.assertEqual("PASS", bundle["normalization_status"])
                self.assertEqual("approved", bundle["production_status"])
                self.assertEqual("RESOLVED", bundle["receptor_preparation_blocker"])
                self.assertEqual(expected_hashes, bundle["production_provenance_sha256"])

    def test_c_and_d_explicit_provenance_hash_mismatch_fails(self) -> None:
        changed = self.changed_registry()
        changed["groups"]["C"]["production_provenance_sha256"][
            "docking_receptor_pdbqt_sha256"
        ] = "0" * 64
        path = self.write_registry(changed, "explicit_provenance_hash")
        with self.assertRaisesRegex(ValueError, "production provenance SHA-256 mismatch"):
            resolve_receptor_bundle(path, COMPETITION_ROOT, "C")

    def test_serorder_diagnostic_is_rejected_as_production(self) -> None:
        changed = self.changed_registry()
        changed["groups"]["C"]["receptor_monomer_pdb"] = (
            "GPT_vina_branch_20260906/receptor_pdb_fixed/"
            "h26d3_rank009_rb_s1701_dom_TfR1_A_SERorder_test.pdb"
        )
        path = self.write_registry(changed, "serorder")
        with self.assertRaisesRegex(ValueError, "SERorder_test"):
            resolve_receptor_bundle(path, COMPETITION_ROOT, "C")

    def test_hash_mismatch_fails(self) -> None:
        changed = self.changed_registry()
        changed["groups"]["A"]["sha256"]["receptor_pdbqt"] = "0" * 64
        path = self.write_registry(changed, "hash")
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            resolve_receptor_bundle(path, COMPETITION_ROOT, "A")

    def test_missing_monomer_pdbqt_or_dimer_fails(self) -> None:
        for field in ("receptor_monomer_pdb", "receptor_pdbqt", "receptor_dimer_pdb"):
            with self.subTest(field=field):
                changed = self.changed_registry()
                template = changed["groups"]["B"]["template_id"]
                changed["groups"]["B"][field] = f"missing/{template}_{field}.pdb"
                path = self.write_registry(changed, field)
                with self.assertRaisesRegex(FileNotFoundError, "missing"):
                    resolve_receptor_bundle(path, COMPETITION_ROOT, "B")

    def test_cross_template_component_is_rejected(self) -> None:
        changed = self.changed_registry()
        changed["groups"]["B"]["receptor_pdbqt"] = changed["groups"]["A"]["receptor_pdbqt"]
        changed["groups"]["B"]["sha256"]["receptor_pdbqt"] = changed["groups"]["A"]["sha256"]["receptor_pdbqt"]
        path = self.write_registry(changed, "cross_template")
        with self.assertRaisesRegex(ValueError, "cross-template"):
            resolve_receptor_bundle(path, COMPETITION_ROOT, "B")

    def test_box_values_are_group_specific_and_match_audited_source(self) -> None:
        actual = {}
        for group in EXPECTED_TEMPLATES:
            bundle = resolve_receptor_bundle(REGISTRY, COMPETITION_ROOT, group)
            actual[group] = (tuple(bundle["vina_box"]["center"]), tuple(bundle["vina_box"]["size"]))
        self.assertEqual(EXPECTED_BOXES, actual)
        self.assertEqual(4, len(set(actual.values())))

    def test_existing_b3_protocol_matches_group_b_registry(self) -> None:
        bundle = resolve_receptor_bundle(REGISTRY, COMPETITION_ROOT, "B")
        protocol = json.loads(
            (
                COMPETITION_ROOT
                / "docking_test/vina_standardized_v04_20260907/B3/protocol.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(bundle["template_id"], protocol["template_id"])
        self.assertEqual(bundle["sha256"]["receptor_pdbqt"], protocol["receptor_pdbqt_sha256"])
        self.assertEqual(tuple(bundle["vina_box"]["center"]), tuple(protocol["box_center_A"]))
        self.assertEqual(tuple(bundle["vina_box"]["size"]), tuple(protocol["box_size_A"]))

    def test_cross_template_analysis_spec_is_rejected(self) -> None:
        bundle = resolve_receptor_bundle(REGISTRY, COMPETITION_ROOT, "B")
        spec = {
            "group": "B",
            "template_id": EXPECTED_TEMPLATES["A"],
            "receptor_id": f"{EXPECTED_TEMPLATES['A']}_TfR1_A",
            "target_residues": bundle["target_residues"],
            "docking_protocol": {
                "box_center_A": bundle["vina_box"]["center"],
                "box_size_A": bundle["vina_box"]["size"],
            },
        }
        with self.assertRaisesRegex(ValueError, "template_id.*registry"):
            validate_spec_against_bundle(spec, bundle)

    def test_analysis_input_hashes_must_match_registry_bundle(self) -> None:
        bundle = resolve_receptor_bundle(REGISTRY, COMPETITION_ROOT, "B")
        spec = {
            "group": "B",
            "template_id": bundle["template_id"],
            "receptor_id": f"{bundle['template_id']}_TfR1_A",
            "target_residues": bundle["target_residues"],
            "docking_protocol": {
                "box_center_A": bundle["vina_box"]["center"],
                "box_size_A": bundle["vina_box"]["size"],
            },
            "input_sha256": {
                "receptor_monomer_path": "0" * 64,
                "receptor_dimer_path": bundle["sha256"]["receptor_dimer_pdb"],
                "docking_receptor_pdbqt_path": bundle["sha256"]["receptor_pdbqt"],
            },
        }
        with self.assertRaisesRegex(ValueError, "receptor_monomer_path.*registry"):
            validate_spec_against_bundle(spec, bundle)

    def changed_registry(self) -> dict[str, object]:
        return json.loads(REGISTRY.read_text(encoding="utf-8"))

    def write_registry(self, value: dict[str, object], label: str) -> Path:
        path = (
            REPOSITORY_ROOT
            / "test_output"
            / f"registry_{label}_{uuid.uuid4().hex}.json"
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
        return path


if __name__ == "__main__":
    unittest.main()
