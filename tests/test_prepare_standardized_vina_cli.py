from __future__ import annotations

import json
import os
import subprocess
import sys
import unittest
import uuid
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
COMPETITION_ROOT = REPOSITORY_ROOT.parent
SITE_PACKAGES = COMPETITION_ROOT / ".venv" / "Lib" / "site-packages"
SCRIPT = REPOSITORY_ROOT / "scripts" / "prepare_standardized_vina.py"
AUDIT_SCRIPT = REPOSITORY_ROOT / "scripts" / "audit_prepared_vina_input.py"
SPEC = REPOSITORY_ROOT / "specs" / "B3_standardized_docking_v04.json"


class PrepareStandardizedVinaCliTests(unittest.TestCase):
    def test_prepares_one_conformer_and_complete_read_only_provenance(self) -> None:
        output_root = REPOSITORY_ROOT / "test_output" / f"b3_prep_{uuid.uuid4().hex}"
        environment = dict(os.environ)
        environment["PYTHONPATH"] = str(SITE_PACKAGES)
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--spec",
                str(SPEC),
                "--data-root",
                str(COMPETITION_ROOT),
                "--output-root",
                str(output_root),
                "--conformer",
                "conf01",
            ],
            cwd=REPOSITORY_ROOT,
            env=environment,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(0, completed.returncode, completed.stderr)
        expected = {
            "source_sdf/B3_conf01.sdf",
            "ligands/B3_conf01.pdbqt",
            "configs/B3_conf01.txt",
            "receptors/h26d3_boltz_msa_s2718_TfR1_A.pdb",
            "receptors/h26d3_boltz_msa_s2718_TfR1_A.pdbqt",
            "receptors/h26d3_boltz_msa_s2718_TfR1_AB.pdb",
            "audit/B3_conf01_input_audit.json",
            "protocol.json",
        }
        actual = {
            path.relative_to(output_root).as_posix()
            for path in output_root.rglob("*")
            if path.is_file()
        }
        self.assertEqual(expected, actual)

        audit = json.loads(
            (output_root / "audit" / "B3_conf01_input_audit.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual("B3_conf01", audit["run_id"])
        self.assertEqual(1, audit["source_sdf"]["formal_charge"])
        self.assertEqual(21, audit["generated_pdbqt"]["torsdof"])
        self.assertFalse(audit["source_sdf"]["coordinates_regenerated"])
        self.assertEqual("MATCH", audit["source_sdf_hash_verification"])
        self.assertEqual("formal_standardized", audit["screening_protocol_status"])
        self.assertEqual("B", audit["receptor_registry"]["group"])
        self.assertEqual("PASS", audit["receptor_registry"]["preflight_status"])

        config = (output_root / "configs" / "B3_conf01.txt").read_text(encoding="utf-8")
        self.assertIn("receptor = receptors/h26d3_boltz_msa_s2718_TfR1_A.pdbqt", config)
        self.assertIn("ligand = ligands/B3_conf01.pdbqt", config)
        self.assertIn("center_x = -15.739", config)
        self.assertIn("cpu = 8", config)

        protocol = json.loads((output_root / "protocol.json").read_text(encoding="utf-8"))
        self.assertEqual("v04_b79103f6b9650da06a0a", protocol["comparison_protocol_id"])
        self.assertEqual("tfr1_receptor_registry_v04", protocol["receptor_registry_id"])

        gate = subprocess.run(
            [
                sys.executable,
                str(AUDIT_SCRIPT),
                "--spec",
                str(SPEC),
                "--data-root",
                str(COMPETITION_ROOT),
                "--root",
                str(output_root),
                "--conformer",
                "conf01",
            ],
            cwd=REPOSITORY_ROOT,
            env=environment,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(0, gate.returncode, gate.stderr)
        gate_payload = json.loads(
            (output_root / "audit/B3_conf01_pre_docking_gate.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual("PASS", gate_payload["gate_status"])
        self.assertEqual("PASS", gate_payload["receptor_registry"]["preflight_status"])


if __name__ == "__main__":
    unittest.main()
