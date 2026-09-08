from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
import sys
import unittest
import uuid
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
COMPETITION_ROOT = REPOSITORY_ROOT.parent
ANALYZER = REPOSITORY_ROOT / "scripts" / "analyze_peptide_ensemble.py"
SPEC = REPOSITORY_ROOT / "specs" / "B3_ensemble_v04.json"


class B3V04RealIntegrationTests(unittest.TestCase):
    def test_real_standardized_b3_is_candidate_agnostic_and_score_audited(self) -> None:
        output = REPOSITORY_ROOT / "test_output" / f"B3_v04_{uuid.uuid4().hex}" / "B3"
        completed = subprocess.run(
            [
                sys.executable,
                str(ANALYZER),
                "--spec",
                str(SPEC),
                "--data-root",
                str(COMPETITION_ROOT),
                "--output-dir",
                str(output),
            ],
            cwd=REPOSITORY_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(0, completed.returncode, completed.stderr)
        with (output / "pose_table.tsv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(60, len(rows))
        self.assertEqual(
            {"conf01": 20, "conf02": 20, "conf03": 20},
            {name: sum(row["conformer"] == name for row in rows) for name in ("conf01", "conf02", "conf03")},
        )
        self.assertTrue(all(row["vina_score_crosscheck_status"] == "MATCH" for row in rows))
        self.assertTrue(
            all(
                row["legacy_basic_geometry_pass_source"]
                == "v04_recalculated_legacy_definitions"
                for row in rows
            )
        )
        summary = json.loads((output / "peptide_summary.json").read_text(encoding="utf-8"))
        self.assertEqual("B3", summary["identity"]["candidate"])
        self.assertEqual("IPVGNCRYMQ", summary["identity"]["sequence"])
        self.assertEqual(10, summary["identity"]["peptide_length"])
        self.assertEqual("formal_standardized", summary["screening_protocol_status"])
        self.assertEqual("v04_b79103f6b9650da06a0a", summary["comparison_protocol_id"])
        self.assertEqual(60, summary["vina_score_provenance"]["crosscheck_status_counts"]["MATCH"])
        self.assertEqual(
            "v04_recalculated_legacy_definitions",
            summary["contact_fraction_audit"]["source"],
        )
        with (output / "candidate_summary.tsv").open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            candidate_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(1, len(candidate_rows))
        candidate = candidate_rows[0]
        self.assertTrue(
            {
                "conformer_count",
                "eligible_pose_fraction",
                "recurrent_direct_residue_count",
                "representative_count",
                "analysis_status",
                "source_spec_sha256",
                "receptor_registry_sha256",
                "receptor_monomer_sha256",
                "docking_receptor_pdbqt_sha256",
                "receptor_dimer_sha256",
                "chemistry_contract_sha256",
            }.issubset(candidate)
        )
        self.assertEqual("3", candidate["conformer_count"])
        self.assertEqual("60", candidate["total_pose_count"])
        self.assertEqual("59", candidate["eligible_pose_count"])
        self.assertAlmostEqual(59 / 60, float(candidate["eligible_pose_fraction"]), places=9)
        self.assertEqual("6", candidate["recurrent_direct_residue_count"])
        self.assertEqual("9", candidate["representative_count"])
        self.assertEqual("PASS", candidate["analysis_status"])
        self.assertEqual(
            hashlib.sha256(SPEC.read_bytes()).hexdigest(),
            candidate["source_spec_sha256"],
        )
        for field in (
            "source_spec_sha256",
            "receptor_registry_sha256",
            "receptor_monomer_sha256",
            "docking_receptor_pdbqt_sha256",
            "receptor_dimer_sha256",
            "chemistry_contract_sha256",
        ):
            self.assertRegex(candidate[field], re.compile(r"^[0-9a-f]{64}$"))


if __name__ == "__main__":
    unittest.main()
