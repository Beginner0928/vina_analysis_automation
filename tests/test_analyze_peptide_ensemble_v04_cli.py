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
TEST_OUTPUT_ROOT = REPOSITORY_ROOT / "test_output"
ANALYZER = REPOSITORY_ROOT / "scripts" / "analyze_peptide_ensemble.py"
SPEC = REPOSITORY_ROOT / "specs" / "A4_ensemble_v04.json"
SCRIPTS_ROOT = REPOSITORY_ROOT / "scripts"
if str(SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_ROOT))

from analyze_peptide_ensemble import (  # noqa: E402
    CANDIDATE_COLUMNS,
    build_candidate_summary_row,
)


CANDIDATE_SUMMARY_REQUIRED_COLUMNS = {
    "summary_schema_version",
    "analysis_id",
    "candidate",
    "group",
    "sequence",
    "peptide_length",
    "template_id",
    "receptor_id",
    "screening_protocol_status",
    "comparison_protocol_id",
    "conformer_count",
    "total_pose_count",
    "eligible_pose_count",
    "eligible_pose_fraction",
    "best_raw_vina_score",
    "best_raw_pose_id",
    "median_raw_vina_score",
    "best_eligible_vina_score",
    "best_eligible_pose_id",
    "median_eligible_vina_score",
    "mean_eligible_vina_score",
    "q25_eligible_vina_score",
    "q75_eligible_vina_score",
    "mean_conformer_median_eligible_vina_score",
    "best_eligible_vina_score_range_across_conformers",
    "shared_recurrent_direct_core",
    "recurrent_direct_residue_count",
    "representative_count",
    "warning_count",
    "analysis_status",
    "source_spec_sha256",
    "receptor_registry_sha256",
    "receptor_monomer_sha256",
    "docking_receptor_pdbqt_sha256",
    "receptor_dimer_sha256",
    "chemistry_contract_sha256",
}
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
REQUIRED_OUTPUTS = {
    "pose_table.tsv",
    "pose_target_contacts.tsv",
    "target_contact_frequency.tsv",
    "conformer_summary.tsv",
    "representative_poses.tsv",
    "vina_score_summary.tsv",
    "candidate_summary.tsv",
    "peptide_summary.md",
    "peptide_summary.json",
    "run_manifest.txt",
}


class A4V04IntegrationTests(unittest.TestCase):
    def test_v04_outputs_primary_score_and_protocol_provenance(self) -> None:
        output_dir = TEST_OUTPUT_ROOT / f"A4_v04_{uuid.uuid4().hex}" / "A4"
        completed = subprocess.run(
            [
                sys.executable,
                str(ANALYZER),
                "--spec",
                str(SPEC),
                "--data-root",
                str(COMPETITION_ROOT),
                "--output-dir",
                str(output_dir),
            ],
            cwd=REPOSITORY_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(0, completed.returncode, completed.stderr)
        self.assertEqual(REQUIRED_OUTPUTS, {path.name for path in output_dir.iterdir()})

        with (output_dir / "pose_table.tsv").open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            poses = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(59, len(poses))
        self.assertTrue(
            all(row["vina_score_primary_source"] == "docking_output_pdbqt_remark" for row in poses)
        )
        self.assertTrue(all(row["vina_score_crosscheck_status"] == "MATCH" for row in poses))

        with (output_dir / "vina_score_summary.tsv").open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            score_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(["ALL", "conf01", "conf02", "conf03"], [row["conformer"] for row in score_rows])
        overall = score_rows[0]
        self.assertEqual("A4_conf02_m1", overall["best_raw_pose_id"])
        self.assertAlmostEqual(-6.470, float(overall["best_raw_vina_score"]), places=6)
        self.assertEqual("A4_conf02_m1", overall["best_eligible_pose_id"])

        with (output_dir / "candidate_summary.tsv").open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            candidate_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(1, len(candidate_rows))
        candidate = candidate_rows[0]
        self.assertTrue(CANDIDATE_SUMMARY_REQUIRED_COLUMNS.issubset(candidate))
        self.assertEqual("legacy_regression_only", candidate["screening_protocol_status"])
        self.assertEqual("docking_output_pdbqt_remark", candidate["vina_score_primary_source"])
        self.assertEqual("0.4", candidate["summary_schema_version"])
        self.assertEqual("3", candidate["conformer_count"])
        self.assertEqual("59", candidate["total_pose_count"])
        self.assertAlmostEqual(
            int(candidate["eligible_pose_count"]) / int(candidate["total_pose_count"]),
            float(candidate["eligible_pose_fraction"]),
            places=9,
        )
        self.assertEqual(
            len(set(candidate["shared_recurrent_direct_core"].split(";"))),
            int(candidate["recurrent_direct_residue_count"]),
        )
        self.assertEqual(
            len(candidate["representative_pose_ids"].split(";")),
            int(candidate["representative_count"]),
        )
        self.assertEqual("0", candidate["warning_count"])
        self.assertEqual("PASS", candidate["analysis_status"])
        for field in (
            "source_spec_sha256",
            "receptor_registry_sha256",
            "receptor_monomer_sha256",
            "docking_receptor_pdbqt_sha256",
            "receptor_dimer_sha256",
            "chemistry_contract_sha256",
        ):
            self.assertRegex(candidate[field], SHA256_RE)
        self.assertEqual(
            hashlib.sha256(SPEC.read_bytes()).hexdigest(),
            candidate["source_spec_sha256"],
        )

        summary = json.loads((output_dir / "peptide_summary.json").read_text(encoding="utf-8"))
        self.assertEqual("0.4", summary["schema_version"])
        self.assertEqual("A4", summary["identity"]["candidate"])
        self.assertEqual(11, summary["identity"]["peptide_length"])
        self.assertEqual("legacy_regression_only", summary["screening_protocol_status"])
        self.assertEqual("docking_output_pdbqt_remark", summary["vina_score_provenance"]["primary_source"])
        self.assertEqual("PASS", summary["receptor_registry"]["preflight_status"])
        self.assertEqual("A", summary["receptor_registry"]["group"])

        markdown = (output_dir / "peptide_summary.md").read_text(encoding="utf-8")
        self.assertLess(
            markdown.index("## Vina score provenance"),
            markdown.index("## Legacy 5 A receptor contact-fraction provenance"),
        )
        self.assertIn(
            "## Legacy 5 A receptor contact-fraction provenance\n\n"
            "fraction of peptide residues",
            markdown,
        )

        manifest = (output_dir / "run_manifest.txt").read_text(encoding="utf-8")
        self.assertIn("vina_score_primary_source=docking_output_pdbqt_remark", manifest)
        self.assertIn("screening_protocol_status=legacy_regression_only", manifest)
        self.assertIn("score_crosscheck.MATCH=59", manifest)
        self.assertIn("receptor_registry.preflight_status=PASS", manifest)

    def test_score_mismatch_fails_without_analysis_package(self) -> None:
        run_id = uuid.uuid4().hex
        spec = json.loads(SPEC.read_text(encoding="utf-8"))
        source_metrics = (
            COMPETITION_ROOT
            / spec["conformers"][0]["pose_metrics_path"]
        )
        changed_metrics = TEST_OUTPUT_ROOT / f"A4_v04_mismatch_{run_id}_metrics.tsv"
        changed_metrics.parent.mkdir(parents=True, exist_ok=True)
        lines = source_metrics.read_text(encoding="utf-8").splitlines()
        fields = lines[1].split("\t")
        fields[1] = str(float(fields[1]) + 0.01)
        lines[1] = "\t".join(fields)
        changed_metrics.write_text("\n".join(lines) + "\n", encoding="utf-8")
        spec["conformers"][0]["pose_metrics_path"] = changed_metrics.relative_to(
            COMPETITION_ROOT
        ).as_posix()
        del spec["input_sha256"]["conformers"]["conf01"]["pose_metrics_path"]
        changed_spec = TEST_OUTPUT_ROOT / f"A4_v04_mismatch_{run_id}.json"
        changed_spec.write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8")
        output_dir = TEST_OUTPUT_ROOT / f"A4_v04_mismatch_{run_id}" / "A4"

        completed = subprocess.run(
            [
                sys.executable,
                str(ANALYZER),
                "--spec",
                str(changed_spec),
                "--data-root",
                str(COMPETITION_ROOT),
                "--output-dir",
                str(output_dir),
            ],
            cwd=REPOSITORY_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertNotEqual(0, completed.returncode)
        self.assertIn("WARNING_SCORE_MISMATCH", completed.stderr)
        self.assertFalse(output_dir.exists())


class CandidateSummaryContractTests(unittest.TestCase):
    @staticmethod
    def _inputs() -> dict[str, object]:
        detail = lambda value: {"sha256": value * 64}
        return {
            "source_spec_path": detail("1"),
            "receptor_registry_path": detail("2"),
            "receptor_monomer_path": detail("3"),
            "docking_receptor_pdbqt_path": detail("4"),
            "receptor_dimer_path": detail("5"),
            "chemistry_contract_path": detail("6"),
        }

    @staticmethod
    def _score_row(total: int, eligible: int) -> dict[str, object]:
        return {
            "scope": "peptide",
            "conformer": "ALL",
            "raw_pose_count": total,
            "eligible_pose_count": eligible,
            "eligible_pose_fraction": eligible / total,
            "best_raw_vina_score": -6.0,
            "best_raw_pose_id": "X_conf01_m1",
            "median_raw_vina_score": -5.0,
            "best_eligible_vina_score": -6.0,
            "best_eligible_pose_id": "X_conf01_m1",
            "median_eligible_vina_score": -5.0,
            "mean_eligible_vina_score": -5.1,
            "q25_eligible_vina_score": -5.5,
            "q75_eligible_vina_score": -4.5,
            "mean_conformer_median_eligible_vina_score": -5.0,
            "best_eligible_vina_score_range_across_conformers": 0.4,
        }

    def _row(self, conformer_count: int, warnings: list[str]) -> dict[str, object]:
        spec = {
            "analysis_id": "X_v04",
            "candidate": "X",
            "sequence": "AAA",
            "peptide_length": 3,
            "group": "A",
            "template_id": "template",
            "receptor_id": "receptor",
            "screening_protocol_status": "formal_standardized",
            "comparison_protocol_id": "protocol",
            "vina_score_primary_source": "docking_output_pdbqt_remark",
        }
        total = conformer_count * 2
        summary = {
            "schema_version": "0.4",
            "pose_counts": {
                "total": total,
                "by_conformer": {
                    f"conf{index:02d}": 2
                    for index in range(1, conformer_count + 1)
                },
            },
            "conformer_consistency": {
                "recurrent_direct_shared_core": [154, 154, 159]
            },
            "representative_poses": [{"pose_id": "X_conf01_m1"}],
            "input_files": self._inputs(),
            "warnings": warnings,
        }
        return build_candidate_summary_row(
            spec, summary, [self._score_row(total, total - 1)]
        )

    def test_conformer_count_uses_actual_ensemble_inventory(self) -> None:
        rows = [self._row(count, []) for count in (1, 2, 4)]
        self.assertTrue(all("conformer_count" in row for row in rows))
        self.assertEqual([1, 2, 4], [row["conformer_count"] for row in rows])

    def test_derived_counts_and_fraction_follow_frozen_definitions(self) -> None:
        row = self._row(2, [])
        self.assertTrue(
            {
                "eligible_pose_fraction",
                "recurrent_direct_residue_count",
                "representative_count",
            }.issubset(row)
        )
        self.assertEqual(3 / 4, row["eligible_pose_fraction"])
        self.assertEqual(2, row["recurrent_direct_residue_count"])
        self.assertEqual(1, row["representative_count"])

    def test_analysis_status_is_a_stable_success_enum(self) -> None:
        clean = self._row(1, [])
        warned = self._row(1, ["EXPECTED_WARNING"])
        self.assertIn("analysis_status", clean)
        self.assertIn("analysis_status", warned)
        self.assertEqual("PASS", clean["analysis_status"])
        self.assertEqual(
            "PASS_WITH_WARNINGS", warned["analysis_status"]
        )

    def test_required_columns_are_frozen(self) -> None:
        self.assertTrue(CANDIDATE_SUMMARY_REQUIRED_COLUMNS.issubset(CANDIDATE_COLUMNS))


if __name__ == "__main__":
    unittest.main()
