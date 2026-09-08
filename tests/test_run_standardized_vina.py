from __future__ import annotations

import sys
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "scripts"))

from run_standardized_vina import build_vina_command, inspect_vina_output_text  # noqa: E402


class RunStandardizedVinaTests(unittest.TestCase):
    def test_command_uses_one_config_and_output(self) -> None:
        command = build_vina_command(
            Path("vina.exe"), Path("configs/B3_conf01.txt"), Path("outputs/B3_conf01_out.pdbqt")
        )
        self.assertEqual(
            [
                "vina.exe",
                "--config",
                str(Path("configs/B3_conf01.txt")),
                "--out",
                str(Path("outputs/B3_conf01_out.pdbqt")),
            ],
            command,
        )

    def test_output_audit_uses_actual_model_labels_and_requires_scores(self) -> None:
        text = (
            "MODEL 1\nREMARK VINA RESULT: -6.100 0 0\nATOM      1  C   UNL     1       0.000   0.000   0.000  1.00  0.00     0.0 C\nENDMDL\n"
            "MODEL 3\nREMARK VINA RESULT: -5.900 1 2\nATOM      1  C   UNL     1       1.000   0.000   0.000  1.00  0.00     0.0 C\nENDMDL\n"
        )
        audit = inspect_vina_output_text(text)
        self.assertEqual(2, audit["actual_model_count"])
        self.assertEqual([1, 3], audit["model_labels"])
        self.assertEqual({1: -6.1, 3: -5.9}, audit["vina_scores"])

        with self.assertRaisesRegex(ValueError, "REMARK VINA RESULT"):
            inspect_vina_output_text("MODEL 1\nATOM      1  C\nENDMDL\n")


if __name__ == "__main__":
    unittest.main()
