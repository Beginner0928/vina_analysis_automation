from __future__ import annotations

import io
import json
import sys
import unittest
import uuid
from pathlib import Path
from unittest.mock import patch


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "scripts"))

from run_standardized_vina import (  # noqa: E402
    build_vina_command,
    inspect_vina_output_text,
    run_one,
)
from standardized_vina_inputs import sha256  # noqa: E402


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

    def _prepared_root(self) -> tuple[Path, Path]:
        root = REPOSITORY_ROOT / "test_output" / f"v05b_fake_vina_{uuid.uuid4().hex}"
        run_id = "X1_conf01"
        config = root / "configs" / f"{run_id}.txt"
        ligand = root / "ligands" / f"{run_id}.pdbqt"
        audit = root / "audit" / f"{run_id}_input_audit.json"
        vina = root / "fake_vina.exe"
        for path, text in (
            (config, "receptor = receptor.pdbqt\n"),
            (ligand, "ATOM      1  C\n"),
            (vina, "synthetic executable\n"),
        ):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding="utf-8", newline="\n")
        audit.parent.mkdir(parents=True, exist_ok=True)
        audit.write_text(
            json.dumps(
                {
                    "run_id": run_id,
                    "config": {"sha256": sha256(config)},
                    "generated_pdbqt": {"sha256": sha256(ligand)},
                }
            )
            + "\n",
            encoding="utf-8",
            newline="\n",
        )
        (root / "protocol.json").write_text(
            json.dumps({"docking_protocol": {"version": "1.2.7"}}) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        return root, vina

    def _run_fake(self, output_text: str | None, exit_code: int):
        root, vina = self._prepared_root()
        calls = []

        class FakeProcess:
            def __init__(self, command, **kwargs):
                calls.append((command, kwargs))
                if output_text is not None:
                    output = Path(kwargs["cwd"]) / command[-1]
                    output.parent.mkdir(parents=True, exist_ok=True)
                    output.write_text(output_text, encoding="ascii", newline="\n")
                self.stdout = io.StringIO("synthetic Vina output\n")

            def wait(self):
                return exit_code

        with patch("run_standardized_vina.query_vina_version", return_value="AutoDock Vina v1.2.7"):
            with patch("run_standardized_vina.subprocess.Popen", side_effect=FakeProcess):
                status_path, complete = run_one(root, "X1", "conf01", vina)
        return json.loads(status_path.read_text(encoding="utf-8")), complete, calls

    def test_fake_vina_success_uses_argv_without_shell_and_records_models(self) -> None:
        output = (
            "MODEL 1\nREMARK VINA RESULT: -6.100 0 0\n"
            "ATOM      1  C   UNL     1       0.000   0.000   0.000  1.00  0.00     0.0 C\nENDMDL\n"
        )
        status, complete, calls = self._run_fake(output, 0)
        self.assertTrue(complete)
        self.assertEqual([1], status["model_labels"])
        self.assertEqual(0, status["exit_code"])
        self.assertNotIn("shell", calls[0][1])
        self.assertIsInstance(calls[0][0], list)

    def test_fake_vina_nonzero_missing_and_malformed_outputs_do_not_complete(self) -> None:
        cases = (
            (None, 7, "exited with code 7"),
            (None, 0, "missing or empty"),
            ("not a PDBQT model\n", 0, "no complete MODEL"),
        )
        for output, exit_code, message in cases:
            with self.subTest(exit_code=exit_code, output=output):
                status, complete, _ = self._run_fake(output, exit_code)
                self.assertFalse(complete)
                self.assertIn(message, status["error"])


if __name__ == "__main__":
    unittest.main()
