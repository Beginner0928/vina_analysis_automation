from __future__ import annotations

import contextlib
import io
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from run_batch_v05 import main  # noqa: E402


BASE = [
    "--manifest", str(ROOT / "specs" / "top40_candidates_v05.json"),
    "--protocol", str(ROOT / "specs" / "screening_protocol_v05.json"),
    "--data-root", str(ROOT.parent),
    "--output-root", str(ROOT / "test_output" / "v05b_cli_synthetic"),
    "--vina", str(ROOT / "test_output" / "never_execute_vina.exe"),
    "--candidates", "A4",
]


class BatchExecutionCliTests(unittest.TestCase):
    def test_no_mode_is_rejected_without_calling_either_path(self) -> None:
        with patch("run_batch_v05.run_dry_run") as dry_run:
            with patch("run_batch_v05.execute_batch") as execute:
                with contextlib.redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit):
                        main(BASE)
        dry_run.assert_not_called()
        execute.assert_not_called()

    def test_explicit_execute_routes_to_executor(self) -> None:
        expected = {
            "outcome": "BATCH_COMPLETE",
            "candidates": [],
            "vina_docking_process_count": 0,
        }
        stdout = io.StringIO()
        with patch("run_batch_v05.execute_batch", return_value=expected) as execute:
            with contextlib.redirect_stdout(stdout):
                code = main(BASE + ["--execute"])
        self.assertEqual(0, code)
        self.assertEqual(expected, json.loads(stdout.getvalue()))
        self.assertFalse(execute.call_args.kwargs["resume"])

    def test_resume_is_forwarded_only_in_execution_mode(self) -> None:
        expected = {
            "outcome": "BATCH_COMPLETE_WITH_FAILURES",
            "candidates": [],
            "vina_docking_process_count": 0,
        }
        with patch("run_batch_v05.execute_batch", return_value=expected) as execute:
            with contextlib.redirect_stdout(io.StringIO()):
                code = main(BASE + ["--resume"])
        self.assertEqual(4, code)
        self.assertTrue(execute.call_args.kwargs["resume"])

    def test_modes_are_pairwise_mutually_exclusive(self) -> None:
        for pair in (
            ["--dry-run", "--execute"],
            ["--dry-run", "--resume"],
            ["--execute", "--resume"],
        ):
            with self.subTest(pair=pair):
                with contextlib.redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit):
                        main(BASE + pair)

    def test_dry_run_preserves_v05a_routing(self) -> None:
        expected = {"outcome": "DRY_RUN_PASS", "vina_docking_process_count": 0}
        with patch("run_batch_v05.run_dry_run", return_value=expected) as dry_run:
            with patch("run_batch_v05.execute_batch") as execute:
                with contextlib.redirect_stdout(io.StringIO()):
                    code = main(BASE + ["--dry-run"])
        self.assertEqual(0, code)
        dry_run.assert_called_once()
        execute.assert_not_called()

    def test_force_option_does_not_exist(self) -> None:
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                main(BASE + ["--force"])


if __name__ == "__main__":
    unittest.main()
