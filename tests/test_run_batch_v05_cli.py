from __future__ import annotations

import contextlib
import io
import json
import sys
import unittest
import uuid
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT.parent
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from run_batch_v05 import main  # noqa: E402


BASE = [
    "--manifest", str(ROOT / "specs" / "top40_candidates_v05.json"),
    "--protocol", str(ROOT / "specs" / "screening_protocol_v05.json"),
    "--data-root", str(DATA_ROOT),
    "--vina", str(DATA_ROOT / "tools" / "AutoDockVina" / "vina.exe"),
    "--dry-run",
]


class BatchCliTests(unittest.TestCase):
    def test_cli_emits_machine_readable_pass_for_explicit_ready_subset(self) -> None:
        output = ROOT / "test_output" / f"v05a_cli_{uuid.uuid4().hex}"
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            code = main(BASE + ["--output-root", str(output), "--candidates", "B3", "A4"])
        result = json.loads(stdout.getvalue())
        self.assertEqual(0, code)
        self.assertEqual("DRY_RUN_PASS", result["outcome"])
        self.assertEqual(["A4", "B3"], result["selected_candidate_ids"])
        self.assertEqual(0, result["vina_docking_process_count"])
        self.assertFalse(output.exists())

    def test_cli_candidate_incompleteness_has_stable_nonzero_exit(self) -> None:
        output = ROOT / "test_output" / f"v05a_cli_ab_{uuid.uuid4().hex}"
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            code = main(BASE + ["--output-root", str(output), "--groups", "A", "B"])
        result = json.loads(stdout.getvalue())
        self.assertEqual(3, code)
        self.assertEqual("DRY_RUN_FAIL_CANDIDATE", result["outcome"])
        self.assertEqual(20, result["selected_candidate_count"])

    def test_cli_requires_exactly_one_selector(self) -> None:
        output = ROOT / "test_output" / f"v05a_cli_selector_{uuid.uuid4().hex}"
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                main(BASE + ["--output-root", str(output)])
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                main(BASE + ["--output-root", str(output), "--groups", "A", "--candidates", "A1"])

    def test_cli_does_not_expose_v05b_flags(self) -> None:
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            with self.assertRaises(SystemExit):
                main(BASE + ["--output-root", str(ROOT / "test_output" / uuid.uuid4().hex), "--candidates", "A4", "--resume"])
        self.assertIn("unrecognized arguments: --resume", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
