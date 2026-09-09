"""Dry-run, execute, or resume a V0.5 candidate subset with an explicit mode."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from batch_preflight_v05 import run_dry_run
from batch_execution_v05 import execute_batch


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--vina", type=Path, required=True)
    selectors = parser.add_mutually_exclusive_group(required=True)
    selectors.add_argument("--groups", nargs="+")
    selectors.add_argument("--candidates", nargs="+")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true")
    mode.add_argument("--execute", action="store_true")
    mode.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)

    function = run_dry_run if args.dry_run else execute_batch
    keywords = {
        "groups": args.groups,
        "candidate_ids": args.candidates,
    }
    if args.execute or args.resume:
        keywords["resume"] = args.resume
    result = function(
        args.manifest,
        args.protocol,
        args.data_root,
        args.output_root,
        args.vina,
        **keywords,
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return {
        "DRY_RUN_PASS": 0,
        "DRY_RUN_FAIL_GLOBAL": 2,
        "DRY_RUN_FAIL_CANDIDATE": 3,
        "BATCH_COMPLETE": 0,
        "BATCH_FAIL_GLOBAL": 2,
        "BATCH_COMPLETE_WITH_FAILURES": 4,
    }[result["outcome"]]


if __name__ == "__main__":
    raise SystemExit(main())
