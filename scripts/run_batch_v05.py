"""Plan and preflight a V0.5 candidate subset without running Vina docking."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from batch_preflight_v05 import run_dry_run


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
    parser.add_argument("--dry-run", action="store_true", required=True)
    args = parser.parse_args(argv)

    result = run_dry_run(
        args.manifest,
        args.protocol,
        args.data_root,
        args.output_root,
        args.vina,
        groups=args.groups,
        candidate_ids=args.candidates,
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return {
        "DRY_RUN_PASS": 0,
        "DRY_RUN_FAIL_GLOBAL": 2,
        "DRY_RUN_FAIL_CANDIDATE": 3,
    }[result["outcome"]]


if __name__ == "__main__":
    raise SystemExit(main())
