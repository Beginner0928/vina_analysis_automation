"""Audit an existing competition conformer package without modifying it."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from competition_conformer_core import audit_generation_package_against_contract


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--protocol-directory", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        audit = audit_generation_package_against_contract(
            args.protocol_directory.resolve(),
            args.protocol.resolve(),
            args.data_root.resolve(),
        )
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print(
        f"Competition conformer audit PASS: candidates={audit['candidate_count']} "
        f"SDFs={audit['sdf_count']} provenance={audit['contract_provenance']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
