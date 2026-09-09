"""Generate the fixed A4/B3 competition-conformer software smoke package (no docking)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from competition_conformer_core import (
    generate_protocol_package,
    load_authoritative_generation_inputs,
)


SMOKE_CANDIDATES = ("A4", "B3")
PROJECT_ROOT = Path(__file__).resolve().parents[1]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument(
        "--smoke-root",
        type=Path,
        required=True,
        help="Fresh external directory under which the protocol package is created.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        smoke_root = args.smoke_root.resolve()
        try:
            smoke_root.relative_to(PROJECT_ROOT)
        except ValueError:
            pass
        else:
            raise ValueError("Smoke output must remain outside the Git repository")
        (
            protocol,
            protocol_sha256,
            candidates,
            manifest_sha256,
            chemistry_contract,
            chemistry_sha256,
        ) = load_authoritative_generation_inputs(args.protocol.resolve(), args.data_root.resolve())
        by_id = {candidate["candidate_id"]: candidate for candidate in candidates}
        selected = [by_id[candidate_id] for candidate_id in SMOKE_CANDIDATES]
        protocol_directory = smoke_root / protocol["protocol_id"]
        manifest = generate_protocol_package(
            selected,
            protocol,
            protocol_sha256,
            manifest_sha256,
            chemistry_contract,
            chemistry_sha256,
            protocol_directory,
        )
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print(f"Generated A4/B3 smoke package: {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
