"""Generate the complete canonical Top40 competition conformer package (no docking)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from competition_conformer_core import (
    generate_protocol_package,
    load_authoritative_generation_inputs,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument(
        "--all",
        action="store_true",
        help="Explicitly generate all 40 canonical candidates into the production protocol directory.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        (
            protocol,
            protocol_sha256,
            candidates,
            manifest_sha256,
            chemistry_contract,
            chemistry_sha256,
        ) = load_authoritative_generation_inputs(args.protocol.resolve(), args.data_root.resolve())
        protocol_directory = (
            args.data_root.resolve()
            / "standardized_ligands"
            / protocol["protocol_id"]
        )
        manifest = generate_protocol_package(
            candidates,
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
    print(f"Generated competition conformers: {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
