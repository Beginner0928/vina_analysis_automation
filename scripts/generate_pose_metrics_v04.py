"""Generate standardized V0.4 pose metrics for one completed Vina output."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pose_metrics_v04 import generate_pose_metrics


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--docking-output", type=Path, required=True)
    parser.add_argument("--source-sdf", type=Path, required=True)
    parser.add_argument("--ligand-pdbqt", type=Path, required=True)
    parser.add_argument("--receptor-monomer", type=Path, required=True)
    parser.add_argument("--receptor-dimer", type=Path, required=True)
    parser.add_argument("--receptor-chain", default="A")
    parser.add_argument("--dimer-chain", default="B")
    parser.add_argument("--sequence", required=True)
    parser.add_argument("--target-residues", required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        rows = generate_pose_metrics(
            docking_output=args.docking_output.resolve(),
            source_sdf=args.source_sdf.resolve(),
            ligand_pdbqt=args.ligand_pdbqt.resolve(),
            receptor_monomer=args.receptor_monomer.resolve(),
            receptor_dimer=args.receptor_dimer.resolve(),
            receptor_chain=args.receptor_chain,
            dimer_chain=args.dimer_chain,
            sequence=args.sequence,
            target_residues={int(value) for value in args.target_residues.split(",")},
            output_tsv=args.output_tsv.resolve(),
        )
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print(f"Pose metrics written: {args.output_tsv.resolve()} ({len(rows)} MODEL rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
