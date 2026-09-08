"""Validate a completed standardized Vina run before the next serial run begins."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

from peptide_ensemble_core import read_pose_metrics_table, validate_model_inventory
from peptide_score_summary import crosscheck_vina_scores, parse_vina_scores
from pose_analysis_core import parse_vina_model_blocks
from standardized_vina_inputs import sha256


def audit_completed_run(root: Path, candidate: str, conformer: str) -> Path:
    run_id = f"{candidate}_{conformer}"
    paths = {
        "config": root / "configs" / f"{run_id}.txt",
        "ligand_pdbqt": root / "ligands" / f"{run_id}.pdbqt",
        "input_audit": root / "audit" / f"{run_id}_input_audit.json",
        "pre_docking_gate": root / "audit" / f"{run_id}_pre_docking_gate.json",
        "log": root / "logs" / f"{run_id}.log",
        "output_pdbqt": root / "outputs" / f"{run_id}_out.pdbqt",
        "run_status": root / "run_status" / f"{run_id}.json",
        "pose_metrics": root / "analysis" / f"{run_id}_pose_metrics.tsv",
    }
    missing = [name for name, path in paths.items() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Completed-run audit missing artifacts: {missing}")
    status = json.loads(paths["run_status"].read_text(encoding="utf-8"))
    if status.get("exit_code") != 0 or status.get("complete") is not True:
        raise ValueError(f"Vina run status is not successful: {status}")
    pre_gate = json.loads(paths["pre_docking_gate"].read_text(encoding="utf-8"))
    if pre_gate.get("gate_status") != "PASS":
        raise ValueError("Pre-docking gate did not pass")
    log_text = paths["log"].read_text(encoding="utf-8", errors="replace")
    if "EXIT_CODE=0" not in log_text:
        raise ValueError("Vina log does not record EXIT_CODE=0")

    blocks = parse_vina_model_blocks(paths["output_pdbqt"])
    scores = parse_vina_scores(blocks)
    metrics_header, metrics = read_pose_metrics_table(paths["pose_metrics"])
    validate_model_inventory(conformer, set(blocks), set(metrics), len(blocks))
    score_checks = crosscheck_vina_scores(scores, metrics, 0.001)
    if len(blocks) != int(status["actual_model_count"]):
        raise ValueError("Run status actual MODEL count differs from parsed output")
    if list(blocks) != [int(value) for value in status["model_labels"]]:
        raise ValueError("Run status MODEL labels differ from parsed output")
    if sha256(paths["output_pdbqt"]) != status["output_sha256"]:
        raise ValueError("Output PDBQT hash changed after run status was written")

    result = {
        "audit_datetime": datetime.now().astimezone().isoformat(timespec="seconds"),
        "gate": "POST_DOCKING_CONF01_PILOT_GATE" if conformer == "conf01" else "POST_DOCKING_GATE",
        "gate_status": "PASS",
        "run_id": run_id,
        "vina_normal_exit": True,
        "pdbqt_parse_status": "PASS",
        "remark_vina_result_status": "PASS",
        "actual_model_count": len(blocks),
        "model_labels": list(blocks),
        "score_parser": "V0.4 docking-output PDBQT REMARK parser",
        "score_crosscheck_status_counts": {
            "MATCH": sum(row["status"] == "MATCH" for row in score_checks.values())
        },
        "score_crosscheck_tolerance_kcal_mol": 0.001,
        "pose_metrics_row_count": len(metrics),
        "pose_metrics_header": metrics_header,
        "hashes": {name: sha256(path) for name, path in paths.items()},
        "configuration_and_provenance_complete": True,
        "run_status": status,
    }
    output = root / "audit" / f"{run_id}_post_docking_gate.json"
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite completed-run audit: {output}")
    with output.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(result, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--conformer", required=True)
    args = parser.parse_args(argv)
    try:
        output = audit_completed_run(
            args.root.resolve(), args.candidate, args.conformer
        )
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print(f"Post-docking gate PASS: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
