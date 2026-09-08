"""Run one prepared standardized Vina job and record its actual output inventory."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from peptide_score_summary import VINA_RESULT_PATTERN
from standardized_vina_inputs import sha256


def build_vina_command(vina: Path, config: Path, output: Path) -> list[str]:
    return [str(vina), "--config", str(config), "--out", str(output)]


def inspect_vina_output_text(text: str) -> dict[str, Any]:
    blocks = re.findall(r"(?ms)^MODEL\s+\d+\s*$.*?^ENDMDL\s*$", text)
    if not blocks:
        raise ValueError("Output PDBQT contains no complete MODEL/ENDMDL blocks")
    labels: list[int] = []
    scores: dict[int, float] = {}
    atom_counts: dict[int, int] = {}
    for block in blocks:
        label = int(block.splitlines()[0].split()[1])
        if label in scores:
            raise ValueError(f"Duplicate MODEL label {label} in output PDBQT")
        matches = [
            match
            for line in block.splitlines()
            if (match := VINA_RESULT_PATTERN.match(line)) is not None
        ]
        if len(matches) != 1:
            raise ValueError(
                f"MODEL {label} must contain exactly one REMARK VINA RESULT; found {len(matches)}"
            )
        atom_count = sum(
            line.startswith(("ATOM  ", "HETATM")) for line in block.splitlines()
        )
        if atom_count == 0:
            raise ValueError(f"MODEL {label} contains no atom records")
        labels.append(label)
        scores[label] = float(matches[0].group(1))
        atom_counts[label] = atom_count
    return {
        "actual_model_count": len(labels),
        "model_labels": labels,
        "vina_scores": scores,
        "atom_record_counts": atom_counts,
        "remark_vina_result_present_for_all_models": True,
        "pdbqt_parse_status": "PASS",
    }


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def query_vina_version(vina: Path) -> str:
    completed = subprocess.run(
        [str(vina), "--version"], capture_output=True, text=True, check=False
    )
    text = (completed.stdout + completed.stderr).strip()
    if completed.returncode != 0 or not text:
        raise RuntimeError(f"Could not query Vina version from {vina}: {text}")
    return text.splitlines()[0]


def run_one(root: Path, candidate: str, conformer: str, vina: Path) -> tuple[Path, bool]:
    run_id = f"{candidate}_{conformer}"
    config = root / "configs" / f"{run_id}.txt"
    input_audit_path = root / "audit" / f"{run_id}_input_audit.json"
    ligand = root / "ligands" / f"{run_id}.pdbqt"
    output = root / "outputs" / f"{run_id}_out.pdbqt"
    log_path = root / "logs" / f"{run_id}.log"
    status_path = root / "run_status" / f"{run_id}.json"
    for path, label in ((config, "config"), (input_audit_path, "input audit"), (ligand, "ligand")):
        if not path.is_file():
            raise FileNotFoundError(f"Prepared {label} is missing: {path}")
    for path in (output, log_path, status_path):
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite existing run artifact: {path}")

    input_audit = load_json(input_audit_path)
    if input_audit.get("run_id") != run_id:
        raise ValueError(f"Input audit run_id does not match {run_id}")
    if sha256(config) != input_audit["config"]["sha256"]:
        raise ValueError(f"Config hash changed after input audit: {config}")
    if sha256(ligand) != input_audit["generated_pdbqt"]["sha256"]:
        raise ValueError(f"Ligand PDBQT hash changed after input audit: {ligand}")

    root.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    status_path.parent.mkdir(parents=True, exist_ok=True)
    vina_version = query_vina_version(vina)
    protocol = load_json(root / "protocol.json")
    expected_version = str(protocol["docking_protocol"]["version"])
    if expected_version not in vina_version:
        raise ValueError(
            f"Vina version mismatch: protocol expects {expected_version}, executable reports {vina_version}"
        )
    command = build_vina_command(
        vina.resolve(), config.relative_to(root), output.relative_to(root)
    )
    started = datetime.now(timezone.utc)
    start_clock = time.monotonic()
    with log_path.open("x", encoding="utf-8", newline="\n") as log:
        log.write(f"START_UTC={started.isoformat()}\n")
        log.write("COMMAND=" + " ".join(command) + "\n")
        log.write(f"VINA_VERSION={vina_version}\n")
        log.flush()
        process = subprocess.Popen(
            command,
            cwd=root,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
        exit_code = process.wait()
        elapsed = time.monotonic() - start_clock
        log.write(f"EXIT_CODE={exit_code}\n")
        log.write(f"ELAPSED_SECONDS={elapsed:.2f}\n")

    status: dict[str, Any] = {
        "run_id": run_id,
        "start_utc": started.isoformat(),
        "end_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": elapsed,
        "command": command,
        "vina_executable": str(vina.resolve()),
        "vina_executable_sha256": sha256(vina),
        "vina_version": vina_version,
        "exit_code": exit_code,
        "config_path": config.relative_to(root).as_posix(),
        "config_sha256": sha256(config),
        "ligand_pdbqt_path": ligand.relative_to(root).as_posix(),
        "ligand_pdbqt_sha256": sha256(ligand),
        "input_audit_path": input_audit_path.relative_to(root).as_posix(),
        "input_audit_sha256": sha256(input_audit_path),
        "log_path": log_path.relative_to(root).as_posix(),
        "log_sha256": sha256(log_path),
        "output_path": output.relative_to(root).as_posix(),
        "complete": False,
    }
    try:
        if exit_code != 0:
            raise RuntimeError(f"Vina exited with code {exit_code}")
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError("Vina output PDBQT is missing or empty")
        output_audit = inspect_vina_output_text(
            output.read_text(encoding="ascii", errors="replace")
        )
        status.update(output_audit)
        status["output_sha256"] = sha256(output)
        status["complete"] = True
    except Exception as exc:
        status["error"] = str(exc)
    with status_path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(status, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return status_path, bool(status["complete"])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--conformer", required=True)
    parser.add_argument("--vina", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        status_path, complete = run_one(
            args.root.resolve(), args.candidate, args.conformer, args.vina.resolve()
        )
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print(f"Run status: {status_path} (complete={complete})")
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
