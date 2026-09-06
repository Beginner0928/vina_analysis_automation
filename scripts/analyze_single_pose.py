"""Headless V0.1 reproduction of the manually verified A4_conf02 model 1 pose."""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Callable

import pymol
from pymol import cmd


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
GOLDEN_PATH = REPOSITORY_ROOT / "templates" / "A4_conf02_m1_golden.json"
SCRIPT_PATH = Path(__file__).resolve()

MEASUREMENT_COLUMNS = (
    "interaction",
    "ligand_residue",
    "ligand_atom",
    "receptor_residue",
    "receptor_atom",
    "distance_A",
    "angle_deg",
    "interpretation_source",
)

OUTPUT_NAMES = (
    "A4_conf02_m1_measurements.tsv",
    "A4_conf02_m1_contacts.tsv",
    "A4_conf02_m1_metrics.json",
    "A4_conf02_m1_summary.md",
    "run_manifest.txt",
    "A4_conf02_m1_interface.pse",
    "A4_conf02_m1_interface.png",
)


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_input_files(golden: dict[str, Any]) -> dict[str, Path]:
    resolved: dict[str, Path] = {}
    for label, raw_path in golden["input_files"].items():
        path = Path(raw_path).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Required input is missing ({label}): {path}")
        resolved[label] = path
    return resolved


def parse_metric(
    row: dict[str, str], header: set[str], field: str, converter: Callable[[str], Any]
) -> Any:
    if field not in header:
        return "NA"
    raw = row.get(field, "").strip()
    if raw == "":
        return "NA"
    try:
        return converter(raw)
    except (TypeError, ValueError):
        return "NA"


def parse_semicolon_integers(raw: str) -> list[int]:
    if not raw.strip():
        return []
    return [int(value) for value in raw.split(";") if value.strip()]


def read_model_metrics(path: Path, model: int) -> tuple[list[str], dict[str, Any]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames
        if not header or len(header) < 2:
            raise ValueError(f"Pose metrics is not a valid TSV: {path}")
        if "model" not in header:
            raise ValueError(f"Pose metrics header has no 'model' field: {header}")
        model_row = None
        for row in reader:
            if row.get("model", "").strip() == str(model):
                model_row = row
                break
    if model_row is None:
        raise ValueError(f"MODEL {model} is absent from pose metrics: {path}")

    fields = set(header)
    metrics = {
        "model": parse_metric(model_row, fields, "model", int),
        "score": parse_metric(model_row, fields, "score", float),
        "contact_fraction": parse_metric(
            model_row, fields, "contact_fraction", float
        ),
        "contacted_peptide_residues": parse_metric(
            model_row, fields, "contacted_peptide_residues", parse_semicolon_integers
        ),
        "receptor_contact_residues": parse_metric(
            model_row, fields, "receptor_contact_residues", parse_semicolon_integers
        ),
        "target_contacts": parse_metric(
            model_row, fields, "target_contacts", parse_semicolon_integers
        ),
        "target_contact_count": parse_metric(
            model_row, fields, "target_contact_count", int
        ),
        "receptor_clash_pairs": parse_metric(
            model_row, fields, "receptor_clash_pairs", int
        ),
        "min_receptor_distance": parse_metric(
            model_row, fields, "min_receptor_distance", float
        ),
        "chain_b_clash_pairs": parse_metric(
            model_row, fields, "chain_b_clash_pairs", int
        ),
        "min_chain_b_distance": parse_metric(
            model_row, fields, "min_chain_b_distance", float
        ),
        "basic_geometry_pass": parse_metric(
            model_row, fields, "basic_geometry_pass", int
        ),
        "raw_model_row": {name: model_row.get(name, "") for name in header},
    }
    return list(header), metrics


def require_single_atom(selection: str) -> None:
    count = cmd.count_atoms(selection)
    if count != 1:
        raise ValueError(
            f"Expected exactly one atom for selection {selection!r}, found {count}"
        )


def load_vina_multimodel_pdbqt(path: Path, object_name: str) -> int:
    """Load every Vina MODEL into one PyMOL object without rewriting the PDBQT."""
    text = path.read_text(encoding="utf-8")
    blocks = re.findall(r"(?ms)^MODEL\s+\d+\s*$.*?^ENDMDL\s*$", text)
    if not blocks:
        raise ValueError(f"No MODEL/ENDMDL blocks found in Vina output: {path}")
    for state, block in enumerate(blocks, start=1):
        cmd.read_pdbstr(block + "\n", object_name, state=state)
    return len(blocks)


def collect_ca_residues(selection: str) -> dict[int, str]:
    observed: list[tuple[str, str]] = []
    cmd.iterate(
        f"({selection}) and name CA",
        "observed.append((resn, resi))",
        space={"observed": observed},
    )
    result: dict[int, str] = {}
    for resn, resi in observed:
        try:
            number = int(resi)
        except ValueError as exc:
            raise ValueError(f"Non-integer receptor residue identifier: {resi}") from exc
        label = f"{resn}{number}"
        if number in result and result[number] != label:
            raise ValueError(f"Ambiguous receptor residue {number}: {result[number]}, {label}")
        result[number] = label
    return result


def calculate_measurements(golden: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    chain = golden["receptor_chain"]
    for expected in golden["measurements"]:
        ligand_selection = (
            f"A4_conf02_m1 and id {int(expected['ligand_atom_id'])}"
        )
        receptor_selection = (
            "TfR1_A and "
            f"chain {chain} and resi {int(expected['receptor_resi'])} "
            f"and name {expected['receptor_atom']}"
        )
        require_single_atom(ligand_selection)
        require_single_atom(receptor_selection)
        distance = cmd.get_distance(ligand_selection, receptor_selection)

        angle = None
        angle_atom_id = expected.get("ligand_angle_atom_id")
        if angle_atom_id is not None:
            angle_selection = f"A4_conf02_m1 and id {int(angle_atom_id)}"
            require_single_atom(angle_selection)
            angle = cmd.get_angle(
                ligand_selection, angle_selection, receptor_selection
            )

        rows.append(
            {
                "interaction": expected["interaction"],
                "ligand_residue": expected["ligand_residue"],
                "ligand_atom": expected["ligand_atom"],
                "receptor_residue": expected["receptor_residue"],
                "receptor_atom": expected["receptor_atom"],
                "distance_A": distance,
                "angle_deg": angle,
                "interpretation_source": golden["interpretation_source"],
            }
        )
    return rows


def create_measurement_objects(golden: dict[str, Any]) -> None:
    chain = golden["receptor_chain"]
    object_names: list[str] = []
    for index, expected in enumerate(golden["measurements"], start=1):
        ligand_selection = (
            f"A4_conf02_m1 and id {int(expected['ligand_atom_id'])}"
        )
        receptor_selection = (
            "TfR1_A and "
            f"chain {chain} and resi {int(expected['receptor_resi'])} "
            f"and name {expected['receptor_atom']}"
        )
        distance_name = f"distance_{index}_{expected['interaction'].replace('-', '_')}"
        cmd.distance(distance_name, ligand_selection, receptor_selection)
        object_names.append(distance_name)
        angle_atom_id = expected.get("ligand_angle_atom_id")
        if angle_atom_id is not None:
            angle_name = f"angle_{index}_{expected['interaction'].replace('-', '_')}"
            cmd.angle(
                angle_name,
                ligand_selection,
                f"A4_conf02_m1 and id {int(angle_atom_id)}",
                receptor_selection,
            )
            object_names.append(angle_name)

    for name in object_names:
        if not name.startswith("distance_1_"):
            cmd.disable(name)


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            output = dict(row)
            if "distance_A" in output:
                output["distance_A"] = f"{float(output['distance_A']):.6f}"
            if "angle_deg" in output:
                value = output["angle_deg"]
                output["angle_deg"] = "" if value is None else f"{float(value):.6f}"
            writer.writerow({column: output.get(column, "") for column in columns})


def format_metric(value: Any, decimals: int = 6) -> str:
    if isinstance(value, float):
        return f"{value:.{decimals}f}"
    if isinstance(value, list):
        return ", ".join(str(item) for item in value)
    return str(value)


def build_summary(
    golden: dict[str, Any],
    state_count: int,
    actual_contacts: list[str],
    metrics: dict[str, Any],
    measurements: list[dict[str, Any]],
) -> str:
    expected_contacts = set(golden["expected_target_contacts"])
    contact_status = "PASS" if set(actual_contacts) == expected_contacts else "FAIL"
    lines = [
        "# A4_conf02 MODEL 1 automated reproduction — V0.1",
        "",
        "## Objective results",
        "",
        f"- Candidate: `{golden['candidate']}`",
        f"- Sequence: `{golden['sequence']}`",
        f"- Conformer: `{golden['conformer']}`",
        f"- Model: `{golden['model']}`",
        f"- State count: `{state_count}`",
        f"- Vina score: `{format_metric(metrics['score'], 3)} kcal/mol`",
        f"- Contact fraction: `{format_metric(metrics['contact_fraction'], 6)}`",
        "- Target contacts: `" + ", ".join(actual_contacts) + "`",
        f"- Golden contact-set comparison: `{contact_status}`",
        f"- Receptor clash pairs: `{format_metric(metrics['receptor_clash_pairs'])}`",
        f"- Minimum receptor distance: `{format_metric(metrics['min_receptor_distance'])} Å`",
        f"- Chain B clash pairs: `{format_metric(metrics['chain_b_clash_pairs'])}`",
        f"- Minimum chain B distance: `{format_metric(metrics['min_chain_b_distance'])} Å`",
        f"- Basic geometry pass: `{format_metric(metrics['basic_geometry_pass'])}`",
        "",
        "## Measurements",
        "",
        "| Interaction | Distance (Å) | Angle (°) | Manual confirmed interpretation |",
        "|---|---:|---:|---|",
    ]
    interpretations = {
        row["interaction"]: row["interpretation"]
        for row in golden["measurements"]
    }
    for row in measurements:
        angle = "" if row["angle_deg"] is None else f"{row['angle_deg']:.6f}"
        lines.append(
            f"| {row['interaction']} | {row['distance_A']:.6f} | {angle} | "
            f"{interpretations[row['interaction']]} |"
        )
    lines.extend(
        [
            "",
            "The interpretations above are copied from the manual golden reference; "
            "they are not inferred from the 4 Å contact calculation.",
            "",
        ]
    )
    return "\n".join(lines)


def build_manifest(
    golden: dict[str, Any],
    inputs: dict[str, Path],
    input_hashes: dict[str, str],
    output_dir: Path,
) -> str:
    version = cmd.get_version()[0]
    lines = [
        f"run_datetime={datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"pymol_version={version}",
        f"python_version={platform.python_version()}",
        f"python_executable={Path(sys.executable).resolve()}",
        f"golden_json={GOLDEN_PATH.resolve()}",
        f"golden_json_sha256={sha256(GOLDEN_PATH)}",
        f"script={SCRIPT_PATH}",
        f"script_sha256={sha256(SCRIPT_PATH)}",
        f"case_id={golden['case_id']}",
    ]
    for label, path in inputs.items():
        lines.append(f"input.{label}.path={path}")
        lines.append(f"input.{label}.sha256={input_hashes[label]}")
    for name in OUTPUT_NAMES:
        lines.append(f"output={output_dir / name}")
    return "\n".join(lines) + "\n"


def configure_scene(golden: dict[str, Any]) -> None:
    cmd.disable("A4_conf02_all")
    cmd.hide("everything", "all")
    cmd.show("cartoon", "TfR1_A")
    cmd.color("gray80", "TfR1_A")
    cmd.show("sticks", "target_site")
    cmd.color("orange", "target_site")
    cmd.show("sticks", "A4_conf02_m1")
    cmd.color("magenta", "A4_conf02_m1")
    cmd.show("sticks", "actual_target_contacts")
    cmd.label("actual_target_contacts and name CA", '"%s%s" % (resn, resi)')
    cmd.set("label_color", "black")
    cmd.set("label_size", 20)
    cmd.set("label_font_id", 7)
    cmd.set("dash_color", "tv_blue")
    cmd.set("dash_width", 4.0)
    cmd.set("dash_gap", 0.25)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    focus = "A4_conf02_m1 or actual_target_contacts"
    cmd.orient(focus)
    cmd.zoom(focus, buffer=7.0, complete=1)
    create_measurement_objects(golden)


def run_analysis() -> Path:
    golden = load_json(GOLDEN_PATH)
    output_dir = Path(golden["output_directory"]).resolve()
    if output_dir.exists():
        raise FileExistsError(
            f"Output directory already exists; refusing to overwrite: {output_dir}"
        )

    inputs = require_input_files(golden)
    input_hashes = {label: sha256(path) for label, path in inputs.items()}
    metrics_header, pose_metrics = read_model_metrics(
        inputs["pose_metrics"], int(golden["model"])
    )

    pymol.finish_launching(["pymol", "-cq"])
    try:
        cmd.load(str(inputs["receptor_chain_a"]), "TfR1_A")
        source_model_count = load_vina_multimodel_pdbqt(
            inputs["vina_output"], "A4_conf02_all"
        )
        state_count = cmd.count_states("A4_conf02_all")
        if source_model_count != int(golden["expected_state_count"]):
            raise ValueError(
                f"Expected {golden['expected_state_count']} MODEL blocks, "
                f"found {source_model_count}"
            )
        if state_count != int(golden["expected_state_count"]):
            raise ValueError(
                f"Expected {golden['expected_state_count']} ligand states, found {state_count}"
            )

        cmd.create("A4_conf02_m1", "A4_conf02_all", 1, 1)
        target_expression = "+".join(str(value) for value in golden["target_residues"])
        target_atom_count = cmd.select(
            "target_site",
            f"TfR1_A and chain {golden['receptor_chain']} and resi {target_expression}",
        )
        if target_atom_count == 0:
            raise ValueError("target_site selection is empty")

        target_labels = collect_ca_residues("target_site")
        missing_target_residues = sorted(
            set(int(value) for value in golden["target_residues"])
            - set(target_labels)
        )
        if missing_target_residues:
            raise ValueError(
                f"Target residues have no CA atom in receptor: {missing_target_residues}"
            )

        cmd.select(
            "actual_target_contacts",
            "byres (target_site within "
            f"{float(golden['contact_cutoff_A']):.3f} of A4_conf02_m1)",
        )
        contact_labels_by_number = collect_ca_residues("actual_target_contacts")
        actual_contacts = [
            contact_labels_by_number[number]
            for number in (int(value) for value in golden["target_residues"])
            if number in contact_labels_by_number
        ]
        measurements = calculate_measurements(golden)

        contact_rows = [
            {
                "target_residue": target_labels[number],
                "within_4A": 1 if number in contact_labels_by_number else 0,
            }
            for number in (int(value) for value in golden["target_residues"])
        ]
        contact_set_status = (
            "PASS"
            if set(actual_contacts) == set(golden["expected_target_contacts"])
            else "FAIL"
        )
        result_metrics = {
            "schema_version": "0.1",
            "case_id": golden["case_id"],
            "candidate": golden["candidate"],
            "conformer": golden["conformer"],
            "model": golden["model"],
            "state_count": state_count,
            "target_site_atom_count": target_atom_count,
            "actual_target_contacts": actual_contacts,
            "golden_target_contacts": golden["expected_target_contacts"],
            "target_contact_set_status": contact_set_status,
            "pose_metrics_header": metrics_header,
            "pose_metrics": pose_metrics,
            "measurements": measurements,
        }

        output_dir.mkdir(parents=False, exist_ok=False)
        write_tsv(
            output_dir / "A4_conf02_m1_measurements.tsv",
            MEASUREMENT_COLUMNS,
            measurements,
        )
        write_tsv(
            output_dir / "A4_conf02_m1_contacts.tsv",
            ("target_residue", "within_4A"),
            contact_rows,
        )
        (output_dir / "A4_conf02_m1_metrics.json").write_text(
            json.dumps(result_metrics, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        (output_dir / "A4_conf02_m1_summary.md").write_text(
            build_summary(
                golden, state_count, actual_contacts, pose_metrics, measurements
            ),
            encoding="utf-8",
            newline="\n",
        )

        configure_scene(golden)
        cmd.save(str(output_dir / "A4_conf02_m1_interface.pse"))
        cmd.png(
            str(output_dir / "A4_conf02_m1_interface.png"),
            width=2400,
            height=1800,
            dpi=300,
            ray=1,
            quiet=0,
        )
        (output_dir / "run_manifest.txt").write_text(
            build_manifest(golden, inputs, input_hashes, output_dir),
            encoding="utf-8",
            newline="\n",
        )
    finally:
        cmd.quit()

    return output_dir


def main() -> int:
    try:
        output_dir = run_analysis()
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print(f"Analysis completed: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
