"""CLI for V0.2 parameterized headless single-pose analysis."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import pymol
from pymol import cmd

from pose_analysis_core import (
    compare_legacy_contacts,
    compute_target_residue_contacts,
    parse_pdbqt_atoms,
    parse_receptor_pdb_atoms,
    parse_vina_model_blocks,
    read_model_metrics,
    require_model_block,
    resolve_inputs,
    sha256,
    validate_spec,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
TEST_OUTPUT_ROOT = (REPOSITORY_ROOT / "test_output").resolve()
SCRIPT_PATH = Path(__file__).resolve()
CORE_PATH = SCRIPT_PATH.with_name("pose_analysis_core.py")
OUTPUT_NAMES = (
    "metrics.json",
    "contacts.tsv",
    "measurements.tsv",
    "summary.md",
    "interface.pse",
    "interface.png",
    "run_manifest.txt",
)
MEASUREMENT_COLUMNS = (
    "interaction",
    "ligand_residue",
    "ligand_atom",
    "receptor_residue",
    "receptor_atom",
    "distance_A",
    "angle_deg",
    "interpretation",
    "interpretation_source",
    "measurement_provenance",
)


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def require_safe_output_directory(path: Path) -> Path:
    resolved = path.resolve()
    if not resolved.is_relative_to(TEST_OUTPUT_ROOT):
        raise ValueError(f"--output-dir must be inside {TEST_OUTPUT_ROOT}: {resolved}")
    if resolved.exists():
        raise FileExistsError(
            f"Output directory already exists; refusing to overwrite: {resolved}"
        )
    return resolved


def require_single_atom(selection: str) -> None:
    count = cmd.count_atoms(selection)
    if count != 1:
        raise ValueError(
            f"Expected exactly one atom for selection {selection!r}, found {count}"
        )


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


def calculate_measurements(
    spec: dict[str, Any], receptor_labels: dict[int, str]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    chain = spec["receptor_chain"]
    for interaction in spec["manual_interactions"]:
        distance_spec = interaction["distance"]
        ligand_selection = f"ligand_pose and id {int(distance_spec['ligand_atom_id'])}"
        receptor_selection = (
            "receptor_monomer and "
            f"chain {chain} and resi {int(distance_spec['receptor_residue'])} "
            f"and name {distance_spec['receptor_atom']}"
        )
        require_single_atom(ligand_selection)
        require_single_atom(receptor_selection)
        distance = cmd.get_distance(ligand_selection, receptor_selection)

        angle_value = None
        angle_spec = interaction.get("angle")
        if angle_spec is not None:
            donor_selection = (
                f"ligand_pose and id {int(angle_spec['ligand_donor_atom_id'])}"
            )
            hydrogen_selection = (
                f"ligand_pose and id {int(angle_spec['ligand_hydrogen_atom_id'])}"
            )
            acceptor = angle_spec["receptor_acceptor"]
            acceptor_selection = (
                "receptor_monomer and "
                f"chain {chain} and resi {int(acceptor['residue'])} "
                f"and name {acceptor['atom']}"
            )
            require_single_atom(donor_selection)
            require_single_atom(hydrogen_selection)
            require_single_atom(acceptor_selection)
            angle_value = cmd.get_angle(
                donor_selection, hydrogen_selection, acceptor_selection
            )

        labels = interaction.get("labels", {})
        receptor_residue_number = int(distance_spec["receptor_residue"])
        rows.append(
            {
                "interaction": interaction["interaction_id"],
                "ligand_residue": labels.get("ligand_residue", ""),
                "ligand_atom": labels.get("ligand_atom", f"id:{distance_spec['ligand_atom_id']}"),
                "receptor_residue": labels.get(
                    "receptor_residue",
                    receptor_labels.get(receptor_residue_number, str(receptor_residue_number)),
                ),
                "receptor_atom": distance_spec["receptor_atom"],
                "distance_A": distance,
                "angle_deg": angle_value,
                "interpretation": interaction.get("interpretation", ""),
                "interpretation_source": interaction.get(
                    "interpretation_source", "manual_spec"
                ),
                "measurement_provenance": "pymol_v02_recomputed_from_manual_spec",
            }
        )
    return rows


def safe_pymol_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_]", "_", value)
    return cleaned[:80] or "interaction"


def create_measurement_objects(spec: dict[str, Any]) -> None:
    chain = spec["receptor_chain"]
    created: list[str] = []
    for index, interaction in enumerate(spec["manual_interactions"], start=1):
        distance_spec = interaction["distance"]
        ligand = f"ligand_pose and id {int(distance_spec['ligand_atom_id'])}"
        receptor = (
            "receptor_monomer and "
            f"chain {chain} and resi {int(distance_spec['receptor_residue'])} "
            f"and name {distance_spec['receptor_atom']}"
        )
        suffix = safe_pymol_name(interaction["interaction_id"])
        distance_name = f"distance_{index}_{suffix}"
        cmd.distance(distance_name, ligand, receptor)
        created.append(distance_name)
        angle_spec = interaction.get("angle")
        if angle_spec is not None:
            acceptor = angle_spec["receptor_acceptor"]
            angle_name = f"angle_{index}_{suffix}"
            cmd.angle(
                angle_name,
                f"ligand_pose and id {int(angle_spec['ligand_donor_atom_id'])}",
                f"ligand_pose and id {int(angle_spec['ligand_hydrogen_atom_id'])}",
                "receptor_monomer and "
                f"chain {chain} and resi {int(acceptor['residue'])} "
                f"and name {acceptor['atom']}",
            )
            created.append(angle_name)
    for name in created:
        if not name.startswith("distance_1_"):
            cmd.disable(name)


def write_tsv(
    path: Path, columns: tuple[str, ...], rows: list[dict[str, Any]]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            output = dict(row)
            if "distance_A" in output:
                output["distance_A"] = f"{float(output['distance_A']):.6f}"
            if "min_distance_A" in output:
                output["min_distance_A"] = f"{float(output['min_distance_A']):.9f}"
            for cutoff_field in ("direct_cutoff_A", "near_cutoff_A"):
                if cutoff_field in output:
                    output[cutoff_field] = f"{float(output[cutoff_field]):.3f}"
            if "angle_deg" in output:
                angle = output["angle_deg"]
                output["angle_deg"] = "" if angle is None else f"{float(angle):.6f}"
            writer.writerow({column: output.get(column, "") for column in columns})


def format_metric(value: Any, decimals: int = 6) -> str:
    if isinstance(value, float):
        return f"{value:.{decimals}f}"
    if isinstance(value, list):
        return ", ".join(str(item) for item in value)
    return str(value)


def build_summary(
    spec: dict[str, Any], recomputed: dict[str, Any], imported: dict[str, Any],
    consistency: dict[str, Any], measurements: list[dict[str, Any]]
) -> str:
    values = imported["values"]
    lines = [
        f"# {spec['analysis_id']} parameterized pose analysis — V0.2",
        "",
        "## Identity",
        "",
        f"- Candidate: `{spec['candidate']}`",
        f"- Sequence: `{spec['sequence']}`",
        f"- Group: `{spec['group']}`",
        f"- Conformer: `{spec['conformer']}`",
        f"- Selected Vina MODEL label: `{spec['model']}`",
        "",
        "## Recomputed by V0.2",
        "",
        f"- MODEL block/state count: `{recomputed['state_count']}`",
        "- Available MODEL labels: `" + format_metric(recomputed["model_labels"]) + "`",
        "- Formal contact definition: ligand/receptor heavy atoms, monomer PDB "
        f"chain `{spec['receptor_chain']}`, raw unrounded distance",
        f"- Direct contacts (d <= {float(spec['direct_contact_cutoff_A']):.3f} Å): `"
        + ", ".join(recomputed["direct_target_contact_labels"])
        + "`",
        "- Near contacts ("
        f"{float(spec['direct_contact_cutoff_A']):.3f} < d <= "
        f"{float(spec['near_contact_cutoff_A']):.3f} Å): `"
        + ", ".join(recomputed["near_target_contact_labels"])
        + "`",
        "- Consistency with legacy heavy-atom <= 5 Å pose-metrics contacts: "
        f"`{consistency['status']}`",
    ]
    if consistency["status"] == "EXPECTED_DEFINITION_DIFFERENCE":
        lines.append(
            "- Legacy-only residues explained by near contacts: `"
            + format_metric(consistency["legacy_only_explained_by_near"])
            + "`"
        )
    if consistency["status"] == "WARNING_MISMATCH":
        lines.extend(
            [
                "- Residues not explained by the direct/near versus legacy 5 Å "
                "definition: `"
                + format_metric(consistency["unexplained_residues"])
                + "`",
            ]
        )
    lines.extend(
        [
            "",
            "## Imported from existing pose_metrics.tsv",
            "",
            f"- Vina score: `{format_metric(values['score'], 3)} kcal/mol`",
            f"- Contact fraction: `{format_metric(values['contact_fraction'])}`",
            f"- Receptor clash pairs: `{format_metric(values['receptor_clash_pairs'])}`",
            "- Minimum receptor distance: "
            f"`{format_metric(values['min_receptor_distance'])} Å`",
            f"- Chain-B clash pairs: `{format_metric(values['chain_b_clash_pairs'])}`",
            "- Minimum chain-B distance: "
            f"`{format_metric(values['min_chain_b_distance'])} Å`",
            f"- Basic geometry pass: `{format_metric(values['basic_geometry_pass'])}`",
            "",
            "Chain-B clash and minimum-distance values are imported from the existing "
            "pose_metrics.tsv; they are not recalculated from the AB dimer by V0.2.",
            "The dimer input is validated and hashed for provenance and future extension.",
            "",
            "## Manual-spec measurements recomputed by PyMOL",
            "",
            "| Interaction | Distance (Å) | Angle (°) | Manual interpretation |",
            "|---|---:|---:|---|",
        ]
    )
    for row in measurements:
        angle = "" if row["angle_deg"] is None else f"{row['angle_deg']:.6f}"
        lines.append(
            f"| {row['interaction']} | {row['distance_A']:.6f} | {angle} | "
            f"{row['interpretation']} |"
        )
    if not measurements:
        lines.append("| _No manual interactions supplied_ |  |  |  |")
    lines.extend(
        [
            "",
            "Manual interpretations are copied from the spec and are not inferred by V0.2.",
            "",
        ]
    )
    return "\n".join(lines)


def build_manifest(
    spec_path: Path,
    spec: dict[str, Any],
    data_root: Path,
    inputs: dict[str, Path],
    input_hashes: dict[str, str],
    output_dir: Path,
    recomputed: dict[str, Any],
    imported: dict[str, Any],
    consistency: dict[str, Any],
) -> str:
    expected_hashes = spec.get("input_sha256", {})
    lines = [
        f"run_datetime={datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"analysis_id={spec['analysis_id']}",
        f"selected_model_label={spec['model']}",
        f"data_root={data_root.resolve()}",
        f"spec={spec_path.resolve()}",
        f"spec_sha256={sha256(spec_path.resolve())}",
        f"script={SCRIPT_PATH}",
        f"script_sha256={sha256(SCRIPT_PATH)}",
        f"core_script={CORE_PATH}",
        f"core_script_sha256={sha256(CORE_PATH)}",
        f"pymol_version={cmd.get_version()[0]}",
        f"python_version={platform.python_version()}",
        f"python_executable={Path(sys.executable).resolve()}",
        "provenance.v02_recomputed=MODEL parsing, state count, raw-coordinate heavy-atom direct/near target contacts, manual-spec PyMOL distances and angles",
        "provenance.pose_metrics_imported=Vina score, contact fraction, receptor clashes/minimum distance, chain-B clashes/minimum distance, basic geometry pass",
        "receptor_dimer_role=input validation, hash/provenance, future extension; not used for chain-B recalculation",
        "contact_definition=heavy atoms; receptor monomer PDB and spec chain; selected Vina MODEL PDBQT coordinates; raw unrounded distance; direct d<=4.0 A; near 4.0<d<=5.0 A",
        f"result.legacy_contact_consistency={consistency['status']}",
        "result.legacy_contacts.explained_by_near="
        + ";".join(str(value) for value in consistency["legacy_only_explained_by_near"]),
        "result.legacy_contacts.unexplained="
        + ";".join(str(value) for value in consistency["unexplained_residues"]),
        f"result.v02_recomputed.direct_target_contacts={';'.join(str(v) for v in recomputed['direct_target_contacts'])}",
        f"result.v02_recomputed.near_target_contacts={';'.join(str(v) for v in recomputed['near_target_contacts'])}",
        f"result.pose_metrics_imported.legacy_5A_target_contacts={format_metric(imported['legacy_5A_target_contacts'])}",
    ]
    for field, path in inputs.items():
        expected = expected_hashes.get(field)
        lines.append(f"input.{field}.spec_path={spec[field]}")
        lines.append(f"input.{field}.resolved_path={path}")
        lines.append(f"input.{field}.sha256={input_hashes[field]}")
        lines.append(f"input.{field}.expected_sha256={expected or 'NOT_PROVIDED'}")
        lines.append(
            f"input.{field}.hash_verification={'MATCH' if expected else 'NOT_PROVIDED'}"
        )
    for name in OUTPUT_NAMES:
        lines.append(f"output={output_dir / name}")
    return "\n".join(lines) + "\n"


def configure_scene(spec: dict[str, Any]) -> None:
    cmd.disable("ligand_all_states")
    cmd.hide("everything", "all")
    cmd.show("cartoon", "receptor_monomer")
    cmd.color("gray80", "receptor_monomer")
    cmd.show("sticks", "target_site")
    cmd.color("orange", "target_site")
    cmd.show("sticks", "ligand_pose")
    cmd.color("magenta", "ligand_pose")
    cmd.show("sticks", "direct_target_contacts")
    cmd.color("orange", "direct_target_contacts")
    cmd.show("sticks", "near_target_contacts")
    cmd.color("yellow", "near_target_contacts")
    cmd.label(
        "(direct_target_contacts or near_target_contacts) and name CA",
        '"%s%s" % (resn, resi)',
    )
    cmd.set("label_color", "black")
    cmd.set("label_size", 20)
    cmd.set("label_font_id", 7)
    cmd.set("dash_color", "tv_blue")
    cmd.set("dash_width", 4.0)
    cmd.set("dash_gap", 0.25)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    focus = "ligand_pose or direct_target_contacts or near_target_contacts"
    cmd.orient(focus)
    cmd.zoom(focus, buffer=7.0, complete=1)
    create_measurement_objects(spec)


def analyze_pose(spec_path: Path, data_root: Path, output_dir: Path) -> Path:
    resolved_spec_path = spec_path.resolve()
    spec = load_json(resolved_spec_path)
    validate_spec(spec)
    safe_output_dir = require_safe_output_directory(output_dir)
    inputs, input_hashes = resolve_inputs(spec, data_root)
    model_blocks = parse_vina_model_blocks(inputs["docking_output_path"])
    model_label = int(spec["model"])
    selected_model_block = require_model_block(model_blocks, model_label)
    metrics_header, pose_metrics = read_model_metrics(
        inputs["pose_metrics_path"], model_label
    )
    ligand_atoms = parse_pdbqt_atoms(selected_model_block)
    receptor_atoms = parse_receptor_pdb_atoms(
        inputs["receptor_monomer_path"], spec["receptor_chain"]
    )
    contact_rows = compute_target_residue_contacts(
        ligand_atoms,
        receptor_atoms,
        spec["target_residues"],
        float(spec["direct_contact_cutoff_A"]),
        float(spec["near_contact_cutoff_A"]),
    )
    direct_rows = [row for row in contact_rows if row["classification"] == "direct"]
    near_rows = [row for row in contact_rows if row["classification"] == "near"]
    direct_numbers = [row["target_residue_number"] for row in direct_rows]
    near_numbers = [row["target_residue_number"] for row in near_rows]
    direct_names = [row["target_residue"] for row in direct_rows]
    near_names = [row["target_residue"] for row in near_rows]
    consistency = compare_legacy_contacts(
        direct_numbers, near_numbers, pose_metrics["target_contacts"]
    )
    if consistency["status"] == "WARNING_MISMATCH":
        print(
            "[WARNING] legacy 5 Å target contacts cannot be reconciled with "
            "V0.2 direct/near heavy-atom contacts; unexplained residues="
            f"{consistency['unexplained_residues']}",
            file=sys.stderr,
        )

    pymol.finish_launching(["pymol", "-cq"])
    try:
        cmd.load(str(inputs["receptor_monomer_path"]), "receptor_monomer")
        label_to_state: dict[int, int] = {}
        for state_index, (label, block) in enumerate(model_blocks.items(), start=1):
            cmd.read_pdbstr(block, "ligand_all_states", state=state_index)
            label_to_state[label] = state_index
        state_count = cmd.count_states("ligand_all_states")
        if state_count != len(model_blocks):
            raise ValueError(
                f"PyMOL state count {state_count} differs from parsed MODEL block count "
                f"{len(model_blocks)}"
            )
        cmd.create(
            "ligand_pose",
            "ligand_all_states",
            label_to_state[model_label],
            1,
        )

        target_expression = "+".join(str(value) for value in spec["target_residues"])
        target_site_atom_count = cmd.select(
            "target_site",
            "receptor_monomer and "
            f"chain {spec['receptor_chain']} and resi {target_expression}",
        )
        if target_site_atom_count == 0:
            raise ValueError("target_site selection is empty")
        target_labels = collect_ca_residues("target_site")
        missing_targets = sorted(set(spec["target_residues"]) - set(target_labels))
        if missing_targets:
            raise ValueError(
                f"Target residues have no CA atom in receptor chain "
                f"{spec['receptor_chain']}: {missing_targets}"
            )

        def select_residues(name: str, residue_numbers: list[int]) -> None:
            if residue_numbers:
                expression = "+".join(str(value) for value in residue_numbers)
                cmd.select(
                    name,
                    "receptor_monomer and "
                    f"chain {spec['receptor_chain']} and resi {expression}",
                )
            else:
                cmd.select(name, "none")

        select_residues("direct_target_contacts", direct_numbers)
        select_residues("near_target_contacts", near_numbers)
        measurements = calculate_measurements(spec, target_labels)

        recomputed = {
            "provenance": "v02_recomputed",
            "field_provenance": {
                "model_mapping_and_state_count": (
                    "in_memory_MODEL_ENDMDL_split_plus_PyMOL_state_construction"
                ),
                "direct_and_near_target_contacts": (
                    "raw_PDB_and_PDBQT_coordinate_heavy_atom_recomputation"
                ),
                "manual_measurements": "PyMOL_get_distance_and_get_angle",
            },
            "model_labels": list(model_blocks),
            "model_label_to_pymol_state": {
                str(label): state for label, state in label_to_state.items()
            },
            "selected_model_label": model_label,
            "selected_pymol_state": label_to_state[model_label],
            "state_count": state_count,
            "parsed_model_block_count": len(model_blocks),
            "pymol_state_count": state_count,
            "target_site_atom_count": target_site_atom_count,
            "contact_definition": {
                "atom_mode": "heavy_atoms_only",
                "receptor_representation": "monomer_PDB_spec_chain",
                "ligand_representation": "selected_Vina_MODEL_PDBQT_coordinates",
                "distance_comparison": "raw_unrounded_float",
                "direct_cutoff_A": float(spec["direct_contact_cutoff_A"]),
                "near_cutoff_A": float(spec["near_contact_cutoff_A"]),
            },
            "direct_target_contacts": direct_numbers,
            "direct_target_contact_labels": direct_names,
            "near_target_contacts": near_numbers,
            "near_target_contact_labels": near_names,
            "target_residue_contacts": contact_rows,
            "measurements": measurements,
        }
        imported_fields = [
            "score",
            "contact_fraction",
            "contacted_peptide_residues",
            "receptor_contact_residues",
            "target_contacts",
            "target_contact_count",
            "receptor_clash_pairs",
            "min_receptor_distance",
            "chain_b_clash_pairs",
            "min_chain_b_distance",
            "basic_geometry_pass",
        ]
        imported = {
            "provenance": "existing_pose_metrics_tsv",
            "source_path": str(inputs["pose_metrics_path"]),
            "header": metrics_header,
            "values": pose_metrics,
            "legacy_5A_target_contacts": pose_metrics["target_contacts"],
            "legacy_target_contact_definition": "heavy atoms, distance <= 5.0 A",
            "field_provenance": {
                field: "existing_pose_metrics_tsv" for field in imported_fields
            },
            "chain_b_note": (
                "Imported from pose_metrics.tsv; receptor dimer was not used for "
                "chain-B recalculation in V0.2."
            ),
        }
        input_details = {
            field: {
                "spec_path": spec[field],
                "resolved_path": str(path),
                "sha256": input_hashes[field],
                "expected_sha256": spec.get("input_sha256", {}).get(field),
                "hash_verification": (
                    "MATCH"
                    if field in spec.get("input_sha256", {})
                    else "NOT_PROVIDED"
                ),
                "role": (
                    "input_validation_hash_provenance_future_extension"
                    if field == "receptor_dimer_path"
                    else "analysis_input"
                ),
            }
            for field, path in inputs.items()
        }
        result_metrics = {
            "schema_version": "0.2",
            "analysis_id": spec["analysis_id"],
            "identity": {
                field: spec[field]
                for field in ("candidate", "sequence", "group", "conformer")
            },
            "input_files": input_details,
            "v02_recomputed": recomputed,
            "pose_metrics_imported": imported,
            "consistency_checks": {"legacy_5A_target_contacts": consistency},
        }

        staging_dir = safe_output_dir.with_name(
            f".{safe_output_dir.name}.incomplete_{hashlib.sha256(str(safe_output_dir).encode()).hexdigest()[:12]}"
        )
        if staging_dir.exists():
            raise FileExistsError(
                f"Staging directory already exists; refusing to overwrite: {staging_dir}"
            )
        staging_dir.mkdir(parents=True, exist_ok=False)
        write_tsv(
            staging_dir / "contacts.tsv",
            (
                "target_residue_number",
                "target_residue",
                "classification",
                "min_distance_A",
                "ligand_atom_id",
                "ligand_atom_name",
                "ligand_atom_type",
                "receptor_atom",
                "direct_cutoff_A",
                "near_cutoff_A",
                "provenance",
            ),
            contact_rows,
        )
        write_tsv(staging_dir / "measurements.tsv", MEASUREMENT_COLUMNS, measurements)
        (staging_dir / "metrics.json").write_text(
            json.dumps(result_metrics, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        (staging_dir / "summary.md").write_text(
            build_summary(spec, recomputed, imported, consistency, measurements),
            encoding="utf-8",
            newline="\n",
        )
        configure_scene(spec)
        cmd.save(str(staging_dir / "interface.pse"))
        cmd.png(
            str(staging_dir / "interface.png"),
            width=2400,
            height=1800,
            dpi=300,
            ray=1,
            quiet=0,
        )
        (staging_dir / "run_manifest.txt").write_text(
            build_manifest(
                resolved_spec_path,
                spec,
                data_root,
                inputs,
                input_hashes,
                safe_output_dir,
                recomputed,
                imported,
                consistency,
            ),
            encoding="utf-8",
            newline="\n",
        )
        staging_dir.replace(safe_output_dir)
    finally:
        cmd.quit()
    return safe_output_dir


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Analyze one Vina MODEL pose with headless PyMOL."
    )
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        output_dir = analyze_pose(args.spec, args.data_root, args.output_dir)
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print(f"Analysis completed: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
