"""Zero-docking V0.5-A batch planner and preflight."""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Any

from candidate_manifest_v05 import (
    load_json_contract,
    select_candidates,
    validate_candidate_manifest,
    validate_ligand_inventory,
    validate_screening_protocol,
)
from environment_preflight_v05 import audit_runtime_environment, parse_requirements_lock
from receptor_registry import resolve_receptor_bundle, validate_spec_against_bundle
from standardized_vina_inputs import (
    audit_source_sdf,
    comparison_protocol_id,
    read_single_sdf,
    sha256,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCHEMA_ROOT = PROJECT_ROOT / "specs"
EXPECTED_CHEMISTRY_RULES = {
    "N_terminus": "protonated",
    "C_terminus": "deprotonated",
    "ASP": "deprotonated",
    "GLU": "deprotonated",
    "LYS": "protonated",
    "ARG": "protonated",
    "HIS": "neutral_HIE_like_NE2_protonated_ND1_unprotonated",
    "CYS": "neutral_free_thiol_unless_explicitly_approved_otherwise",
    "TYR": "neutral",
    "amino_acid_stereochemistry": "standard_L",
    "artificial_disulfide": False,
}


def failure(
    phase: str, error_code: str, message: str, candidate: str | None = None
) -> dict[str, Any]:
    return {
        "candidate": candidate,
        "phase": phase,
        "error_code": error_code,
        "message": message,
    }


def _resolve_relative(base: Path, raw: str, label: str) -> Path:
    path = Path(raw)
    if path.is_absolute():
        raise ValueError(f"{label} must be relative")
    base_resolved = base.resolve()
    resolved = (base_resolved / path).resolve()
    if not resolved.is_relative_to(base_resolved):
        raise ValueError(f"{label} resolves outside its declared base")
    return resolved


def _verify_file(path: Path, expected_hash: str, label: str) -> str:
    if not path.is_file():
        raise FileNotFoundError(f"{label} is missing: {path}")
    actual = sha256(path)
    if actual.lower() != expected_hash.lower():
        raise ValueError(
            f"SHA-256 mismatch for {label}: expected {expected_hash.lower()}, actual {actual}"
        )
    return actual


def _tracked_contract_audit(
    path: Path,
    schema_path: Path,
    *,
    expected_hash: str | None = None,
) -> dict[str, Any]:
    actual = sha256(path)
    if expected_hash is not None and actual.lower() != expected_hash.lower():
        raise ValueError(
            f"SHA-256 mismatch for {path.name}: expected {expected_hash.lower()}, actual {actual}"
        )
    return {
        "path": str(path),
        "sha256": actual,
        "expected_sha256": expected_hash.lower() if expected_hash is not None else None,
        "schema_path": str(schema_path),
        "schema_sha256": sha256(schema_path),
        "schema_status": "PASS",
        "status": "PASS",
    }


def inspect_vina_identity(path: Path, expected_hash: str, expected_version: str) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Vina executable is missing: {path}")
    actual = sha256(path)
    if actual.lower() != expected_hash.lower():
        raise ValueError(
            f"Vina executable SHA-256 mismatch: expected {expected_hash.lower()}, actual {actual}"
        )
    return {
        "path": str(path.resolve()),
        "sha256": actual,
        "version": expected_version,
        "version_verification": "verified_by_locked_executable_sha256_without_process_execution",
        "process_executed": False,
    }


def _git_directory(repository: Path) -> Path:
    marker = repository / ".git"
    if marker.is_dir():
        return marker
    if marker.is_file():
        content = marker.read_text(encoding="utf-8").strip()
        if not content.startswith("gitdir: "):
            raise ValueError(f"Malformed .git file: {marker}")
        value = Path(content[8:])
        return value.resolve() if value.is_absolute() else (repository / value).resolve()
    raise FileNotFoundError(f"Git metadata is missing from {repository}")


def inspect_git_release(repository: Path) -> dict[str, Any]:
    git_dir = _git_directory(repository)
    head_text = (git_dir / "HEAD").read_text(encoding="ascii").strip()
    branch: str | None = None
    if head_text.startswith("ref: "):
        ref = head_text[5:]
        branch = ref.removeprefix("refs/heads/")
        ref_file = git_dir / ref
        if ref_file.is_file():
            commit = ref_file.read_text(encoding="ascii").strip()
        else:
            packed = (git_dir / "packed-refs").read_text(encoding="ascii")
            matches = [line.split()[0] for line in packed.splitlines() if line.endswith(f" {ref}")]
            if len(matches) != 1:
                raise ValueError(f"Cannot resolve Git HEAD ref {ref}")
            commit = matches[0]
    else:
        commit = head_text
    if re.fullmatch(r"[0-9a-f]{40}", commit) is None:
        raise ValueError(f"Invalid Git HEAD commit: {commit!r}")
    return {
        "repository": str(repository.resolve()),
        "branch": branch,
        "head_commit": commit,
        "scientific_semantics_baseline_tag": "v0.4",
    }


def expected_formal_charge(sequence: str) -> int:
    return sequence.count("R") + sequence.count("K") - sequence.count("D") - sequence.count("E")


def _hydrogen_neighbors(atom: Any) -> int:
    return sum(neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors())


def audit_frozen_chemistry(path: Path, candidate: str, sequence: str) -> dict[str, Any]:
    from rdkit import Chem

    charge = expected_formal_charge(sequence)
    base_audit = audit_source_sdf(path, candidate, sequence, charge)
    molecule = read_single_sdf(path)
    template = Chem.MolFromSequence(sequence)
    if template is None:
        raise ValueError(f"Cannot build peptide template for {sequence}")
    source_heavy = [atom for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
    template_heavy = [atom for atom in template.GetAtoms() if atom.GetAtomicNum() > 1]
    sites: dict[tuple[int, str, str], Any] = {}
    for source_atom, template_atom in zip(source_heavy, template_heavy):
        info = template_atom.GetPDBResidueInfo()
        if info is not None:
            sites[(info.GetResidueNumber(), info.GetResidueName().strip(), info.GetName().strip())] = source_atom
    expected_sites = {
        "D": ("ASP", "OD2", -1, 0),
        "E": ("GLU", "OE2", -1, 0),
        "K": ("LYS", "NZ", 1, 3),
        "R": ("ARG", "NH1", 1, 2),
        "C": ("CYS", "SG", 0, 1),
        "Y": ("TYR", "OH", 0, 1),
    }
    checks: dict[str, Any] = {}
    for residue_number, code in enumerate(sequence, 1):
        if code in expected_sites:
            residue_name, atom_name, formal_charge, hydrogen_count = expected_sites[code]
            atom = sites[(residue_number, residue_name, atom_name)]
            actual = (atom.GetFormalCharge(), _hydrogen_neighbors(atom))
            if actual != (formal_charge, hydrogen_count):
                raise ValueError(
                    f"Frozen chemistry mismatch at {candidate} {residue_name}{residue_number} {atom_name}: "
                    f"expected charge/H {(formal_charge, hydrogen_count)}, actual {actual}"
                )
            checks[f"{residue_name}{residue_number}_{atom_name}"] = {
                "formal_charge": actual[0],
                "hydrogen_neighbors": actual[1],
            }
        elif code == "H":
            nd1 = sites[(residue_number, "HIS", "ND1")]
            ne2 = sites[(residue_number, "HIS", "NE2")]
            actual = (
                nd1.GetFormalCharge(),
                _hydrogen_neighbors(nd1),
                ne2.GetFormalCharge(),
                _hydrogen_neighbors(ne2),
            )
            if actual != (0, 0, 0, 1):
                raise ValueError(
                    f"Frozen HIE-like histidine mismatch at {candidate} HIS{residue_number}: {actual}"
                )
            checks[f"HIS{residue_number}"] = {
                "ND1_formal_charge": actual[0],
                "ND1_hydrogen_neighbors": actual[1],
                "NE2_formal_charge": actual[2],
                "NE2_hydrogen_neighbors": actual[3],
            }
    return {**base_audit, "expected_formal_charge": charge, "side_chain_checks": checks}


def _read_audit(path: Path, expected_hash: str, candidate: str, conformer: str, sdf_hash: str) -> dict[str, Any]:
    _verify_file(path, expected_hash, f"{candidate}_{conformer} audit")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Audit root must be an object: {path}")
    run_id = f"{candidate}_{conformer}"
    if value.get("candidate") != candidate:
        raise ValueError(f"Audit candidate mismatch for {run_id}")
    if "run_id" in value and value["run_id"] != run_id:
        raise ValueError(f"Audit run_id mismatch for {run_id}")
    if value.get("status") not in (None, "PASS"):
        raise ValueError(f"Audit status is not PASS for {run_id}")
    source = value.get("source_sdf")
    if isinstance(source, dict) and source.get("sha256") != sdf_hash:
        raise ValueError(f"Audit source SDF hash mismatch for {run_id}")
    if isinstance(value.get("conformers"), list):
        matches = [row for row in value["conformers"] if row.get("conformer") == conformer]
        if len(matches) != 1 or matches[0].get("standardized_sha256") != sdf_hash:
            raise ValueError(f"Conversion audit conformer/hash mismatch for {run_id}")
    return value


def _formal_chemistry_spec(sequence: str) -> dict[str, Any]:
    return {
        "formal_charge": expected_formal_charge(sequence),
        "chemistry_contract_id": "v04_formal_standardized_peptide_chemistry",
        "n_terminus": "protonated",
        "c_terminus": "deprotonated",
        "amino_acid_stereochemistry": "standard_L",
        "artificial_disulfide": False,
    }


def project_candidate_paths(
    candidate_id: str, conformer_names: list[str]
) -> dict[str, Any]:
    candidate_root = f"candidates/{candidate_id}"
    return {
        "candidate_root": candidate_root,
        "runtime_spec": f"{candidate_root}/specs/{candidate_id}_standardized_docking_v05.json",
        "conformers": {
            conformer: {
                "ligand_pdbqt": f"{candidate_root}/ligands/{candidate_id}_{conformer}.pdbqt",
                "config": f"{candidate_root}/configs/{candidate_id}_{conformer}.txt",
                "vina_output": f"{candidate_root}/outputs/{candidate_id}_{conformer}_out.pdbqt",
                "vina_log": f"{candidate_root}/logs/{candidate_id}_{conformer}.log",
                "pose_metrics": f"{candidate_root}/analysis/{candidate_id}_{conformer}_pose_metrics.tsv",
            }
            for conformer in conformer_names
        },
    }


def _protocol_basis(
    candidate: dict[str, Any], protocol: dict[str, Any], bundle: dict[str, Any]
) -> dict[str, Any]:
    return {
        "ligand_preparation": protocol["ligand_preparation"],
        "group": candidate["group"],
        "template_id": bundle["template_id"],
        "receptor_pdbqt_sha256": bundle["sha256"]["receptor_pdbqt"],
        "box_center_A": tuple(bundle["vina_box"]["center"]),
        "box_size_A": tuple(bundle["vina_box"]["size"]),
        "docking_protocol": protocol["docking_protocol"],
    }


def _runtime_spec(
    candidate: dict[str, Any],
    slots: list[dict[str, Any]],
    protocol: dict[str, Any],
    bundle: dict[str, Any],
) -> tuple[dict[str, Any], str]:
    basis = _protocol_basis(candidate, protocol, bundle)
    protocol_id = comparison_protocol_id(basis)
    registry_contract = protocol["contracts"]["receptor_registry"]
    spec = {
        "schema_version": "0.4",
        "analysis_id": f"{candidate['candidate_id']}_standardized_docking_v05",
        "candidate": candidate["candidate_id"],
        "sequence": candidate["sequence"],
        "peptide_length": candidate["peptide_length"],
        "group": candidate["group"],
        "template_id": bundle["template_id"],
        "receptor_id": f"{bundle['template_id']}_TfR1_A",
        "receptor_registry_path": registry_contract["path"],
        "receptor_registry_sha256": registry_contract["sha256"],
        "screening_protocol_status": protocol["screening_protocol_status"],
        "formal_chemistry": _formal_chemistry_spec(candidate["sequence"]),
        "ligand_preparation": protocol["ligand_preparation"],
        "docking_protocol": protocol["docking_protocol"],
        "conformers": [
            {
                "conformer_name": slot["conformer_name"],
                "source_sdf_path": slot["source_sdf_path"],
                "source_sdf_sha256": slot["source_sdf_sha256"],
            }
            for slot in slots
        ],
    }
    validate_spec_against_bundle(spec, bundle)
    return spec, protocol_id


def _error_code(exc: Exception, default: str) -> str:
    message = str(exc).lower()
    if "hash" in message or "sha-256" in message:
        return f"{default}_HASH_MISMATCH"
    if "missing" in message or isinstance(exc, FileNotFoundError):
        return f"{default}_MISSING"
    return f"{default}_INVALID"


def run_dry_run(
    manifest_path: Path,
    protocol_path: Path,
    data_root: Path,
    output_root: Path,
    vina_path: Path,
    *,
    groups: list[str] | None = None,
    candidate_ids: list[str] | None = None,
) -> dict[str, Any]:
    manifest_path = manifest_path.resolve()
    protocol_path = protocol_path.resolve()
    data_root = data_root.resolve()
    output_root = output_root.resolve()
    failures: list[dict[str, Any]] = []
    global_checks: dict[str, Any] = {}

    try:
        manifest = load_json_contract(
            manifest_path, SCHEMA_ROOT / "top40_candidates_v05.schema.json"
        )
        candidates = validate_candidate_manifest(manifest)
    except Exception as exc:
        return {
            "outcome": "DRY_RUN_FAIL_GLOBAL",
            "selected_candidate_count": 0,
            "selected_candidate_ids": [],
            "plan": [],
            "global_checks": global_checks,
            "failures": [failure("global_preflight", "CANDIDATE_MANIFEST_INVALID", str(exc))],
            "vina_docking_process_count": 0,
        }
    try:
        protocol = load_json_contract(
            protocol_path, SCHEMA_ROOT / "screening_protocol_v05.schema.json"
        )
        protocol = validate_screening_protocol(protocol)
        selected = select_candidates(candidates, groups=groups, candidate_ids=candidate_ids)
    except Exception as exc:
        return {
            "outcome": "DRY_RUN_FAIL_GLOBAL",
            "selected_candidate_count": 0,
            "selected_candidate_ids": [],
            "plan": [],
            "global_checks": global_checks,
            "failures": [failure("global_preflight", "PROTOCOL_OR_SELECTOR_INVALID", str(exc))],
            "vina_docking_process_count": 0,
        }
    selected_ids = [row["candidate_id"] for row in selected]
    global_checks["tracked_contracts"] = {
        "screening_protocol": _tracked_contract_audit(
            protocol_path, SCHEMA_ROOT / "screening_protocol_v05.schema.json"
        )
    }

    try:
        manifest_expected = protocol["contracts"]["candidate_manifest"]["sha256"]
        global_checks["candidate_manifest"] = {
            "path": str(manifest_path),
            "sha256": _verify_file(manifest_path, manifest_expected, "candidate manifest"),
            "status": "PASS",
        }
        global_checks["tracked_contracts"]["candidate_manifest"] = _tracked_contract_audit(
            manifest_path,
            SCHEMA_ROOT / "top40_candidates_v05.schema.json",
            expected_hash=manifest_expected,
        )
    except Exception as exc:
        failures.append(failure("global_preflight", "CANDIDATE_MANIFEST_HASH_MISMATCH", str(exc)))

    inventory: dict[str, Any] | None = None
    slots_by_key: dict[tuple[str, str], dict[str, Any]] = {}
    try:
        inventory_contract = protocol["contracts"]["ligand_inventory"]
        inventory_path = _resolve_relative(
            protocol_path.parent, inventory_contract["path"], "ligand inventory path"
        )
        _verify_file(inventory_path, inventory_contract["sha256"], "ligand inventory")
        inventory = load_json_contract(
            inventory_path, SCHEMA_ROOT / "top40_ligand_inventory_v05.schema.json"
        )
        if inventory.get("candidate_manifest_sha256") != sha256(manifest_path):
            raise ValueError("Ligand inventory candidate manifest hash mismatch")
        slots_by_key = validate_ligand_inventory(inventory, candidates, protocol)
        global_checks["ligand_inventory"] = {
            "path": str(inventory_path),
            "sha256": sha256(inventory_path),
            "slot_count": len(slots_by_key),
            "status": "PASS",
        }
        global_checks["tracked_contracts"]["ligand_inventory"] = _tracked_contract_audit(
            inventory_path,
            SCHEMA_ROOT / "top40_ligand_inventory_v05.schema.json",
            expected_hash=inventory_contract["sha256"],
        )
    except Exception as exc:
        failures.append(failure("global_preflight", _error_code(exc, "LIGAND_INVENTORY"), str(exc)))

    registry_path: Path | None = None
    bundles: dict[str, dict[str, Any]] = {}
    try:
        contract = protocol["contracts"]["receptor_registry"]
        registry_path = _resolve_relative(data_root, contract["path"], "receptor registry path")
        _verify_file(registry_path, contract["sha256"], "receptor registry")
        for group in "ABCD":
            bundles[group] = resolve_receptor_bundle(registry_path, data_root, group)
        global_checks["receptor_registry"] = {
            "path": str(registry_path),
            "sha256": sha256(registry_path),
            "resolved_groups": list(bundles),
            "status": "PASS",
        }
    except Exception as exc:
        failures.append(failure("global_preflight", _error_code(exc, "RECEPTOR_REGISTRY"), str(exc)))

    chemistry: dict[str, Any] | None = None
    try:
        contract = protocol["contracts"]["chemistry_contract"]
        chemistry_path = _resolve_relative(data_root, contract["path"], "chemistry contract path")
        _verify_file(chemistry_path, contract["sha256"], "chemistry contract")
        chemistry = json.loads(chemistry_path.read_text(encoding="utf-8"))
        if chemistry.get("contract_id") != contract["contract_id"]:
            raise ValueError("Chemistry contract ID mismatch")
        if chemistry.get("status") != "FROZEN" or chemistry.get("rules") != EXPECTED_CHEMISTRY_RULES:
            raise ValueError("Chemistry contract rules/status differ from the frozen V0.4 contract")
        if chemistry.get("top40_coverage_audit", {}).get("status") != "PASS":
            raise ValueError("Chemistry contract Top40 coverage audit is not PASS")
        global_checks["chemistry_contract"] = {
            "path": str(chemistry_path),
            "sha256": sha256(chemistry_path),
            "status": "PASS",
        }
    except Exception as exc:
        failures.append(failure("global_preflight", _error_code(exc, "CHEMISTRY_CONTRACT"), str(exc)))

    try:
        contract = protocol["contracts"]["dependency_lock"]
        if contract.get("path_base") != "project_root":
            raise ValueError("dependency lock path_base must be project_root")
        lock_path = _resolve_relative(PROJECT_ROOT, contract["path"], "dependency lock path")
        _verify_file(lock_path, contract["sha256"], "dependency lock")
        lock = parse_requirements_lock(lock_path.read_text(encoding="utf-8"))
        global_checks["runtime_environment"] = audit_runtime_environment(lock)
        global_checks["runtime_environment"]["lock_path"] = str(lock_path)
        global_checks["runtime_environment"]["lock_sha256"] = sha256(lock_path)
    except Exception as exc:
        failures.append(failure("global_preflight", _error_code(exc, "RUNTIME_ENVIRONMENT"), str(exc)))

    try:
        contract = protocol["contracts"]["analysis_schema"]
        if contract.get("path_base") != "project_root":
            raise ValueError("analysis schema path_base must be project_root")
        analysis_schema_path = _resolve_relative(
            PROJECT_ROOT, contract["path"], "analysis schema path"
        )
        global_checks["analysis_schema"] = {
            "path": str(analysis_schema_path),
            "sha256": _verify_file(
                analysis_schema_path, contract["sha256"], "analysis schema"
            ),
            "status": "PASS",
        }
    except Exception as exc:
        failures.append(
            failure("global_preflight", _error_code(exc, "ANALYSIS_SCHEMA"), str(exc))
        )

    try:
        global_checks["vina"] = inspect_vina_identity(
            vina_path.resolve(),
            protocol["vina_executable"]["sha256"],
            protocol["vina_executable"]["version"],
        )
    except FileNotFoundError as exc:
        failures.append(failure("global_preflight", "VINA_EXECUTABLE_MISSING", str(exc)))
    except Exception as exc:
        failures.append(failure("global_preflight", "VINA_EXECUTABLE_HASH_MISMATCH", str(exc)))

    try:
        if output_root.exists():
            raise FileExistsError(f"Output root already exists: {output_root}")
        parent = output_root.parent
        if not parent.is_dir() or not os.access(parent, os.W_OK):
            raise PermissionError(f"Output parent is not an existing writable directory: {parent}")
        global_checks["output_root"] = {
            "path": str(output_root),
            "collision": False,
            "parent_write_access": True,
            "status": "PASS",
        }
    except FileExistsError as exc:
        failures.append(failure("global_preflight", "OUTPUT_PATH_COLLISION", str(exc)))
    except Exception as exc:
        failures.append(failure("global_preflight", "OUTPUT_PATH_INVALID", str(exc)))

    try:
        global_checks["git"] = inspect_git_release(PROJECT_ROOT)
    except Exception as exc:
        failures.append(failure("global_preflight", "GIT_RELEASE_INSPECTION_FAILED", str(exc)))

    if failures:
        return {
            "outcome": "DRY_RUN_FAIL_GLOBAL",
            "selected_candidate_count": len(selected),
            "selected_candidate_ids": selected_ids,
            "plan": [],
            "global_checks": global_checks,
            "failures": failures,
            "vina_docking_process_count": 0,
        }

    assert inventory is not None and registry_path is not None and chemistry is not None
    plan: list[dict[str, Any]] = []
    projected: set[str] = set()
    candidate_failures: list[dict[str, Any]] = []
    conformer_names = protocol["formal_conformers"]["required_conformer_names"]
    for candidate in selected:
        candidate_id = candidate["candidate_id"]
        slots = [slots_by_key[(candidate_id, name)] for name in conformer_names]
        bundle = bundles[candidate["group"]]
        projected_paths = project_candidate_paths(candidate_id, conformer_names)
        flattened_paths = [
            projected_paths["candidate_root"],
            projected_paths["runtime_spec"],
            *[
                path
                for conformer_paths in projected_paths["conformers"].values()
                for path in conformer_paths.values()
            ],
        ]
        duplicates = sorted(path for path in flattened_paths if path in projected)
        if duplicates:
            candidate_failures.append(
                failure(
                    "candidate_preflight",
                    "PROJECTED_OUTPUT_COLLISION",
                    f"Duplicate projected paths: {duplicates}",
                    candidate_id,
                )
            )
        projected.update(flattened_paths)
        comparison_id = comparison_protocol_id(_protocol_basis(candidate, protocol, bundle))
        plan_entry: dict[str, Any] = {
            "candidate_id": candidate_id,
            "group": candidate["group"],
            "sequence": candidate["sequence"],
            "peptide_length": candidate["peptide_length"],
            "template_id": bundle["template_id"],
            "receptor_id": f"{bundle['template_id']}_TfR1_A",
            "receptor_paths": bundle["source_paths"],
            "target_residues": bundle["target_residues"],
            "vina_box": bundle["vina_box"],
            "comparison_protocol_id": comparison_id,
            "conformers": [],
            "runtime_spec": None,
            "projected_paths": projected_paths,
            "projected_output_path": projected_paths["candidate_root"],
            "preflight_status": "FAIL",
        }
        validated_slots: list[dict[str, Any]] = []
        ready_slots: list[dict[str, Any]] = []
        for slot in slots:
            conformer = slot["conformer_name"]
            status = slot["readiness_status"]
            if status == "not_prepared":
                candidate_failures.append(
                    failure(
                        "candidate_preflight",
                        "LIGAND_NOT_PREPARED",
                        f"Formal ligand slot is not prepared: {candidate_id} {conformer}",
                        candidate_id,
                    )
                )
                continue
            try:
                sdf = _resolve_relative(data_root, slot["source_sdf_path"], "source SDF path")
                audit = _resolve_relative(data_root, slot["audit_path"], "ligand audit path")
                _verify_file(sdf, slot["source_sdf_sha256"], f"{candidate_id}_{conformer} SDF")
                _read_audit(
                    audit,
                    slot["audit_sha256"],
                    candidate_id,
                    conformer,
                    slot["source_sdf_sha256"],
                )
                chemistry_audit = audit_frozen_chemistry(
                    sdf, candidate_id, candidate["sequence"]
                )
                validated = {
                    **slot,
                    "resolved_source_sdf": str(sdf),
                    "resolved_audit": str(audit),
                    "artifact_validation_status": "PASS",
                    "chemistry_status": "PASS",
                    "formal_charge": chemistry_audit["formal_charge"],
                }
                validated_slots.append(validated)
                if status == "approved_for_formal_screening":
                    ready_slots.append(validated)
                else:
                    candidate_failures.append(
                        failure(
                            "candidate_preflight",
                            "LIGAND_NOT_APPROVED",
                            f"Formal ligand slot is not approved for screening: {candidate_id} {conformer} ({status})",
                            candidate_id,
                        )
                    )
            except Exception as exc:
                validated_slots.append(
                    {
                        **slot,
                        "artifact_validation_status": "FAIL",
                        "artifact_validation_error": str(exc),
                    }
                )
                candidate_failures.append(
                    failure(
                        "candidate_preflight",
                        _error_code(exc, "LIGAND"),
                        str(exc),
                        candidate_id,
                    )
                )
        validated_by_name = {
            slot["conformer_name"]: slot for slot in validated_slots
        }
        plan_entry["conformers"] = [
            {
                **slot,
                **validated_by_name.get(
                    slot["conformer_name"],
                    {"artifact_validation_status": "NOT_AVAILABLE"},
                ),
                "projected_paths": projected_paths["conformers"][slot["conformer_name"]],
            }
            for slot in slots
        ]
        if len(ready_slots) != len(conformer_names) or duplicates:
            plan.append(plan_entry)
            continue
        try:
            runtime_spec, comparison_id = _runtime_spec(candidate, ready_slots, protocol, bundle)
            plan_entry["comparison_protocol_id"] = comparison_id
            plan_entry["runtime_spec"] = runtime_spec
            plan_entry["preflight_status"] = "PASS"
        except Exception as exc:
            candidate_failures.append(
                failure(
                    "candidate_preflight",
                    "RUNTIME_SPEC_OR_REGISTRY_MISMATCH",
                    str(exc),
                    candidate_id,
                )
            )
        plan.append(plan_entry)

    outcome = "DRY_RUN_PASS" if not candidate_failures else "DRY_RUN_FAIL_CANDIDATE"
    return {
        "outcome": outcome,
        "selected_candidate_count": len(selected),
        "selected_candidate_ids": selected_ids,
        "planned_ready_candidate_count": sum(
            row["preflight_status"] == "PASS" for row in plan
        ),
        "plan": plan,
        "global_checks": global_checks,
        "failures": candidate_failures,
        "vina_docking_process_count": 0,
    }
