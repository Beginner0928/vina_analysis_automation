"""Load and validate the V0.5 Top40 identity, ligand, and protocol contracts."""

from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path
from typing import Any


STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
FORMAL_GROUPS = ("A", "B", "C", "D")
FORMAL_TOP40_IDS = tuple(
    f"{group}{index}" for group in FORMAL_GROUPS for index in range(1, 11)
)
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
READINESS_STATUSES = frozenset(
    {
        "approved_for_formal_screening",
        "approved_for_exploratory_screening",
        "conversion_audit_only",
        "not_prepared",
        "conformer_generation_qc_fail_not_docked",
    }
)
FROZEN_DOCKING_PROTOCOL = {
    "engine": "AutoDock Vina",
    "version": "1.2.7",
    "scoring_function": "vina",
    "exhaustiveness": 32,
    "num_modes": 20,
    "energy_range": 5,
    "seed": 1701,
    "cpu": 8,
}
FROZEN_LIGAND_PREPARATION = {
    "protocol_id": "gpt_manual_standardized_ph74_meeko080_backbone_rigid_v1",
    "source_coordinates": "preselected_standardized_SDF_no_regeneration",
    "software": "Meeko 0.8.0",
    "rigidified_bonds": "peptide N-CA and CA-C",
}


def _type_matches(value: Any, expected: str) -> bool:
    if expected == "object":
        return isinstance(value, dict)
    if expected == "array":
        return isinstance(value, list)
    if expected == "string":
        return isinstance(value, str)
    if expected == "integer":
        return isinstance(value, int) and not isinstance(value, bool)
    if expected == "number":
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    if expected == "boolean":
        return isinstance(value, bool)
    if expected == "null":
        return value is None
    raise ValueError(f"Unsupported JSON schema type {expected!r}")


def validate_json_schema(value: Any, schema: dict[str, Any], path: str = "$") -> None:
    """Validate the deliberately small JSON-Schema subset used by V0.5 contracts."""
    if "const" in schema and value != schema["const"]:
        raise ValueError(f"{path} must equal {schema['const']!r}")
    if "enum" in schema and value not in schema["enum"]:
        raise ValueError(f"{path} must be one of {schema['enum']!r}")

    expected_type = schema.get("type")
    if expected_type is not None:
        accepted = [expected_type] if isinstance(expected_type, str) else expected_type
        if not any(_type_matches(value, item) for item in accepted):
            raise ValueError(f"{path} has wrong JSON type; expected {accepted!r}")

    if isinstance(value, dict):
        required = schema.get("required", [])
        missing = [field for field in required if field not in value]
        if missing:
            raise ValueError(f"{path} is missing required fields: {missing}")
        properties = schema.get("properties", {})
        additional = schema.get("additionalProperties")
        if additional is False:
            extras = sorted(set(value) - set(properties))
            if extras:
                raise ValueError(f"{path} has unexpected fields: {extras}")
        for key, item in value.items():
            if key in properties:
                validate_json_schema(item, properties[key], f"{path}.{key}")
            elif isinstance(additional, dict):
                validate_json_schema(item, additional, f"{path}.{key}")

    if isinstance(value, list):
        if len(value) < schema.get("minItems", 0):
            raise ValueError(f"{path} has fewer than {schema['minItems']} items")
        if "maxItems" in schema and len(value) > schema["maxItems"]:
            raise ValueError(f"{path} has more than {schema['maxItems']} items")
        if schema.get("uniqueItems"):
            serialized = [json.dumps(item, sort_keys=True, separators=(",", ":")) for item in value]
            if len(serialized) != len(set(serialized)):
                raise ValueError(f"{path} must contain unique items")
        item_schema = schema.get("items")
        if isinstance(item_schema, dict):
            for index, item in enumerate(value):
                validate_json_schema(item, item_schema, f"{path}[{index}]")

    if isinstance(value, str):
        if len(value) < schema.get("minLength", 0):
            raise ValueError(f"{path} is shorter than {schema['minLength']}")
        pattern = schema.get("pattern")
        if pattern and re.fullmatch(pattern, value) is None:
            raise ValueError(f"{path} does not match {pattern!r}")

    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if "minimum" in schema and value < schema["minimum"]:
            raise ValueError(f"{path} must be >= {schema['minimum']}")


def load_json_contract(path: Path, schema_path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot load JSON contract {path}: {exc}") from exc
    try:
        schema = json.loads(schema_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot load JSON schema {schema_path}: {exc}") from exc
    if not isinstance(value, dict) or not isinstance(schema, dict):
        raise ValueError("Contract and schema roots must be JSON objects")
    validate_json_schema(value, schema)
    return value


def validate_candidate_manifest(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    if manifest.get("schema_version") != "0.5":
        raise ValueError("Candidate manifest schema_version must be 0.5")
    candidates = manifest.get("candidates")
    if not isinstance(candidates, list):
        raise ValueError("Candidate manifest candidates must be an array")
    if len(candidates) != 40 or manifest.get("expected_candidate_count") != 40:
        raise ValueError("Formal Top40 manifest must contain exactly 40 candidates")
    if manifest.get("expected_candidate_ids") != list(FORMAL_TOP40_IDS):
        raise ValueError("Formal Top40 expected_candidate_ids are not the canonical A1-D10 inventory")
    expected_counts = {group: 10 for group in FORMAL_GROUPS}
    if manifest.get("expected_group_counts") != expected_counts:
        raise ValueError("Formal Top40 expected_group_counts must be A/B/C/D = 10")

    seen_ids: set[str] = set()
    seen_rows: set[tuple[str, str]] = set()
    counts: Counter[str] = Counter()
    for index, row in enumerate(candidates):
        if not isinstance(row, dict):
            raise ValueError(f"Candidate row {index} must be an object")
        candidate_id = row.get("candidate_id")
        group = row.get("group")
        sequence = row.get("sequence")
        peptide_length = row.get("peptide_length")
        if candidate_id in seen_ids:
            raise ValueError(f"duplicate candidate_id: {candidate_id}")
        seen_ids.add(candidate_id)
        if group not in FORMAL_GROUPS:
            raise ValueError(f"invalid group for {candidate_id}: {group!r}")
        if not isinstance(candidate_id, str) or not candidate_id.startswith(group):
            raise ValueError(f"ID-prefix/group mismatch for {candidate_id!r} and {group!r}")
        if not isinstance(sequence, str) or not sequence:
            raise ValueError(f"empty sequence for {candidate_id}")
        invalid = sorted(set(sequence) - STANDARD_AMINO_ACIDS)
        if invalid:
            raise ValueError(
                f"Sequence for {candidate_id} contains non-standard amino-acid codes: {invalid}"
            )
        if isinstance(peptide_length, bool) or not isinstance(peptide_length, int):
            raise ValueError(f"peptide_length for {candidate_id} must be an integer")
        if peptide_length != len(sequence):
            raise ValueError(
                f"peptide_length mismatch for {candidate_id}: {peptide_length} != {len(sequence)}"
            )
        identity_row = (group, sequence)
        if identity_row in seen_rows:
            raise ValueError(f"duplicate unintended candidate row for group {group}: {sequence}")
        seen_rows.add(identity_row)
        counts[group] += 1

    actual_ids = [row["candidate_id"] for row in candidates]
    if actual_ids != list(FORMAL_TOP40_IDS):
        missing = sorted(set(FORMAL_TOP40_IDS) - set(actual_ids))
        extra = sorted(set(actual_ids) - set(FORMAL_TOP40_IDS))
        raise ValueError(
            f"Formal Top40 candidate order/inventory mismatch; missing={missing}, extra={extra}"
        )
    if dict(counts) != expected_counts:
        raise ValueError(f"Formal Top40 group counts mismatch: {dict(counts)}")
    return candidates


def select_candidates(
    candidates: list[dict[str, Any]],
    groups: list[str] | None = None,
    candidate_ids: list[str] | None = None,
) -> list[dict[str, Any]]:
    if (groups is None) == (candidate_ids is None):
        raise ValueError("exactly one of groups or candidate_ids is required")
    if groups is not None:
        if len(groups) != len(set(groups)):
            raise ValueError("groups contains a duplicate")
        invalid = sorted(set(groups) - set(FORMAL_GROUPS))
        if invalid:
            raise ValueError(f"invalid group selection: {invalid}")
        selected_groups = set(groups)
        return [row for row in candidates if row["group"] in selected_groups]

    assert candidate_ids is not None
    if len(candidate_ids) != len(set(candidate_ids)):
        raise ValueError("candidate_ids contains a duplicate")
    known = {row["candidate_id"] for row in candidates}
    unknown = sorted(set(candidate_ids) - known)
    if unknown:
        raise ValueError(f"unknown candidate IDs: {unknown}")
    requested = set(candidate_ids)
    return [row for row in candidates if row["candidate_id"] in requested]


def validate_screening_protocol(
    protocol: dict[str, Any], *, enforce_frozen_values: bool = True
) -> dict[str, Any]:
    if protocol.get("schema_version") != "0.5":
        raise ValueError("Screening protocol schema_version must be 0.5")
    conformers = protocol.get("formal_conformers")
    if not isinstance(conformers, dict):
        raise ValueError("formal_conformers must be an object")
    count = conformers.get("required_conformer_count")
    names = conformers.get("required_conformer_names")
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise ValueError("required_conformer_count must be a positive integer")
    if not isinstance(names, list) or len(names) != count or len(set(names)) != count:
        raise ValueError("required_conformer_names must be unique and match required_conformer_count")
    if any(not isinstance(name, str) or not name for name in names):
        raise ValueError("required_conformer_names must contain non-empty strings")
    concurrency = protocol.get("batch_execution", {}).get("vina_process_concurrency")
    if concurrency != 1:
        raise ValueError("V0.5 production policy requires vina_process_concurrency = 1")
    if enforce_frozen_values:
        if protocol.get("ligand_preparation") != FROZEN_LIGAND_PREPARATION:
            raise ValueError(
                "Screening protocol differs from frozen V0.4 ligand preparation values"
            )
        if protocol.get("docking_protocol") != FROZEN_DOCKING_PROTOCOL:
            raise ValueError("Screening protocol differs from frozen V0.4 formal Vina parameters")
        if count != 3 or names != ["conf01", "conf02", "conf03"]:
            raise ValueError("Formal V0.5 protocol requires the approved three conformers")
        vina = protocol.get("vina_executable", {})
        if vina.get("version") != "AutoDock Vina v1.2.7":
            raise ValueError("Frozen Vina version must be AutoDock Vina v1.2.7")
        if not isinstance(vina.get("sha256"), str) or not SHA256_PATTERN.fullmatch(
            vina["sha256"].lower()
        ):
            raise ValueError("vina_executable.sha256 must be a SHA-256 value")
    return protocol


def validate_ligand_inventory(
    inventory: dict[str, Any],
    candidates: list[dict[str, Any]],
    protocol: dict[str, Any],
) -> dict[tuple[str, str], dict[str, Any]]:
    if inventory.get("schema_version") != "0.5":
        raise ValueError("Ligand inventory schema_version must be 0.5")
    names = protocol["formal_conformers"]["required_conformer_names"]
    expected = [(row["candidate_id"], name) for row in candidates for name in names]
    slots = inventory.get("slots")
    if not isinstance(slots, list):
        raise ValueError("Ligand inventory slots must be an array")
    if inventory.get("expected_slot_count") != len(expected) or len(slots) != len(expected):
        raise ValueError(
            f"Ligand slot inventory must contain exactly {len(expected)} protocol-derived slots"
        )
    by_key: dict[tuple[str, str], dict[str, Any]] = {}
    for slot in slots:
        if not isinstance(slot, dict):
            raise ValueError("Each ligand slot must be an object")
        key = (slot.get("candidate_id"), slot.get("conformer_name"))
        if key in by_key:
            raise ValueError(f"duplicate ligand slot: {key}")
        if slot.get("readiness_status") not in READINESS_STATUSES:
            raise ValueError(f"Invalid readiness status for ligand slot {key}")
        status = slot["readiness_status"]
        artifact_fields = (
            "source_sdf_path",
            "source_sdf_sha256",
            "audit_path",
            "audit_sha256",
            "ligand_pdbqt_path",
            "ligand_pdbqt_sha256",
            "pdbqt_audit_path",
            "pdbqt_audit_sha256",
        )
        if status in {"not_prepared", "conformer_generation_qc_fail_not_docked"}:
            if any(field in slot for field in artifact_fields):
                raise ValueError(f"Unprepared ligand slot {key} must not fabricate artifact provenance")
            if status == "conformer_generation_qc_fail_not_docked" and slot.get(
                "exclusion_reason"
            ) != "CONFORMER_GENERATION_QC_FAIL":
                raise ValueError(f"QC-excluded ligand slot {key} has the wrong exclusion reason")
        else:
            required_fields = list(artifact_fields[:4])
            if status == "approved_for_exploratory_screening":
                required_fields.extend(artifact_fields[4:])
                required_fields.append("generation_provenance")
            missing = [field for field in required_fields if field not in slot]
            if missing:
                label = "PDBQT/provenance" if status == "approved_for_exploratory_screening" else "provenance"
                raise ValueError(f"Ligand slot {key} is missing {label} fields: {missing}")
            for field in (
                "source_sdf_path",
                "audit_path",
                "ligand_pdbqt_path",
                "pdbqt_audit_path",
            ):
                if field not in slot:
                    continue
                if Path(slot[field]).is_absolute():
                    raise ValueError(f"Ligand slot {key} {field} must be relative to --data-root")
            for field in (
                "source_sdf_sha256",
                "audit_sha256",
                "ligand_pdbqt_sha256",
                "pdbqt_audit_sha256",
            ):
                if field not in slot:
                    continue
                if not SHA256_PATTERN.fullmatch(str(slot[field]).lower()):
                    raise ValueError(f"Ligand slot {key} {field} must be SHA-256")
        by_key[key] = slot
    if list(by_key) != expected:
        missing = [key for key in expected if key not in by_key]
        extra = [key for key in by_key if key not in set(expected)]
        raise ValueError(f"Ligand slot inventory/order mismatch; missing={missing}, extra={extra}")
    return by_key


def load_candidate_manifest(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    value = load_json_contract(path, path.with_name("top40_candidates_v05.schema.json"))
    return value, validate_candidate_manifest(value)


def load_screening_protocol(path: Path) -> dict[str, Any]:
    value = load_json_contract(path, path.with_name("screening_protocol_v05.schema.json"))
    return validate_screening_protocol(value)


def load_ligand_inventory(
    path: Path,
    candidates: list[dict[str, Any]],
    protocol: dict[str, Any],
) -> tuple[dict[str, Any], dict[tuple[str, str], dict[str, Any]]]:
    value = load_json_contract(path, path.with_name("top40_ligand_inventory_v05.schema.json"))
    return value, validate_ligand_inventory(value, candidates, protocol)
