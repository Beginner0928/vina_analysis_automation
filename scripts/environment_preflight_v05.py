"""Validate the frozen Windows CPython 3.12 runtime used by V0.5."""

from __future__ import annotations

import importlib.metadata
import platform
import re
import sys
from typing import Any


REQUIREMENT = re.compile(r"^([A-Za-z0-9_.-]+)==([^\s;]+)$")


def _normalized_name(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def parse_requirements_lock(text: str) -> dict[str, str]:
    locked: dict[str, str] = {}
    for number, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        match = REQUIREMENT.fullmatch(line)
        if match is None:
            raise ValueError(f"Lock line {number} is not an exact name==version pin: {raw!r}")
        name = _normalized_name(match.group(1))
        if name in locked:
            raise ValueError(f"Duplicate package in dependency lock: {name}")
        locked[name] = match.group(2)
    if not locked:
        raise ValueError("Dependency lock contains no packages")
    return locked


def audit_runtime_environment(
    locked: dict[str, str],
    *,
    python_version: tuple[int, int, int] | None = None,
    installed_versions: dict[str, str] | None = None,
) -> dict[str, Any]:
    version = python_version or (
        sys.version_info.major,
        sys.version_info.minor,
        sys.version_info.micro,
    )
    if version[:2] != (3, 12):
        raise ValueError(f"V0.5 requires Python 3.12.x; actual={version}")
    if installed_versions is None:
        installed_versions = {}
        for package in locked:
            try:
                installed_versions[package] = importlib.metadata.version(package)
            except importlib.metadata.PackageNotFoundError as exc:
                raise ValueError(f"Required package is missing: {package}") from exc
    mismatches = {
        package: {"expected": expected, "actual": installed_versions.get(package)}
        for package, expected in locked.items()
        if installed_versions.get(package) != expected
    }
    if mismatches:
        raise ValueError(f"Locked package version mismatch: {mismatches}")
    return {
        "status": "PASS",
        "required_python_minor": "3.12",
        "python_version": ".".join(str(item) for item in version),
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "packages": {name: installed_versions[name] for name in sorted(locked)},
    }
