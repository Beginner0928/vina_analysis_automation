# Top40 Competition Conformer Generation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a deterministic, fail-closed RDKit workflow that can generate exactly three standardized competition conformers per canonical Top40 peptide without invoking Vina or reusing historical conformers.

**Architecture:** A version-controlled JSON protocol freezes scientific parameters and provenance requirements. A focused core module owns chemistry construction/audit, deterministic embedding and MMFF94s optimization, RMSD-based selection, SDF serialization, and external manifests; thin generation, audit, and A4/B3 smoke CLIs call that core. The existing Top40 manifest remains the only candidate identity authority, and generated structures always live outside Git.

**Tech Stack:** Python 3.12, RDKit 2026.03.6, JSON Schema contracts, `unittest`, TSV/JSON/SDF outputs.

**Spec:** User-approved protocol `tfr1_top40_sdfgen_competition_v1_b3derived` in the active task request.

## Global Constraints

- Do not run Vina or create PDBQT/docking artifacts.
- Do not modify V0.1-V0.5-B scientific or execution semantics.
- Do not overwrite historical A4/B3/C2/D8 conformers or an existing production protocol directory.
- Use canonical `specs/top40_candidates_v05.json` as the only candidate identity source.
- Use RDKit 2026.03.6, ETKDGv3, 30 requested source conformers, `useRandomCoords=True`, `pruneRmsThresh=0.25`, `enforceChirality=True`, and `numThreads=1`.
- Use MMFF94s with 1000 iterations; only convergence status 0 is eligible and no fallback is permitted.
- Select exactly three conformers with the fixed heavy-atom aligned RMSD threshold ladder 2.0, 1.5, 1.0, 0.5 A.
- Keep generated SDF and smoke-test artifacts outside Git; do not commit, push, or tag.

---

### Task 1: Freeze the protocol contract

**Files:**
- Create: `specs/competition_conformer_protocol_v1.json`
- Create: `specs/competition_conformer_protocol_v1.schema.json`
- Test: `tests/test_competition_conformer_contract.py`

**Interfaces:**
- Consumes: canonical Top40 manifest hash and frozen external chemistry-contract identity/hash.
- Produces: `load_generation_protocol(path: Path) -> dict[str, Any]` expectations for the core.

- [ ] Write tests that validate the checked-in protocol/schema, reject parameter drift, reject absolute paths, and prove protocol bytes affect the protocol SHA-256.
- [ ] Run the contract tests and observe failure because the new files/core loader do not exist.
- [ ] Add the JSON protocol and schema with the exact frozen values and limitations.
- [ ] Add the minimal loader/semantic validator in the core module during Task 2 and run the contract tests green.

### Task 2: Implement standardized peptide chemistry and audit

**Files:**
- Create: `scripts/competition_conformer_core.py`
- Test: `tests/test_competition_conformer_core.py`

**Interfaces:**
- Produces: `derive_candidate_seed`, `expected_formal_charge`, `build_standardized_peptide`, `audit_standardized_peptide`, `bond_graph_hash`, and `stereochemistry_hash`.
- Consumes: canonical candidate rows and the frozen chemistry contract.

- [ ] Write failing tests for literal deterministic seeds, sequence/length/charge, termini, Asp/Glu/Lys/Arg, HIE-like histidine, free Cys-SH, Tyr-OH, L stereochemistry, and unintended S-S rejection.
- [ ] Run the tests and confirm failure from missing production APIs.
- [ ] Implement construction from `Chem.MolFromSequence`, apply only the frozen formal-charge/hydrogen conventions, add explicit hydrogens, and validate every chemistry invariant without repair-on-audit.
- [ ] Run the focused tests green and refactor only duplicated atom-site lookup code.

### Task 3: Implement deterministic generation and selection

**Files:**
- Modify: `scripts/competition_conformer_core.py`
- Test: `tests/test_competition_conformer_core.py`

**Interfaces:**
- Produces: `build_etkdg_parameters`, `generate_candidate_conformers`, `heavy_atom_aligned_rmsd`, and `select_formal_conformers`.
- Selection consumes energy/status rows and one RMSD callback; it returns original RDKit IDs, energies, the highest successful threshold, and all pairwise RMSDs.

- [ ] Write failing tests for all ETKDG/MMFF settings, convergence-only filtering, energy/source-ID tie breaking, preserved source IDs, deterministic greedy selection, highest successful threshold, and fail-closed behavior below three usable/diverse conformers.
- [ ] Run the tests and confirm each failure reflects missing behavior.
- [ ] Implement ETKDGv3 with a candidate-specific seed, serial embedding/optimization, fixed-index heavy-atom Kabsch alignment via `rdMolAlign.GetAlignmentTransform`, and a fresh greedy attempt at each threshold.
- [ ] Run the focused tests green.

### Task 4: Implement deterministic external outputs and audit

**Files:**
- Modify: `scripts/competition_conformer_core.py`
- Create: `scripts/generate_competition_conformers.py`
- Create: `scripts/audit_competition_conformers.py`
- Test: `tests/test_competition_conformer_cli.py`

**Interfaces:**
- Produces one candidate directory containing `candidate_conf01.sdf` through `candidate_conf03.sdf`, plus protocol-level `preparation_manifest.json`, `generation_summary.tsv`, and `SHA256SUMS.tsv`.
- Audit CLI reads SDFs and external hashes, reconstructs expected chemistry from the canonical manifest, and fails on any mismatch.

- [ ] Write failing tests for names, required SDF properties, external self-hashes, manifest/TSV structure, audit success/failure, and refusal to overwrite an existing protocol directory.
- [ ] Run the CLI tests and observe expected missing-script/API failures.
- [ ] Implement exclusive file creation with LF JSON/TSV, deterministic property ordering, canonical relative paths, and read-back audit before publishing each manifest row.
- [ ] Run CLI tests green and confirm no raw structure fixture is written inside the repository.

### Task 5: Implement dedicated A4/B3 smoke support

**Files:**
- Create: `scripts/smoke_competition_conformers.py`
- Modify: `tests/test_competition_conformer_cli.py`

**Interfaces:**
- Consumes: protocol, canonical manifest, data root, and a fresh external smoke root.
- Produces: an A4/B3-only protocol package; refuses any candidate set other than canonical A4 and B3 and refuses an existing target.

- [ ] Write failing CLI tests proving the fixed A4/B3 subset, external-path requirement, no overwrite, and absence of any Vina import/process call.
- [ ] Run tests red, implement the thin smoke wrapper, then run tests green.

### Task 6: Documentation, smoke reproducibility, and full regression

**Files:**
- Create: `docs/v0.5_competition_conformers.md`
- Test: all new and existing test suites.

**Interfaces:**
- Documents the exact chemistry, seed, RMSD, selection, output, audit, limitation, and recovery contracts.

- [ ] Document the workflow without claiming exhaustive or experimentally validated sampling.
- [ ] Run all new unit/CLI tests in the project chemistry environment.
- [ ] Generate A4 and B3 in two fresh external smoke roots and audit both packages.
- [ ] Compare selected source IDs, energies, RMSDs, protocol/SDF metadata, and raw SDF hashes between runs; diagnose any mismatch without relaxing scientific rules.
- [ ] Run V0.5-B, V0.5-A, chemistry/project, PyMOL, and V0.1-V0.4 no-docking regressions.
- [ ] Verify the feature branch, empty index, `git diff --check`, LF policy, absence of secrets/absolute production paths, absence of tracked raw science data, zero Vina processes, and unchanged historical A4/B3/C2/D8 hashes.
