# V0.5-B Batch Execution Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add safe serialized batch execution, immutable candidate attempts, failure isolation, and provenance-aware resume without running real docking during development.

**Architecture:** Keep V0.5-A as the global and candidate planning gate. Add one orchestration module with append-only attempt state and injected execution adapters, reuse V0.4 preparation/docking/metrics/analysis functions, and extend the CLI with execution plus `--resume` while retaining exact dry-run behavior.

**Tech Stack:** Python 3.12, `unittest`, pathlib, JSON/TSV, existing RDKit/Meeko V0.4 modules, mocked/fake subprocess execution.

**Spec:** `docs/v0.5b_batch_execution.md`

## Global Constraints

- Do not run real Vina, C2/D8 docking, Top40 docking, or conformer generation.
- Do not redefine V0.4/V0.5-A scientific semantics or protocol values.
- Execute candidates and conformers serially; one Vina process maximum per workstation.
- Create files exclusively; never delete, overwrite, reuse an invalid attempt, or expose `--force`.
- Do not commit, push, merge, or tag during implementation.

---

### Task 1: State and identity primitives

**Files:**
- Create: `scripts/batch_execution_v05.py`
- Create: `tests/test_batch_execution_v05.py`

**Interfaces:**
- Produces: `CandidateState`, `transition_state(attempt_root, state, detail)`, `build_candidate_identity(...)`, `make_attempt_id(identity, created_at)`.
- Consumes: candidate plan entries and raw-byte contract hashes from V0.5-A preflight.

- [ ] **Step 1: Write failing state-transition and identity tests**

```python
def test_complete_cannot_follow_pending(self):
    with self.assertRaises(ValueError):
        transition_state(self.attempt, CandidateState.COMPLETE, {})

def test_source_sdf_change_changes_identity_digest(self):
    self.assertNotEqual(identity_digest(self.left), identity_digest(self.right))
```

- [ ] **Step 2: Run the focused test and verify RED because the module is absent**

Run: `python -m unittest discover -s tests -p test_batch_execution_v05.py -v`
Expected: import failure for `batch_execution_v05`.

- [ ] **Step 3: Implement the enum, transition graph, canonical JSON hashing, lineage-aware identity, and timestamp-plus-digest attempt ID**

```python
class CandidateState(str, Enum):
    PENDING = "PENDING"
    PREFLIGHT_PASS = "PREFLIGHT_PASS"
    RUNNING = "RUNNING"
    COMPLETE = "COMPLETE"
    FAILED_PREFLIGHT = "FAILED_PREFLIGHT"
    FAILED_DOCKING = "FAILED_DOCKING"
    FAILED_ANALYSIS = "FAILED_ANALYSIS"
```

- [ ] **Step 4: Run the focused test and verify GREEN**

### Task 2: Immutable package and resume validation

**Files:**
- Modify: `scripts/batch_execution_v05.py`
- Modify: `tests/test_batch_execution_v05.py`

**Interfaces:**
- Produces: `create_attempt(...)`, `validate_completed_attempt(...)`, `find_resumable_attempt(...)`, `publish_completion(...)`.
- Consumes: the identity and state primitives from Task 1.

- [ ] **Step 1: Add failing tests for collisions, incomplete/corrupt attempts, every required provenance mismatch, and valid resume**

```python
def test_resume_rejects_corrupt_artifact(self):
    package = make_complete_package(self.root)
    package.output.write_text("changed", encoding="utf-8")
    result = validate_completed_attempt(package.root, package.identity)
    self.assertFalse(result.valid)
```

- [ ] **Step 2: Run and confirm each new test fails at the missing validator behavior**
- [ ] **Step 3: Implement exclusive JSON writes, required-artifact hashes, completion publication, and exact identity comparison**
- [ ] **Step 4: Run and verify all package/resume tests pass**

### Task 3: Serialized execution and failure isolation

**Files:**
- Modify: `scripts/batch_execution_v05.py`
- Modify: `scripts/analyze_peptide_ensemble.py`
- Modify: `tests/test_batch_execution_v05.py`
- Modify: `tests/test_analyze_peptide_ensemble_v04_cli.py`

**Interfaces:**
- Produces: `ExecutionAdapters`, `execute_candidate(...)`, `execute_batch(...)`, and keyword-only safe output-root support for internal batch analysis.
- Consumes: V0.4 `prepare_one`, `audit_one`, `run_one`, `generate_pose_metrics`, and `analyze_peptide_ensemble`.

- [ ] **Step 1: Add failing tests for canonical candidate order, sequential conformers, maximum one fake process active, isolated preflight/docking/analysis failures, global fail-closed, and scientific weakness reaching COMPLETE**
- [ ] **Step 2: Add a failing safety test proving the legacy analyzer CLI still rejects output outside `test_output` while the internal batch call accepts only its declared attempt root**
- [ ] **Step 3: Run and verify RED**
- [ ] **Step 4: Implement the injected adapter pipeline and narrowly scoped analyzer output-root override; keep the CLI default unchanged**
- [ ] **Step 5: Run and verify GREEN, then rerun V0.3/V0.4 ensemble tests**

### Task 4: Production Vina wrapper integrity

**Files:**
- Modify: `tests/test_batch_execution_v05.py`
- Modify: `tests/test_run_standardized_vina.py`
- Modify: `scripts/batch_execution_v05.py`

**Interfaces:**
- Produces: batch-level validation around existing `run_standardized_vina.run_one`.
- Consumes: exact argv builder, run status, PDBQT MODEL parser, and expected protocol `num_modes`.

- [ ] **Step 1: Add failing fake-process tests for success, non-zero exit, missing/malformed output, missing MODEL, and unexpected MODEL inventory**
- [ ] **Step 2: Run and verify RED for batch-level MODEL enforcement**
- [ ] **Step 3: Implement the minimal wrapper validation without `shell=True`, retries, or parameter changes**
- [ ] **Step 4: Run and verify GREEN**

### Task 5: CLI execution and resume

**Files:**
- Modify: `scripts/run_batch_v05.py`
- Create: `tests/test_run_batch_v05_execution_cli.py`

**Interfaces:**
- Produces: required, mutually exclusive `--dry-run`, `--execute`, and `--resume` modes.
- Consumes: `run_dry_run` and `execute_batch`.

- [ ] **Step 1: Add failing tests that preserve dry-run, reject missing/conflicting modes and `--force`, route only explicit `--execute` to new execution, and pass `--resume` only in resume mode**
- [ ] **Step 2: Run and verify RED because execution arguments are unavailable**
- [ ] **Step 3: Implement CLI routing and stable exit codes without future flags**
- [ ] **Step 4: Run V0.5-A and V0.5-B CLI tests and verify GREEN**

### Task 6: Documentation and full release audit

**Files:**
- Modify: `docs/v0.5b_batch_execution.md`
- Test: all `tests/test_*.py`

**Interfaces:**
- Consumes: final observed state, package contract, and CLI behavior.
- Produces: review-ready documentation and verification evidence.

- [ ] **Step 1: Reconcile documentation with final names and on-disk layout**
- [ ] **Step 2: Run all new V0.5-B tests and V0.5-A regressions**
- [ ] **Step 3: Run the  chemistry/project suite, separate PyMOL suite, V0.1 standalone validator, V0.3 legacy A4, V0.4 formal A4, and B3 integration without docking**
- [ ] **Step 4: Audit diff, frozen files, raw/scientific data, secrets, EOL, Git status, and staged area**
- [ ] **Step 5: Perform a requirement-by-requirement review and report `READY_FOR_V05B_REVIEW` only with fresh evidence**
