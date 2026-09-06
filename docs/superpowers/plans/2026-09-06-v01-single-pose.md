# V0.1 Single-Pose Reproduction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reproduce the manually verified `A4_conf02_out.pdbqt` MODEL 1 analysis headlessly and validate it against immutable golden data.

**Architecture:** A golden JSON is the single regression contract. `analyze_single_pose.py` reads immutable source files, uses the PyMOL API for structural measurements/rendering, and writes a machine-readable result bundle into one new output directory; `validate_single_pose.py` independently compares that bundle with the golden contract and returns a non-zero exit code on any mismatch.

**Tech Stack:** Python 3.12 standard library, open-source PyMOL 3.1.0 Python API, TSV, JSON, Markdown, Git.

**Spec:** `C:/Users/Roy/.codex/attachments/689b1c13-d150-4b3f-884f-3ca3ba13b304/pasted-text.txt`

## Global Constraints

- Do not run AutoDock Vina.
- Do not modify receptor, dimer, ligand, PDBQT, pose-metrics, or the manual `report_package/A4_conf02_m1` package.
- Write generated test artifacts only beneath `test_output/A4_conf02_m1_v01`.
- Refuse to overwrite an existing output directory.
- Do not alter ligand protonation or infer hydrogen bonds from a 4 Å contact.
- Preserve the manually confirmed PDBQT atom IDs exactly.
- Do not commit or push until validation reports `FINAL RESULT: PASS`.

---

### Task 1: Golden contract and failing validator

**Files:**
- Create: `templates/A4_conf02_m1_golden.json`
- Create: `scripts/validate_single_pose.py`

**Interfaces:**
- Consumes: immutable expected metadata, measurements, contacts, tolerances, and output path from the golden JSON.
- Produces: `validation_report.txt` when a result bundle exists; exits non-zero when required analysis outputs are absent or any comparison fails.

- [ ] **Step 1: Create the golden JSON**

Store candidate `A4`, sequence `SPSNFITMYDW`, conformer `conf02`, model `1`, state count `19`, score `-6.470`, target residues, expected contacts, atom IDs `102/103/108/114`, three distances, two angles, and tolerances `0.02 Å`, `0.20°`, and `0.001 kcal/mol`.

- [ ] **Step 2: Write the validator before the analyzer**

The validator must parse JSON and true tab-separated TSV with the standard library, compare sets and numeric values, write `[PASS]`/`[FAIL]` entries plus `FINAL RESULT`, and return `0` only for a complete pass.

- [ ] **Step 3: Prove the missing-results failure**

```powershell
F:\aistuff\Conda\envs\pymol\python.exe scripts\validate_single_pose.py
```

Expected: non-zero exit and a clear missing-results diagnostic; no source data changes.

### Task 2: Headless single-pose analyzer

**Files:**
- Create: `scripts/analyze_single_pose.py`

**Interfaces:**
- Consumes: `templates/A4_conf02_m1_golden.json` and the four immutable source files named there.
- Produces: measurements TSV, contacts TSV, metrics JSON, summary Markdown, run manifest, PSE, and PNG in the configured output directory.

- [ ] **Step 1: Implement input and output safety**

Resolve paths from the golden file, require every input to exist, calculate SHA-256, and fail before launching PyMOL if `test_output/A4_conf02_m1_v01` already exists.

- [ ] **Step 2: Implement PyMOL state extraction and contacts**

Launch `pymol.finish_launching(["pymol", "-cq"])`, load `TfR1_A` and `A4_conf02_all`, require `cmd.count_states(...) == 19`, create `A4_conf02_m1` from state 1, and compute target contacts with a `byres` selection within 4.0 Å.

- [ ] **Step 3: Implement exact atom-ID measurements**

Use PyMOL selectors with ligand IDs `108`, `114`, `102`, and `103`; use receptor chain/residue/atom names exactly as specified; call `cmd.get_distance` and `cmd.get_angle` directly.

- [ ] **Step 4: Parse model-1 pose metrics**

Read the TSV header with `csv.DictReader(delimiter="\t")`, locate model 1, preserve the available score/contact/clash fields, and represent unavailable fields as `NA` without inventing values.

- [ ] **Step 5: Write the result bundle and render**

Write valid UTF-8 tab-separated TSV files, JSON metrics, objective Markdown summary, and manifest; create named measurement objects; show only the Trp11–Gln160 distance by default; save the PSE and render a white-background 2400×1800 ray-traced PNG focused on the pose and actual target contacts.

### Task 3: Integration run and regression validation

**Files:**
- Generate only under: `test_output/A4_conf02_m1_v01/`

- [ ] **Step 1: Run the analyzer headlessly**

```powershell
F:\aistuff\Conda\envs\pymol\python.exe scripts\analyze_single_pose.py
```

Expected: exit `0` and all required result files present.

- [ ] **Step 2: Run regression validation**

```powershell
F:\aistuff\Conda\envs\pymol\python.exe scripts\validate_single_pose.py
```

Expected: `validation_report.txt` ends with `FINAL RESULT: PASS`; if not, investigate code/data provenance without changing golden values.

- [ ] **Step 3: Verify generated artifacts**

Reload the PSE headlessly, decode the PNG, parse both TSV files as tab-separated data, parse JSON, and verify the manifest SHA-256 values against the current immutable inputs.

### Task 4: User documentation and final repository audit

**Files:**
- Create: `docs/v0.1_single_pose.md`

- [ ] **Step 1: Document V0.1**

Describe purpose, fixed inputs, generated outputs, golden contract, exact run/validation commands, overwrite protection, and explicitly unsupported biological or generalized-docking conclusions.

- [ ] **Step 2: Audit Git without committing**

```powershell
git status --short
git diff --check
```

Expected: only intended source/docs/template files appear because `test_output/` remains ignored; no commit or push is performed.
