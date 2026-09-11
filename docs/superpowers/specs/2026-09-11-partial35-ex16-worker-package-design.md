# PARTIAL35 ex16 portable Windows worker package design

Status: approved with the 2026-09-11 variable-model-count amendment; implementation is authorized only in the isolated worktree defined below.

Date: 2026-09-11

## 1. Purpose and review boundary

This design specifies a portable Windows worker for the existing `partial35 competition exploratory ex16` docking workflow:

```text
audited PDBQT + frozen receptor + frozen ex16 protocol
  -> portable worker
  -> explicit candidate subset
  -> serial AutoDock Vina
  -> immutable worker results
  -> host-side isolated import and audit
```

The package is a transport and execution layer for an already frozen scientific protocol. It does not define a new scientific protocol, regenerate conformers, run Meeko, modify ligand or receptor bytes, analyze or rank candidates, or infer work allocation from host completion state.

The initial design received written approval subject to the variable-model-count amendment in this revision. Implementation must occur in an independent Git worktree. No implementation step may stop, pause, kill, restart, resume, modify, or otherwise interfere with the active host ex16 batch.

## 2. Frozen scientific identity

The following values are exact, mandatory, and fail-closed. They are copied into the package manifest and independently validated from the packaged bytes.

| Field | Frozen value |
|---|---|
| `protocol_id` | `tfr1_partial35_competition_exploratory_ex16_v05` |
| Protocol raw-byte SHA-256 | `c94218bce29ff52f8b052ba04f7a6ef5701440ce31c6e76d2d727be7630b37ee` |
| Engine/version | `AutoDock Vina 1.2.7` |
| `vina.exe` source | `E:\research_related\competition\tools\AutoDockVina\vina.exe` |
| `vina.exe` SHA-256 | `e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5` |
| `vina.exe` size | `1,233,920` bytes |
| PE machine | `0x8664`, Windows x86_64 / AMD64 |
| Scoring | `vina` |
| `exhaustiveness` | `16` |
| `num_modes` | `20` |
| `energy_range` | `5` |
| `seed` | `1701` |
| `cpu` | `8` |
| Vina process concurrency | `1` |

`num_modes = 20` is the frozen maximum number of modes requested from Vina, not a requirement that every successful output contain exactly 20 models. A valid output contains `actual_model_count` models where `1 <= actual_model_count <= 20`, with labels exactly `1..actual_model_count` and a parseable `REMARK VINA RESULT` in every model block. This interpretation does not change the scientific parameter or protocol identity.

The current binary is copied byte-for-byte to `vina/vina.exe`; it is never downloaded, updated, patched, or replaced. The worker may not accept a `--vina` override. Package and runtime validation both reject any binary drift.

The package contains the complete canonical allowlist, in canonical order:

```text
A1 A2 A3 A4 A5 A6 A7 A9 A10
B2 B3 B4 B5 B6 B7 B8 B9 B10
C1 C2 C3 C4 C5 C6 C7 C8 C9 C10
D1 D2 D4 D5 D6 D8 D9
```

Each candidate maps to exactly `conf01`, `conf02`, and `conf03`, for 105 ligand PDBQT identities. The package records no host-completion snapshot and has no code that queries host completion to choose candidates.

The four receptor PDBQT identities and boxes are those in the frozen receptor registry:

| Group | Template | Receptor PDBQT SHA-256 | Center (A) | Size (A) |
|---|---|---|---|---|
| A | `h26d3_af3_s2718_m2` | `f95e4f3eceb2774a711466c3d5690e14e256122f85d42df6dc459ea820f81ca1` | `[-11.144, 23.985, -11.463]` | `[28.44, 36.605, 33.415]` |
| B | `h26d3_boltz_msa_s2718` | `ee9b2d8a33ed0dd24c23df7c8cf751f3b280e513734444b1a25c7bdd4529038f` | `[-15.739, 16.133, 16.584]` | `[38.023, 26.643, 27.375]` |
| C | `h26d3_rank009_rb_s1701_dom` | `efa6f988bbc2384526137b9a963c195b557322529045d84b460dcf11a519af95` | `[-3.026, 24.096, -5.173]` | `[30.501, 37.743, 27.644]` |
| D | `h26d3_rank016_mdref_s1701_r1` | `561ef56c9691cd339593b2ffda6439d551e6b245e93a0271470ed940b04f38d2` | `[28.128, 3.3, -9.748]` | `[40.255, 29.623, 30.63]` |

## 3. Non-goals and immutability rules

The worker must not:

- regenerate conformers or require RDKit;
- prepare ligands or require Meeko;
- edit ligand, receptor, protocol, inventory, audit, script, or Vina payload files;
- accept scientific parameter overrides;
- choose candidates automatically, including by reading any host batch state;
- split a candidate's conformers across machines as a normal workflow;
- run more than one Vina process at a time;
- treat partial output as complete;
- overwrite an attempt, output, log, manifest, state event, import, or comparison report;
- import directly into the active host `batch_run`;
- analyze, summarize, rank, or compare ex16 scores scientifically;
- scan or consume the historical ex32 batch as screening data.

The historical ex32 directory may be named only in provenance/readme text explaining its exclusion. No ex32 file, hash, result, completion state, or score is packaged or accepted by the importer.

## 4. Architecture and isolation

### 4.1 Development isolation

After written approval, implementation uses an independent worktree at:

`E:\research_related\competition\vina_analysis_automation_worker_dev`

on branch:

`feat/partial35-ex16-worker-package`

The active checkout remains on its current branch. Implementation must not run `checkout`, `reset`, `clean`, branch switching, or code replacement in the active checkout. The worktree is created from the reviewed design commit. Tests use temporary directories and fake Vina adapters; they do not target the live batch root.

### 4.2 Build isolation

The builder reads the frozen protocol, candidate manifest, ligand inventory, 105 audited ligand PDBQTs, 105 audit records, receptor registry, four receptor PDBQTs, and the exact host `vina.exe`. It writes only to a new staging directory outside both the Git repository and active batch output. All destinations use exclusive creation. Existing destinations cause failure; no force or cleanup mode exists.

The builder validates every source hash immediately before copying and every destination hash immediately after copying. A source that changes between checks causes failure. It never reads candidate completion to determine package contents.

### 4.3 Runtime isolation

The extracted payload is treated as read-only. Runtime files are confined to an explicitly selected output root, defaulting to `results/` beneath the extracted package. A normalized-path guard rejects any path escape through `..`, absolute paths in manifests, symlinks/reparse points, or unexpected drive changes.

The host importer writes only beneath a new exclusive directory:

`<HOST_EX16_ROOT>/incoming_worker_results/<import_id>/`

It treats `<HOST_EX16_ROOT>/batch_run/` as read-only comparison evidence and never publishes into it. Failed or incomplete imports are retained for audit; automated cleanup and deletion are out of scope.

### 4.4 Host A5 interruption provenance

The original host A5 attempt is classified `TECHNICAL_INCOMPLETE_NOT_SCIENTIFIC_FAILURE`. Its `conf01` and `conf02` evidence remains retained, but A5 has no valid candidate completion because `conf03` was interrupted before a run-status/output completion record existed. The approved host resume records A5 as `FAILED_PREFLIGHT` / `RESUME_NO_VALID_COMPLETE_ATTEMPT`, does not rerun `A5_conf03`, and continues from A6. This is an operational provenance fact only: it is not a binding-energy result, biological failure, completed-candidate claim, or package-creation-time allocation rule. A5 remains in the full worker allowlist and may be allocated only by a later explicit human `-Candidates` decision.

## 5. Proposed package layout

```text
partial35_ex16_worker_v1/
  README_WORKER.md
  package_manifest.json
  PACKAGE_MANIFEST_SHA256.txt
  protocol/
    screening_protocol_partial35_ex16_v05.json
    top40_candidates_v05.json
    top40_ligand_inventory_partial35_v05.json
    receptor_registry_v04.json
  ligands/
    A1/A1_conf01.pdbqt
    ... exactly 105 PDBQT files ...
    D9/D9_conf03.pdbqt
  audits/
    A1_conf01_ligand_audit.json
    ... exactly 105 audit records ...
    D9_conf03_ligand_audit.json
  receptors/
    A/h26d3_af3_s2718_m2_TfR1_A.pdbqt
    B/h26d3_boltz_msa_s2718_TfR1_A.pdbqt
    C/h26d3_rank009_rb_s1701_dom_TfR1_A.pdbqt
    D/h26d3_rank016_mdref_s1701_r1_TfR1_A.pdbqt
  vina/
    vina.exe
  scripts/
    preflight_worker.ps1
    run_worker.ps1
    export_worker_results.ps1
  worker/
    __main__.py
    cli.py
    hashing.py
    manifest.py
    pe.py
    preflight.py
    state.py
    runner.py
    output_audit.py
    export.py
  schemas/
    package_manifest.schema.json
    worker_run_manifest.schema.json
    conformer_result.schema.json
    result_bundle_manifest.schema.json
    import_audit.schema.json
  licenses/
    AutoDock-Vina-APACHE-2.0.txt
  THIRD_PARTY_NOTICES.md
```

The release builder produces two external artifacts:

```text
partial35_ex16_worker_v1.zip
partial35_ex16_worker_v1.zip.sha256
```

The ZIP checksum is necessarily external because a ZIP cannot contain its own final hash without self-reference. `worker_package_sha256` is a separate stable payload identity defined below; `archive_sha256` is the raw-byte SHA-256 of the final ZIP and is recorded in the sidecar after ZIP creation.

The `.zip.sha256` sidecar is UTF-8 text containing the lowercase 64-hex digest, two spaces, and the ZIP filename, followed by one LF. The ZIP and sidecar must be transferred together and the ZIP hash must be checked before extraction.

The package includes the complete Apache License 2.0 text and a third-party notice with upstream project `ccsb-scripps/AutoDock-Vina`, version `1.2.7`, license `Apache License 2.0`, packaged path, binary size, PE architecture, and frozen binary SHA-256. Existing upstream copyright or license text is preserved verbatim when present. No license file is removed or rewritten.

## 6. Package identity and manifest contract

### 6.1 Avoiding recursive hashes

`worker_package_sha256` is the SHA-256 of a canonical UTF-8 JSON identity projection containing:

- package schema and package version;
- frozen protocol and Vina identities;
- canonical candidate allowlist;
- sorted payload entries with normalized relative path, role, byte size, and raw-byte SHA-256.

The projection excludes `package_manifest.json`, `PACKAGE_MANIFEST_SHA256.txt`, runtime `results/`, and the final ZIP. This makes the package identity deterministic and independently recomputable without a self-referential hash. `PACKAGE_MANIFEST_SHA256.txt` stores the raw-byte SHA-256 of `package_manifest.json`. The final ZIP SHA-256 is stored only in the external `.zip.sha256` sidecar.

### 6.2 Required package manifest fields

```json
{
  "schema_version": "1.0",
  "package_id": "tfr1_partial35_ex16_worker_v1",
  "package_version": "partial35-ex16-worker-v1",
  "worker_package_sha256": "<canonical payload identity SHA-256>",
  "scientific_identity": {
    "protocol_id": "tfr1_partial35_competition_exploratory_ex16_v05",
    "protocol_sha256": "c94218bce29ff52f8b052ba04f7a6ef5701440ce31c6e76d2d727be7630b37ee",
    "docking_protocol": {
      "engine": "AutoDock Vina",
      "version": "1.2.7",
      "scoring_function": "vina",
      "exhaustiveness": 16,
      "num_modes": 20,
      "num_modes_semantics": "maximum_requested_modes",
      "energy_range": 5,
      "seed": 1701,
      "cpu": 8,
      "process_concurrency": 1
    }
  },
  "candidate_allowlist": ["A1", "...", "D9"],
  "ligands": [
    {
      "candidate_id": "A1",
      "conformer_name": "conf01",
      "path": "ligands/A1/A1_conf01.pdbqt",
      "sha256": "<64 lowercase hex>",
      "size": 0,
      "audit_path": "audits/A1_conf01_ligand_audit.json",
      "audit_sha256": "<64 lowercase hex>"
    }
  ],
  "receptors": [
    {
      "group": "A",
      "receptor_id": "h26d3_af3_s2718_m2_TfR1_A",
      "path": "receptors/A/h26d3_af3_s2718_m2_TfR1_A.pdbqt",
      "sha256": "f95e4f3eceb2774a711466c3d5690e14e256122f85d42df6dc459ea820f81ca1",
      "size": 0,
      "box_center_A": [-11.144, 23.985, -11.463],
      "box_size_A": [28.44, 36.605, 33.415]
    }
  ],
  "vina": {
    "path": "vina/vina.exe",
    "upstream_project": "ccsb-scripps/AutoDock-Vina",
    "version": "1.2.7",
    "license": "Apache-2.0",
    "sha256": "e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5",
    "size": 1233920,
    "pe_machine": "0x8664"
  },
  "scripts": [{"path": "scripts/run_worker.ps1", "sha256": "<hash>", "size": 0}],
  "licenses": [{"path": "licenses/AutoDock-Vina-APACHE-2.0.txt", "sha256": "<hash>", "size": 0}],
  "package_contents": [{"path": "<relative path>", "role": "<role>", "sha256": "<hash>", "size": 0}],
  "excluded_historical_inputs": {
    "exhaustiveness_32": "HISTORICAL_RUNTIME_CALIBRATION_ONLY_NOT_PACKAGED"
  }
}
```

Angle-bracket strings and zero sizes in this design snippet are typed metavariables illustrating the schema shape, not release values. The generated manifest contains only actual lowercase hashes, measured byte sizes, complete arrays, and relative paths; it contains no ellipsis or metavariable text.

The release manifest must contain exactly 35 allowed candidates, exactly 105 distinct `(candidate_id, conformer_name)` ligand rows, exactly three canonical conformers per candidate, exactly four group receptor rows, the packaged Vina row, every executable/runtime script with a hash, all license/attribution files, and every payload file. Duplicate or unlisted payload paths fail validation.

## 7. Proposed command-line interface

All public commands are PowerShell entry points that invoke only Python 3.12 standard-library modules included in the package.

### 7.1 Preflight only

```powershell
.\scripts\preflight_worker.ps1 `
  -Candidates C1,C2,C3,C4,C5 `
  -WorkerId lab-worker-02
```

Preflight never starts a docking job. It may execute `vina/vina.exe --version` with a short timeout solely to verify the frozen version.

### 7.2 New worker run

```powershell
.\scripts\run_worker.ps1 `
  -Candidates C1,C2,C3,C4,C5 `
  -WorkerId lab-worker-02
```

Optional `-OutputRoot` changes only the runtime result location. The default is `<PACKAGE_ROOT>/results`. Candidate spelling and case must exactly match the canonical allowlist. Empty, unknown, repeated, or noncanonical values fail before Vina starts.

### 7.3 Resume an interrupted run

```powershell
.\scripts\run_worker.ps1 `
  -Candidates C1,C2,C3,C4,C5 `
  -WorkerId lab-worker-02 `
  -ResumeRunId 20260911T120000Z_6f96166e-6c71-4c58-81bf-6d5a30cf3071
```

Resume always requires the same explicit candidate list, in the same order, plus the same `WorkerId`. It never infers allocation. It requires the same package identity and machine fingerprint. A different allocation, package, worker identity, or machine fails closed and requires a separate new run and later manual conflict handling.

### 7.4 Export completed or interrupted evidence

```powershell
.\scripts\export_worker_results.ps1 `
  -RunId 20260911T120000Z_6f96166e-6c71-4c58-81bf-6d5a30cf3071
```

Export creates a new immutable result bundle and SHA-256 sidecar. It never changes the run directory and may export a partial run only when explicitly requested with `-IncludeIncomplete`; incomplete artifacts remain visibly non-results.

### 7.5 Planned host-side commands

These are development-repository tools, not worker ZIP runtime dependencies:

```powershell
python scripts/import_partial35_ex16_worker_results.py `
  --result-bundle <worker-result-bundle.zip> `
  --host-ex16-root E:\research_related\competition\docking_test\partial35_exploratory_ex16_v05
```

The builder is likewise host-only:

```powershell
python scripts/build_partial35_ex16_worker_package.py `
  --data-root E:\research_related\competition `
  --output-root <new-exclusive-release-directory>
```

Neither host command has `--force`, `--overwrite`, `--resume`, or a mode that writes into `batch_run`.

## 8. Worker preflight contract

Every new run and every resume executes the complete preflight before any docking process starts. Before each conformer, the worker revalidates the protocol, selected ligand, selected receptor, Vina binary, and current run identity.

Preflight gates run in this order:

1. **Platform gate**
   - Windows x86_64/AMD64 only.
   - CPython major/minor exactly 3.12.
   - Standard library imports only; no RDKit, Meeko, NumPy, Vina Python package, pip, conda, or network access.
2. **Path and lock gate**
   - Package and output roots normalize to allowed locations.
   - Payload paths are relative, nonescaping, and not reparse-point indirections.
   - The selected run directory lock is acquired with an OS-exclusive file handle.
   - Failure to acquire the lock stops the run.
3. **Package-integrity gate**
   - Manifest schema/version/package identity are exact.
   - `PACKAGE_MANIFEST_SHA256.txt` matches the manifest bytes.
   - Canonical `worker_package_sha256` recomputes exactly.
   - All listed payload hashes/sizes match and no required or executable payload is unlisted.
4. **Vina gate**
   - `vina/vina.exe` exists at the fixed path.
   - Size is exactly 1,233,920 bytes.
   - SHA-256 is exactly `e0c4...3cce5`.
   - DOS/PE signatures are valid and PE machine is exactly `0x8664`.
   - `vina.exe --version`, invoked without any config/output argument and with a timeout, reports exactly AutoDock Vina 1.2.7.
   - Any failure stops execution; no alternate binary is searched.
5. **Protocol gate**
   - Protocol file raw-byte hash is exactly `c942...b37ee`.
   - Parsed `protocol_id`, engine, version, scoring, five numeric parameters, process concurrency, conformer names, inventory hash, receptor-registry hash, and candidate-manifest hash match the manifest.
6. **Allocation gate**
   - `-Candidates` is present and nonempty.
   - Every ID exactly matches the allowlist, with no duplicates.
   - For every selected candidate, inventory membership is exactly `conf01/conf02/conf03` and all three remain assigned to this run.
7. **Input gate**
   - Each selected ligand and audit file exists and matches size/hash.
   - Each selected group receptor exists and matches size/hash.
   - Candidate-to-group mapping and group-to-receptor/box mapping are exact.
8. **Output gate**
   - New-run output directory does not exist before exclusive creation.
   - Resume run exists, is not terminal `COMPLETE`, matches allocation and identities, and contains no contradictory COMPLETE records.
   - Sufficient filesystem write access is available. A disk-space warning is recorded; inability to create a small probe file fails before Vina starts. The probe is created in a unique preflight-attempt directory and retained, not deleted.

A preflight report is immutable and records every gate as PASS/FAIL plus observed values. A failure runs zero docking processes.

## 9. State machine and interruption/resume semantics

### 9.1 Storage hierarchy

```text
results/runs/<worker_run_id>/
  run_identity.json
  run.lock
  lock_sessions/<session_id>.json
  sessions/<session_id>/session_start.json
  sessions/<session_id>/session_end.json
  states/<ordinal>_<event_id>_<state>.json
  preflight/<preflight_id>/report.json
  candidates/<candidate_id>/
    states/<ordinal>_<event_id>_<state>.json
    conformers/<confNN>/
      attempts/<attempt_id>/
        attempt_identity.json
        states/<ordinal>_<event_id>_<state>.json
        config.txt
        vina.log.partial
        output.pdbqt.partial
        vina.log
        output.pdbqt
        result_manifest.json
      complete_result.json
  worker_run_manifest.json
```

All identity, state, attempt, and result records use exclusive creation. Every conformer retry gets a new `attempt_id`; an old attempt is never reused, truncated, renamed over, or deleted. Temporary/partial files remain in their attempt directory as interruption evidence.

### 9.2 Run states

```text
ALLOCATED
  -> PREFLIGHT_PASS
  -> RUNNING
  -> COMPLETE
  -> INTERRUPTED
  -> PAUSED_FAILED
  -> FAILED_PREFLIGHT
```

- `FAILED_PREFLIGHT` is terminal for that invocation and starts zero docking processes.
- `INTERRUPTED` and `PAUSED_FAILED` are resumable after a fresh full preflight.
- A power loss may leave the last durable state as `RUNNING`; resume records a new `RECOVERED_INCOMPLETE_SESSION` event and continues without editing the old event.
- `COMPLETE` requires all allocated candidates COMPLETE and a sealed `worker_run_manifest.json`.

Each runner invocation has a new immutable session identity. `session_start.json` records its start timestamp and `session_end.json` records its observed end timestamp and outcome. A hard power loss may legitimately leave a session without `session_end.json`; the next session records that fact in a new recovery event rather than synthesizing or backdating the missing record. Thus every observed end is retained without requiring a mutable run-status file.

### 9.3 Candidate states

```text
PENDING -> RUNNING -> COMPLETE
                   -> INTERRUPTED
                   -> PAUSED_FAILED
```

A candidate becomes COMPLETE only when exactly its three canonical conformers have valid `complete_result.json` records under the same run/package/protocol/receptor/ligand identities. The next candidate cannot start until the current candidate is COMPLETE or the run pauses. This preserves candidate as the minimum allocation unit and prevents normal cross-machine conformer splitting.

### 9.4 Conformer attempt states

```text
PENDING
  -> INPUTS_REVALIDATED
  -> VINA_RUNNING
  -> VINA_EXITED
  -> OUTPUT_AUDIT_PASS
  -> COMPLETE

VINA_RUNNING -> INTERRUPTED_INCOMPLETE
VINA_RUNNING -> FAILED_DOCKING
VINA_EXITED  -> FAILED_OUTPUT_AUDIT
```

Vina writes only to the attempt's `.partial` paths. A conformer is COMPLETE only when all conditions hold:

- process exit code is zero;
- output exists and is nonempty;
- `actual_model_count` is between 1 and the frozen requested maximum `num_modes = 20`, inclusive;
- model blocks are well formed and their labels, in file order, are exactly `1..actual_model_count`, without duplicates, gaps, or extra labels;
- every model block has a parseable `REMARK VINA RESULT` score;
- input, config, executable, log, and output hashes are recorded;
- the result manifest is valid and exclusively created;
- a fresh `complete_result.json` seal points to that immutable attempt.

Promotion from `.partial` uses a rename to a previously nonexistent unique final artifact path. Existing destinations are never replaced. A partial output is never referenced by `complete_result.json`, never exported as a valid result, and never accepted by the importer.

### 9.5 Interruption cases

- **Ctrl+C:** the runner forwards termination to its child Vina process, waits briefly, records `INTERRUPTED_INCOMPLETE` when possible, and exits nonzero. It does not remove partial files.
- **Power loss or forced shutdown:** durable state and partial bytes remain. On resume, a previously `VINA_RUNNING` attempt without a valid completion seal is classified by a new event as incomplete; it is not resumed in place.
- **Sleep/hibernate:** if the process survives, it continues in the same attempt. If it does not, recovery follows the power-loss rule. Wall-clock and monotonic elapsed observations are both recorded when available.
- **Vina nonzero exit or timeout:** the attempt is sealed as failed with exit evidence; output remains partial/nonvalid. Resume creates a new attempt.
- **Runner crash after valid output but before seal:** resume audits the old attempt but does not retroactively modify it or silently promote it. A new explicit recovery result event may seal it only if every immutable piece and identity validates; otherwise a new attempt is created. The recovery decision is recorded.

Already COMPLETE conformers are hash-verified and skipped on resume. They are never overwritten or rerun. Resume revalidates the entire package/protocol identity plus every allocated ligand/receptor before inspecting progress.

### 9.6 Concurrency lock

`run.lock` is opened with an OS-level exclusive, nonshared file handle for the entire runner session. The file may persist after a crash; liveness is represented by the open handle, not by file existence. A new process may reopen the persistent file exclusively only after the old handle has been released by the operating system. Each successful lock acquisition creates a new immutable `lock_sessions/<session_id>.json`. Failure to acquire the handle stops before state changes or Vina execution.

## 10. Worker run and result manifest schemas

### 10.1 Privacy-preserving machine fingerprint

The worker records `machine_fingerprint_sha256`, not raw machine identifiers. It is exactly SHA-256 over the UTF-8 bytes of `partial35-ex16-worker-v1\n<lowercase Windows MachineGuid>\nAMD64`. The MachineGuid is read with the Python standard-library `winreg` module. If it cannot be read, preflight fails rather than inventing an unstable fallback. The raw MachineGuid is held only in memory and is not written. Windows username, computer name, home directory, absolute paths, IP address, MAC address, and user profile are excluded from scientific identity and manifests.

`worker_id` is an operator-supplied nonsecret label matching a conservative character pattern. It is provenance, not a scientific parameter.

### 10.2 Worker run manifest

```json
{
  "schema_version": "1.0",
  "worker_package_version": "partial35-ex16-worker-v1",
  "worker_package_sha256": "<canonical payload identity SHA-256>",
  "worker_run_id": "<UTC basic timestamp>_<UUIDv4>",
  "worker_id": "lab-worker-02",
  "machine_fingerprint_sha256": "<64 lowercase hex>",
  "machine_platform": {"os": "Windows", "architecture": "AMD64", "python": "3.12.x"},
  "candidate_allocation": ["C1", "C2", "C3", "C4", "C5"],
  "protocol_id": "tfr1_partial35_competition_exploratory_ex16_v05",
  "protocol_sha256": "c94218bce29ff52f8b052ba04f7a6ef5701440ce31c6e76d2d727be7630b37ee",
  "vina_sha256": "e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5",
  "start_timestamp_utc": "<RFC3339 UTC>",
  "end_timestamp_utc": "<RFC3339 UTC or null until sealed>",
  "state": "COMPLETE",
  "conformer_results": ["<relative complete-result paths>"],
  "state_history_sha256": "<canonical ordered event digest>",
  "result_contents": [{"path": "<relative path>", "sha256": "<hash>", "size": 0}]
}
```

The sealed final manifest is written once when the run reaches COMPLETE. Before that, `run_identity.json`, immutable session start/end records, and append-only state events are authoritative; no mutable JSON status file is required. An incomplete export contains a separately named `worker_run_snapshot.json` computed from those records and never presents that snapshot as a completion manifest.

### 10.3 Conformer result

Each complete conformer result records:

- worker package/run/worker/machine identities;
- candidate, group, conformer, and attempt IDs;
- protocol ID/hash and exact Vina settings;
- Vina version/binary hash/size/PE machine;
- receptor ID/path/hash and box;
- ligand path/hash and source audit hash;
- generated config bytes/hash;
- start/end timestamps, elapsed observation, command arguments with package-relative paths, and exit code;
- output/log paths, sizes, hashes, model labels, mandatory `actual_model_count`, and Vina scores;
- output audit status and completion seal hash.

No command record contains an absolute local path. Paths in manifests are normalized package- or run-relative paths.

The output-audit portion has the following required shape; `num_modes_requested` remains the frozen scientific setting and is never inferred from the observed count:

```json
{
  "output_audit": {
    "status": "PASS",
    "num_modes_requested": 20,
    "actual_model_count": 17,
    "model_labels": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
    "vina_result_scores": [-7.1, -6.9, -6.8, -6.7, -6.6, -6.5, -6.4, -6.3, -6.2, -6.1, -6.0, -5.9, -5.8, -5.7, -5.6, -5.5, -5.4]
  }
}
```

`actual_model_count` equals the number of well-formed `MODEL`/`ENDMDL` blocks, the length of `model_labels`, and the length of `vina_result_scores`. It must satisfy `1 <= actual_model_count <= num_modes_requested`; labels must equal the integer sequence `1..actual_model_count`. Zero blocks, more than 20 blocks, duplicate/out-of-order/gapped labels, malformed blocks, or a block without a parseable `REMARK VINA RESULT` fail the audit and cannot produce a completion seal.

### 10.4 Result bundle manifest

The exported bundle contains only the selected run's identities, state journal, configs, logs, complete outputs, incomplete evidence when explicitly requested, and a machine-readable bundle manifest. The bundle manifest includes every file path/hash/size and clearly separates:

- `valid_complete_results`;
- `incomplete_attempt_evidence`;
- `failed_attempt_evidence`.

Only `valid_complete_results` may be considered for later manual host acceptance.

## 11. Host import and conflict contract

The importer performs validation before comparison:

1. Safely stage the result bundle under a new `incoming_worker_results/<import_id>/` directory.
2. Reject traversal, links/reparse points, unlisted files, duplicate paths, invalid schemas, or content-hash mismatches.
3. Verify package, protocol, Vina, candidate allocation, machine, ligand, receptor, config, output, and state-history identities.
4. Re-audit each claimed COMPLETE output for exit success, `1 <= actual_model_count <= 20`, exact contiguous labels `1..actual_model_count`, one or more parseable Vina result remarks in every corresponding model block as required by the output schema, and hashes.
5. Reject incomplete/failed evidence as valid screening results while retaining it for audit.
6. Compare each imported candidate against the active ex16 host batch and prior incoming imports. Never inspect ex32 results for completion or comparison.
7. Write one immutable import audit and one candidate comparison report per candidate. Do not write into `batch_run`.

Candidate import classifications are:

- `READY_FOR_MANUAL_ACCEPTANCE`: no host ex16 result or earlier incoming claim exists and all imported evidence passes.
- `DUPLICATE_OR_CONFLICT`: the host ex16 batch or an earlier incoming import already has the candidate. This top-level classification is mandatory even when bytes match.
- `INVALID_WORKER_RESULT`: identity, schema, output audit, allocation, or hash validation fails.
- `INCOMPLETE_WORKER_RESULT`: fewer than the three required COMPLETE conformers are valid.

For `DUPLICATE_OR_CONFLICT`, the report compares at minimum:

- protocol ID and SHA-256;
- candidate allocation and candidate identity;
- all three ligand hashes;
- receptor identity and hash;
- Vina identity and settings;
- each conformer output hash, `actual_model_count`, exact model-label vector, and score vector;
- worker/package/run/machine identities and timestamps.

The report may add a nonauthoritative detail such as `BYTE_IDENTICAL_DUPLICATE` or `CONTENT_CONFLICT`, but the top-level state remains `DUPLICATE_OR_CONFLICT`. The importer never selects the newer result, never resolves a conflict, never marks host completion, and never overwrites either copy. Human acceptance into any scientific result set is a later, separately authorized operation.

## 12. Threat and error model

| Threat/error | Fail-closed behavior |
|---|---|
| Wrong or modified protocol | Raw-byte hash and parsed-value validation fail before Vina. |
| Replaced/damaged Vina | Fixed path, size, SHA-256, PE machine, and version checks fail; no fallback or download. |
| Modified ligand/receptor/audit | Manifest hash mismatch fails before the affected conformer and at every resume. |
| Candidate typo, outsider, duplicate, or lowercase alias | Exact allowlist/allocation validation fails. |
| Attempt to split one candidate | Allocation/state contract requires all three conformers in one run and machine identity. |
| Two runners use one run directory | OS-exclusive lock acquisition fails for the second runner. |
| Crash, Ctrl+C, shutdown, sleep, or Vina failure | Attempt remains immutable and incomplete; no completion seal; resume creates or explicitly recovers an attempt after full revalidation. |
| Truncated/parseable-looking output | Exit, well-formed block, `1..actual_model_count` label, per-model remark, size, and hash gates prevent completion. |
| Disk full or permission loss | Exclusive write fails; run pauses with retained evidence and no overwrite/cleanup. |
| Malicious or malformed ZIP/path traversal | Normalized relative-path checks and safe extraction reject the bundle. |
| Clock drift | UUID provides uniqueness; timestamps are provenance only and never establish result precedence. |
| Username/path leakage | Raw username and absolute paths are excluded; fingerprint inputs are hashed in memory. |
| Same-seed cross-machine numerical differences | Import classifies duplicate candidate results for human comparison; it never assumes byte identity. |
| Antivirus quarantine or blocked executable | Vina preflight fails and instructs the operator to report; no replacement is accepted. |
| Extra DLL/executable/script injection | Unlisted executable/runtime files and package content drift fail package integrity checks. |
| ex32 contamination | Package builder and importer use explicit ex16 paths/identity allowlists; tests assert no historical result enters payload or completion logic. |
| License/attribution omission | Builder requires and hashes Apache-2.0 text and third-party notice before release. |

The package is designed to detect accidental corruption and protocol drift and to constrain normal operator error. It is not a sandbox against an administrator who deliberately modifies both code and manifests; the external ZIP checksum and trusted host-side import audit are the release trust boundary.

## 13. Planned tests

Tests must not launch a production docking job or write to the live batch. All mutable fixtures use temporary directories. Subprocess tests use a fake Vina adapter unless the test is the bounded `--version` metadata check.

### 13.1 Package builder and manifest

- Freeze-test every scientific identity value listed in Section 2.
- Verify the builder selects the exact 35 candidates and 105 audited PDBQTs, independent of host completion.
- Verify exactly three conformers per candidate and four receptor groups.
- Verify source-before-copy and destination-after-copy hashes.
- Verify Vina source path, byte size, SHA-256, `0x8664` PE machine, destination path, and byte identity.
- Verify all scripts, licenses, notices, inputs, and receptors appear once in package contents with hashes.
- Verify deterministic canonical `worker_package_sha256` and manifest raw-byte hash.
- Verify final ZIP raw-byte SHA-256 sidecar after creation.
- Reject existing output/staging/ZIP destinations; verify there is no force path.
- Assert no ex32 artifact path or result is included.

### 13.2 Worker preflight

- PASS on a complete pristine fixture.
- Reject missing/damaged/wrong-version/wrong-size/wrong-PE Vina independently.
- Reject protocol hash or parsed-setting drift, including any one of the five Vina settings.
- Reject unknown, duplicate, empty, noncanonical-case, or malformed candidates.
- Reject any candidate with fewer/more/different conformers.
- Reject ligand, audit, receptor, box, inventory, script, manifest, or license drift.
- Reject Python other than 3.12 and non-Windows/non-AMD64 release execution.
- Verify preflight starts zero docking processes.
- Verify every resume reruns the complete identity/input preflight.

### 13.3 Serial execution and output audit

- Fake-Vina integration test proves maximum concurrent docking processes equals one.
- Prove candidate order and conformer order are deterministic.
- Prove a candidate's three conformers remain in one run/machine allocation.
- Accept exit zero with exactly 20 contiguous models and a parseable Vina result remark in every model block.
- Accept exit zero with exactly 17 contiguous models and a parseable Vina result remark in every model block.
- Accept exit zero with exactly 1 model and a parseable Vina result remark in its model block.
- Reject 0 models and reject 21 models.
- Reject label gaps such as `1,2,4`, duplicate labels, out-of-order labels, malformed/truncated model blocks, and any model missing a parseable `REMARK VINA RESULT`.
- Verify every successful conformer result and host import audit records the exact `actual_model_count` while retaining frozen `num_modes = 20` as the requested scientific setting.
- Verify config contains exact receptor/ligand relative paths, box, and frozen settings.
- Verify outputs/logs/configs and manifests are exclusive-create and hash recorded.

### 13.4 State, interruption, resume, and locking

- Simulate Ctrl+C during Vina and retain partial output without a completion seal.
- Simulate power loss after each state boundary, including after output creation and before sealing.
- Simulate Vina nonzero exit, malformed output, and runner exception.
- Verify COMPLETE conformers are hash-checked, skipped, and never rerun on resume.
- Verify incomplete attempts remain unchanged and a retry uses a new attempt ID.
- Verify resume rejects allocation, worker, machine, package, protocol, ligand, receptor, or Vina drift.
- Verify two processes cannot acquire the same run lock concurrently.
- Verify a persistent lock file is recoverable after the owning handle exits without deletion.
- Verify no state/result file is overwritten during any recovery path.

### 13.5 Result export and host import

- Verify bundle manifest hashes every exported file and separates valid results from incomplete evidence.
- Reject unsafe archive paths, links, duplicate entries, unlisted entries, hash drift, and schema drift.
- Verify importer writes only to a fresh `incoming_worker_results/<import_id>` tree.
- Verify importer never writes to or changes `batch_run`.
- Verify a new candidate becomes only `READY_FOR_MANUAL_ACCEPTANCE`, not host-complete.
- Verify any existing host/imported candidate becomes `DUPLICATE_OR_CONFLICT`.
- Verify both identical and differing outputs retain both copies and produce field-by-field comparison reports.
- Verify protocol, ligand, receptor, Vina, and output comparisons are present.
- Verify incomplete/failed results cannot count as valid screening results.
- Verify ex32 directories are neither scanned nor used by import/completion tests.

### 13.6 Release verification

- Run targeted unit and integration tests in the isolated implementation worktree.
- Run the full repository no-docking test suite.
- Run `git diff --check` and a changed-file review.
- Build into a new external directory only after all tests pass.
- Recompute and audit every extracted payload hash from the produced ZIP.
- Compare packaged `vina/vina.exe` byte-for-byte with the frozen host binary.
- Record ZIP path, size, SHA-256, package identity, test counts, and build timestamp.
- Confirm the active host Python/Vina PIDs were not controlled or restarted and live batch files were not modified.

## 14. README requirements

`README_WORKER.md` must clearly separate:

**Worker runtime requirements**

- Windows x86_64/AMD64;
- Python 3.12;
- extract the ZIP;
- run preflight with an explicit candidate subset;
- run or resume with the same explicit allocation.

**Not required on the worker**

- separate AutoDock Vina installation;
- pip/conda environment creation;
- RDKit, Meeko, NumPy, or host analysis dependencies;
- conformer generation or ligand preparation;
- network access.

The README must state that the package includes fixed AutoDock Vina 1.2.7, replacement is prohibited, and any missing file, damage, version mismatch, or hash mismatch requires stopping and reporting. It must explain candidate-level allocation, serial execution, safe resume, immutable attempts, result export, and the fact that host import is manual and conflict-aware.

## 15. Acceptance criteria

The design is implementable when all of the following are demonstrable without touching the active host batch:

1. The package binds exactly the frozen ex16 identity and contains no alternate scientific parameters.
2. The ZIP is offline-capable with Python 3.12 standard library and its bundled exact Vina binary.
3. Every run requires an explicit valid candidate subset and maps each selected candidate to exactly three serial conformers.
4. Interruption and resume preserve all evidence, skip validated COMPLETE conformers, and never overwrite.
5. Package, run, conformer, bundle, and import identities are machine-readable and hash-auditable.
6. Concurrent access to one run directory fails closed.
7. Import is isolated, never auto-merges, and flags every duplicate candidate as `DUPLICATE_OR_CONFLICT`.
8. ex32 history is excluded from payload, completion, and ranking semantics.
9. Vina licensing and attribution are included without modifying the binary.
10. All tests and release verification run outside the active batch worktree/output.

## 16. Resolved design decisions and open questions

Resolved decisions:

- Use a purpose-built standard-library Python worker, not the complete V0.5 batch runner.
- Include all 35 candidates and 105 audited ligands; allocation remains a required runtime argument.
- Candidate is the minimum allocation unit.
- Resume retains incomplete attempts and creates a new attempt unless a fully valid unsealed result can be explicitly recovered.
- Use a canonical payload identity for `worker_package_sha256` and an external sidecar for the final ZIP raw-byte hash.
- Keep host import isolated and nonauthoritative.

Open questions: none required for implementation. Any later request to change a scientific value, split candidates across machines, accept a different Vina binary, add automatic host merge, or include ex32 data is a protocol/scope change and requires a new reviewed design rather than an implementation-time choice.
