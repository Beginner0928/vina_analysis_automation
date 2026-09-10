# PARTIAL35 Exploratory Screening Implementation Plan

> Execute gate by gate. Stop downstream work on any candidate-set, protocol, chemistry, artifact-count, hash, receptor, Vina, collision, or provenance failure.

## Scope and invariants

- Generate and screen exactly the canonical 35-candidate set in `specs/partial35_exploratory_package_v05.json`.
- Preserve `tfr1_top40_sdfgen_competition_v1_b3derived` and its SHA-256 byte-for-byte.
- Keep A8, B1, D3, D7, and D10 classified only as `CONFORMER_GENERATION_QC_FAIL_NOT_DOCKED`.
- Never write the partial package into the reserved complete-Top40 protocol directory.
- Use the existing competition conformer scientific primitive, frozen ligand preparation, receptor registry, Vina settings, V0.5-B serial execution, and V0.4 analysis.
- Require explicit `--execute` or `--resume`; dry-run executes zero Vina processes.

## Phase 1: package contract and generation wrapper

Files:

- Add `specs/partial35_exploratory_package_v05.json` and schema.
- Add `scripts/partial35_contract_v05.py`.
- Add `scripts/generate_partial35_conformers.py`.
- Add `scripts/audit_partial35_conformers.py`.
- Add contract and CLI tests.

Implementation:

1. Validate the package contract against its schema.
2. Load canonical candidates only through the frozen competition protocol authority.
3. Reject duplicates, unknown IDs, wrong included/excluded sets, non-canonical order, or a completeness claim.
4. Compute the candidate-set digest from canonical JSON identity rows in canonical manifest order.
5. Generate each candidate sequentially with `generate_candidate_conformers` and publish through the existing immutable package writer.
6. Add a separate immutable `partial_package_manifest.json` containing package identity, parent protocol/hash, manifest/hash, chemistry contract/hash, included/excluded sets, generator Git commit, runtime versions, and the explicit `NOT_A_COMPLETE_TOP40_SCREEN` statement.
7. Audit both the frozen generation package and partial-package contract.

## Phase 2: PDBQT package and V0.5-B bridge

Files:

- Add `scripts/prepare_partial35_ligands.py` and `scripts/audit_partial35_ligands.py`.
- Extend the ligand-inventory schema and validator additively.
- Extend V0.5-B preflight to verify optional preprepared PDBQT bytes and preparation audits.
- Extend standardized preparation so an explicitly declared, hash-locked PDBQT is copied and re-audited instead of regenerated.
- Include preprepared PDBQT identity in the immutable candidate/resume identity.
- Add unit/integration tests, including N-terminal Pro chemistry audit coverage.

Implementation:

1. Audit Gate B before preparation.
2. Prepare exactly 105 PDBQT files with the frozen Meeko function; do not touch the source SDFs.
3. Publish one immutable audit record per conformer, a hash inventory, a full 120-slot derived ligand inventory (105 approved exploratory slots and 15 excluded slots), and a derived screening protocol.
4. Require all 105 PDBQT files to pass readability, mapping, coordinate, source-SDF, preparation-protocol, package, and hash checks.
5. During V0.5-B preflight, require external PDBQT validation for exploratory-ready slots.
6. During execution, copy the already-audited PDBQT into the immutable attempt and rerun the normal pre-docking audit. Never invoke Meeko again for these slots.

## Phase 3: group-wise summary

Files:

- Add `scripts/summarize_partial35_screen.py` and tests.

Implementation:

1. Read only valid V0.5-B COMPLETE attempts and their one-row candidate summaries.
2. Require the exact approved 35-candidate identity and exact excluded set.
3. Produce per-group TSV/JSON/Markdown summaries with screened/canonical counts, score distributions, eligible fraction, recurrent contacts, representatives, warnings, and status.
4. Never emit a cross-group global rank and label Vina scores as predictions.

## Phase 4: patch freeze

1. Run targeted tests first, then all authorized no-docking regressions in the chemistry and PyMOL environments.
2. Verify the frozen protocol hash, frozen historical artifacts, LF policy, raw-data exclusions, and `git diff --check`.
3. Commit exactly the patch as `feat: add partial35 exploratory screening workflow` and push the feature branch.

## Phase 5: authorized production gates

1. Gate A: generate 35 candidates into `standardized_ligands/partial_packages/competition_v1_partial35_mmff_qc_pass`.
2. Gate B: require 35 candidates and 105 audited SDF files.
3. Prepare and audit 105 PDBQT inputs in a distinct external partial35 Vina-input package.
4. Gate C: require 105/105 PDBQT PASS and publish the derived inventory/protocol.
5. Run exact partial35 `--dry-run`; require zero Vina processes.
6. Run exact partial35 `--execute` serially. If interrupted, validate and use `--resume`; never restart completed attempts.
7. Validate all technical outputs and produce the group-wise exploratory summary.
