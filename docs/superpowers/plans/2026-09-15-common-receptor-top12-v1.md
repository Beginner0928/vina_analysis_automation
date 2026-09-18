# Common-receptor Top12 dry-run implementation plan

> **For Codex:** Implement inline with fail-closed checks. This supersedes the multireceptor plan as the formal docking path. Do not delete or extend the retained alignment work. Do not commit or push.

**Goal:** Prepare an isolated 36-job dry-run in which all Top12 ligands use one locked complete group-A TfR1 chain-A receptor and one fixed box, with zero Vina processes.

**Minimal architecture:** One checked-in protocol JSON, one checked-in selection JSON, one thin planner script, and one focused unittest file. Reuse `receptor_registry.resolve_receptor_bundle` and `standardized_vina_inputs.write_vina_config`; read existing ligand inventory without changing chemistry or conformers.

## Frozen inputs

- Receptor PDBQT: `GPT_vina_branch_20260906/receptor_pdbqt/h26d3_af3_s2718_m2_TfR1_A.pdbqt`, locked by registry hash.
- Receptor PDB: `GPT_vina_branch_20260906/receptor_pdb/h26d3_af3_s2718_m2_TfR1_A.pdb`.
- Dimer PDB: `GPT_vina_branch_20260906/receptor_dimer/h26d3_af3_s2718_m2_TfR1_AB.pdb`.
- Center: `[-10.827158, 23.888659, -11.674443]`.
- Size: `[30.591256, 37.457546, 33.838885]`.
- Vina: exhaustiveness 16, num_modes 20, energy_range 5, seed 1701, cpu 8, concurrency 1.
- Selection: A3/A9/A4, B3/B7/B4, C2/C7/C8, D8/D4/D6; each `conf01..conf03`.
- Output root: `E:/research_related/competition/docking_test/common_receptor_top12_v1`.

## Task 1: Contracts and tests

Create `specs/common_receptor_top12_v1.json`, `specs/top12_common_receptor_selection_v1.json`, and `tests/test_common_receptor_top12_v1.py`. Tests assert exact receptor identities/hashes, box, Vina parameters, candidate sequences/order, and 12 x 3 coverage. Verification: `python -m unittest tests.test_common_receptor_top12_v1 -v`.

## Task 2: Thin dry-run planner

Create `scripts/plan_common_receptor_top12_v1.py` with:

```python
def load_contract(path: Path) -> dict[str, Any]: ...
def resolve_jobs(*, protocol: Mapping[str, Any], selection: Mapping[str, Any], inventory: Mapping[str, Any], data_root: Path, output_root: Path) -> list[dict[str, Any]]: ...
def validate_jobs(jobs: Sequence[Mapping[str, Any]]) -> dict[str, Any]: ...
def materialize_dry_run(...) -> dict[str, Any]: ...
def main(argv: Sequence[str] | None = None) -> int: ...
```

The planner resolves registry group A once, verifies the exact user-specified PDB/PDBQT/dimer paths and hashes, verifies candidate sequences against `top40_candidates_v05.json`, verifies all 36 ligand PDBQT paths/hashes from the existing partial35 inventory, and writes 36 configs with `write_vina_config`. It never imports or invokes a Vina runner. It uses exclusive writes and refuses a nonempty/existing output root.

Dry-run report must assert: job count 36; unique receptor PDBQT path/hash count 1; unique center count 1; unique size count 1; unique Vina settings count 1; candidate/conformer coverage 12 x 3; Vina process count 0; state `DRY_RUN_COMPLETE/PREPARED`.

## Task 3: Ranged MODEL validation compatibility

Extend `batch_execution_v05.validate_vina_result` backward-compatibly: old explicit-label callers keep exact validation; the common-receptor protocol can require actual count 1-20, labels exactly `1..N`, no duplicates, and a parseable Vina score/remark for every model. Add focused tests including 19 valid models and invalid zero/21/gaps/duplicates/missing score evidence.

## Task 4: Materialize and verify

Before writing, require output root absent and Vina process count 0. Run the planner with the exact repo/data/output paths. Re-read every config and independently verify all uniqueness counts. Re-hash source receptor and ligand inputs after generation. Run focused existing regressions for registry/config/batch execution and `git diff --check`. Do not run Vina.

Write `md_decision_summary.md` under the output root listing the exact receptor PDB/PDBQT/dimer, fixed center/size, 12 sequences, three conformers, and a note that MD prioritization requires completed docking poses. Final report ends with `READY_FOR_REAL_DOCKING = YES/NO` and `vina_docking_process_count = 0`.
