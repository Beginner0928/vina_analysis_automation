# TfR1 receptor registry and C/D normalization record

Date: 2026-09-07

## Production mapping

`specs/receptor_registry_v04.json` is the only machine-readable V0.4 source
for group switching. A group resolves one indivisible coordinate bundle:
monomer PDB, docking receptor PDBQT, dimer PDB, Vina box, target residues, and
locked hashes.

| Group | Template ID | Monomer SHA-256 | PDBQT SHA-256 | Dimer SHA-256 |
|---|---|---|---|---|
| A | `h26d3_af3_s2718_m2` | `fdc9398d178ebbcd1cd8d33446ca31236e86a7ffa1a413d35281151968ee2bc4` | `f95e4f3eceb2774a711466c3d5690e14e256122f85d42df6dc459ea820f81ca1` | `dd7e1b1eaa1be2fa2dca88da8f51edb2915d1a02891cfa4362f72c6d70c7ed27` |
| B | `h26d3_boltz_msa_s2718` | `52fd892a2bcd296834e8825a46eff7a250f7887d4453708975c4b698880bf21f` | `ee9b2d8a33ed0dd24c23df7c8cf751f3b280e513734444b1a25c7bdd4529038f` | `ee96275e47fc64b1fbba460b8f62791908e43c86985c8a88d5042edd2db19551` |
| C | `h26d3_rank009_rb_s1701_dom` | `ff73d6851c20a4d47c70247a100a5c019372b28f4223afbdd4fed6d6604a99b6` | `efa6f988bbc2384526137b9a963c195b557322529045d84b460dcf11a519af95` | `a51400e0843812a5bfe8fe1fc8cf43079f293e46f697326ec0101a51ba857935` |
| D | `h26d3_rank016_mdref_s1701_r1` | `329a16214e0c0fafafeadc2e7fd9911ebf681f47d171edc0304c9c8e166ada53` | `561ef56c9691cd339593b2ffda6439d551e6b245e93a0271470ed940b04f38d2` | `45b6e859bef755ea9c2a7898f19f2ebd8cc939ac4aaa62d9ca7ca5221acc9fec` |

All paths in the registry are relative to the competition data root. The
audited box source is
`GPT_vina_branch_20260906/config/vina_boxes.tsv`, SHA-256
`e9a1b47227a37ec226cb78c2af8f88104a019b9b66ffd1a10d93d7ab8419fe9b`.

## C/D diagnosis and normalization

The original C/D receptor PDBs failed current Meeko polymer/template matching
for all 48 SER residues after missing element fields were repaired. Heavy-atom
composition and local geometry were intact. Compared with the working B
template, the only reproducible difference was SER ATOM-record order:
`N,CA,CB,OG,C,O` versus `N,CA,C,O,CB,OG`.

The accepted normalization repairs element fields and reorders each SER record
to `N,CA,C,O,CB,OG`. It does not alter coordinates, atom serials/names,
residue identities/numbers, chains, occupancy, or B-factor. It does not delete
any residue. `--delete_bad_res` and `--allow_bad_res` are forbidden because the
affected set includes target residues SER151 and SER159.

Production monomers are therefore exactly:

- C: `GPT_vina_branch_20260906/receptor_pdb_fixed/h26d3_rank009_rb_s1701_dom_TfR1_A_normalized.pdb`
- D: `GPT_vina_branch_20260906/receptor_pdb_fixed/h26d3_rank016_mdref_s1701_r1_TfR1_A_normalized.pdb`

Files containing `SERorder_test` are diagnostic provenance only and are
rejected by production preflight.

Manual review approved both normalized C/D production bundles. The registry
therefore records `normalization_status = PASS`, `production_status = approved`,
and `receptor_preparation_blocker = RESOLVED`. It also records explicit,
independently named SHA-256 provenance for the normalized monomer PDB, docking
receptor PDBQT, and receptor dimer PDB. Preflight requires those provenance
hashes to match both the locked compatibility hashes and the files on disk.

## Acceptance audit

For A/B/C/D, the monomer, receptor PDBQT, dimer, and Vina box were re-resolved
from real local files and hash-checked. Target residues
`150,151,154,158,159,160,161,163,385` are present. Monomer chain A matches
dimer chain A exactly. Monomer/PDBQT heavy-atom coordinates match within
PDBQT serialization precision (<0.001 Å). Group A has an equivalent terminal
residue O/OXT naming swap; coordinate-set matching is exact after that explicit
two-atom equivalence and no other name remapping is allowed.

Registry resolution fails on unknown groups, missing files, hash mismatches,
cross-template components, non-normalized C/D monomers, diagnostic filenames,
or box values that differ from the audited source. There is no silent fallback.
