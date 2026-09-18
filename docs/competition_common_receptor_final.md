# Competition common-receptor final baseline

This competition baseline uses one complete TfR1 chain A receptor for every
completed candidate:

`E:\research_related\competition\GPT_vina_branch_20260906\receptor_pdbqt\h26d3_af3_s2718_m2_TfR1_A.pdbqt`

All jobs use the same docking box and AutoDock Vina settings:

- center: `[-10.827158, 23.888659, -11.674443]`
- size: `[30.591256, 37.457546, 33.838885]`
- AutoDock Vina 1.2.7, scoring `vina`
- exhaustiveness `16`, num_modes `20`, energy_range `5`, seed `1701`, cpu `8`
- direct contact: heavy-atom distance `<= 4.0 Å`

Of the 40 original candidates, 35 completed common-receptor docking with three
conformers each. The Top12 and remaining23 runs comprise 105 completed jobs.
The five excluded candidates (`A8`, `B1`, `D3`, `D7`, and `D10`) are blocked by
upstream conformer/input preparation QC; they are not weak-binder calls.

The unified 35-candidate summary is stored at
`results/common_receptor_all35_candidate_summary.csv`. It is the docking
baseline for subsequent immunogenicity, stability, synthesizability, and final
candidate selection. Older A/B/C/D group-specific receptor results are not a
formal basis for cross-group ranking.
