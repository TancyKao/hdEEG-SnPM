# Analysis catalog — the source of truth

One row per analysis: what it answers, what to import, required columns/roles, the
statistic, its test fixture, and what it outputs. Verified end-to-end by
`matlab -batch "test_all"` (stat-correctness in `test_glm_snpm` / `test_lmm_snpm` /
`test_legacy_snpm`; per-analysis outputs + spectral + edge cases in `test_all`).

## Input sources
Two ways to feed any analysis:
- **CSV files** — wide (one row per subject, channel columns `E1..E256`/`Cz`) for the
  legacy t-tests/correlation/circular; a single combined table for the GLM presets; long
  (one row per trial) for `mixedmodel`. Excel (`.xlsx`+sheet) or `.csv`.
- **Spectral folder** (multi-folder) — one EEG_processor folder per cohort/condition; each
  folder is a level of the factor. Covers `pairedT/onesampleT/unpairedT/anova1/rmanova`
  only (see `scripts/README_scripts.md`).

Channel columns must use real montage labels (`E1..E256`/`Cz` for EGI, `Fp1..` for
Compumedics); set Recording system accordingly. Anything not a channel label is metadata.

## Catalog

| Analysis (key) | Research question | Import | Required cols / roles | Statistic (ref) | Fixture |
|---|---|---|---|---|---|
| **Paired t** `pairedT` | Same people, 2 conditions differ? | 2 wide files (A, B), row-matched | — | paired `ttest` | `synthetic_gui/paired_condA,B.csv` |
| **One-sample** `onesampleT` | Within-subject A−B ≠ 0? | 2 wide files (A, B) | — | `ttest`(A−B) | `synthetic_gui/onesample_condA,B.csv` |
| **Unpaired t** `unpairedT` | 2 independent groups differ? | 2 wide files (A, B) | — | `ttest2` | `synthetic_gui/unpaired_groupA,B.csv` |
| **Correlation P** `correlationP` | Linear EEG↔measure assoc? | 2 wide files, each w/ `Subject` | Subject col | Pearson `corr` | `synthetic_gui/corr_eeg,behavior.csv` |
| **Correlation S** `correlationS` | Monotonic assoc (robust)? | as above | Subject col | Spearman `corr` | as above |
| **Circular** `circ_wheeler_watson_Test`, `circ_WatsonsU2Test` | Phase/angle differs? | 2 wide files, angles in **radians** | — | CircStat `circ_*` | `synthetic_gui/circ_condA,B.csv` |
| **One-way ANOVA** `anova1` | 3+ groups differ? | 1 file | `group_col` (≥2 lvls) | F vs `anova1` | `synthetic_gui/glm_anova1.csv` |
| **ANCOVA** `ancova` | Groups differ adjusting covariates? | 1 file | `group_col` + `covariate_cols` | t (2 groups) / F (≥3, Freedman–Lane) | `synthetic_gui/glm_ancova.csv` (2 groups → t); `synthetic_gui/glm_ancova3.csv` (3 groups → F + pairwise post-hoc) |
| **Regression** `regression` | Scales with a predictor? | 1 file | `predictor_col` (+ `covariate_cols`) | t on slope vs `fitlm` | `synthetic_gui/glm_regression.csv` |
| **RM-ANOVA** `rmanova` | 3+ within conditions differ? | 1 file, long | `subject_col` + `condition_col` | F (within perm) | `synthetic_gui/glm_rmanova.csv` |
| **Two-way mixed ANOVA** `mixed2way` | Condition effect differs between groups? | 1 file, long | `group_col` + `subject_col` + `condition_col` | F on interaction | `synthetic_gui/glm_mixed2way.csv` |
| **Mixed model** `mixedmodel` | Trial-level effect (many obs/subject)? | 1 long file | `Subject` + DV + fixed + effect | `fitlme` t/F (Stephan 2021) | `synthetic_gui/lmm_long.csv` |

## Permutation p-values and sampling

- **Minimum-bias p everywhere.** All tiers (legacy t-tests/correlation, GLM, LMM) report the
  permutation p as `(b + 1) / (N + 1)`, where `b` = permutations at least as extreme as the
  observed statistic and `N` = permutation count. The p is never exactly 0; the smallest
  reportable value is `~1/(N+1)` (the legacy tiers force the observed labeling into the null,
  so their floor is `~2/(N+1)`). A result cannot clear α unless `N` is large enough — e.g.
  `N = 100` floors p near 0.01–0.02, so use `permutations ≥ 5000–10000` for final inference.
- **Exact vs Monte-Carlo (legacy `unpairedT` and correlation).** When the requested
  `permutations` meets or exceeds the exchangeability-group size — `nchoosek(nSubj, nGrp)` for
  the unpaired two-group relabeling, `nSubj!` for correlation — the engine runs an **exact**
  test over the full set of labelings (denominator = the true group size); otherwise it samples
  with replacement (Monte-Carlo). This also removes a former hang when more permutations than
  distinct labelings were requested. *Small-n caveat:* an exact test has a coarse p-grid — e.g.
  unpaired `nGrp = 3, nSubj = 6` → only 20 permutations → minimum p ≈ 0.048 — so significance
  can be unreachable with very small samples regardless of the requested count.
- **Repeated-measures guard.** The between-subject GLM presets (`anova1`, `ancova`,
  `regression`, which use `free` permutation) **error** (`core_snpm_glm:repeatedMeasures`) if a
  subject appears in more than one row, and direct you to `rmanova` or `mixedmodel`. Feed them
  one row per subject.

## Outputs (every analysis, via `core_snpm_analysis`)
- `<base>_<timestamp>.mat` — `results_struct` with `T` (incl. `real_T`), `p`, `Clusters`,
  `uncorrsigch`, `correctTFCEsigch` (TFCE-corrected), `SnPMsigch` (cluster-corrected), `chanlocs`.
- `<base>_<timestamp>.xlsx` — significant-channel table (`func_genSnpmTable`).
- **Correlation runs also emit** each channel's effective N (the pairwise-deletion matched-pair
  count) as `results_struct.per_channel_n` and a guarded **`effectiveN`** sheet in the `.xlsx`;
  the partial-correlation global p uses `df = n − 2 − k` (`k` = number of covariates).
- `<base>_<timestamp>_report.html` — one self-contained report with a **TFCE/Cluster toggle**, stat-aware (t-map for the t-tests, **F-map + post-hoc** table for the omnibus ANOVA presets), the mean topographies and both significance maps. (Replaced the former separate `_TFCE_report.html` / `_Cluster_report.html`.)
- Topoplot PNGs.

**Reading the cluster output:** the GLM/LMM presets score a cluster by **mean Wald**
(`sum(Wald)/n_channels`, height-weighted) while the legacy t-tests score by **extent** (size),
so cluster statistics are **not comparable across the two tiers**; and a significant cluster is
a **regional** claim, not a per-electrode one. See
[Interpreting cluster results](../explanation/interpreting-cluster-results.md).

Specialized reports: GLM faceted report `export_glm_report.m` (anova/t templates in
`templates/`); spectral dashboard `export_report.m`; spectral sweep `SWEEP_grid_*.csv`.

## Choosing an analysis

Design trade-offs (categorical vs continuous predictors, unequal-n / dropout, and the
LMM-vs-GLM cost decision) are discussed in the explanation
[Choosing an analysis: design, missing data, and cost](../explanation/choosing-an-analysis.md).
That page also carries the sleep-question → preset cheat-sheet and the method-grounding references.

## Known gaps (flagged by the tests)
- **Circular analyses are non-functional**: the CircStat functions (`circ_wheeler_watson_test`,
  …) are **not bundled** — they must be added to `dependencies/` for `circ_*` to run.
  `test_legacy_snpm` SKIPs the circular check and says so. **Backlog (pre-fix code path):** the
  two circular engines (`snpm_single_threshold_with_TFCE_circ.m`, `snpm_cluster_analysis_circ.m`)
  still carry the phantom `snpm_enumerate_combinations` reference and the uniqueness-loop hang
  that were fixed in the other tiers. They are CircStat-gated (SKIP) so cannot fire today, but
  the exact-vs-Monte-Carlo / no-hang fix must be extended to them before circular support is
  enabled.
- **Event-bridge group recovery (~25%)**: `test_event_group` is a known XFAIL — its
  `planted_truth` channel indices don't match the analysed 178-subset (a fixture issue; the
  core `anova1` engine is correct per `test_glm_snpm`). Revisit with the TurtleWave importer.
