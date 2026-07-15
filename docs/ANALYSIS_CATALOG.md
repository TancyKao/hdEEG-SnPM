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

## Outputs (every analysis, via `core_snpm_analysis`)
- `<base>_<timestamp>.mat` — `results_struct` with `T` (incl. `real_T`), `p`, `Clusters`,
  `uncorrsigch`, `correctTFCEsigch` (TFCE-corrected), `SnPMsigch` (cluster-corrected), `chanlocs`.
- `<base>_<timestamp>.xlsx` — significant-channel table (`func_genSnpmTable`).
- `<base>_<timestamp>_report.html` — one self-contained report with a **TFCE/Cluster toggle**, stat-aware (t-map for the t-tests, **F-map + post-hoc** table for the omnibus ANOVA presets), the mean topographies and both significance maps. (Replaced the former separate `_TFCE_report.html` / `_Cluster_report.html`.)
- Topoplot PNGs.

Specialized reports: GLM faceted report `export_glm_report.m` (anova/t templates in
`templates/`); spectral dashboard `export_report.m`; spectral sweep `SWEEP_grid_*.csv`.

## Choosing an analysis: design, missing data & cost

**It's one model.** The t-tests, ANOVA, ANCOVA and regression are all special cases of the
linear model `Y = Xβ + ε` — the tool runs them through a single OLS engine (`snpm_glm_stat`).
Only the *coding* of the predictors differs, so pick the preset that matches your variable:
- **Categorical** variable of interest (groups) → two-sample t-test / one-way ANOVA / ANCOVA.
- **Continuous** variable of interest → linear regression.
(A two-group t-test = regression on one 0/1 group dummy; k-group ANOVA = k−1 dummies; ANCOVA =
group dummies + a continuous covariate. The presets do the dummy-coding and choose t vs F for you,
so you don't need to hand-build a regression for group comparisons.)

**Unequal n across conditions (e.g. medical dropout).**
- **Between-subject** designs — two-sample t-test, one-way ANOVA, ANCOVA, regression, correlation —
  tolerate unequal group sizes **by construction**. Dropout just means a smaller group; use freely.
- **Within-subject** designs — paired / one-sample t-test, RM-ANOVA, two-way mixed ANOVA — need the
  *same* subjects across conditions. A subject missing a condition is dropped (**complete-case** →
  fewer subjects, less power). Pairing is by `sub-XX` (spectral folder) or by row order (CSV).

**For substantial within-subject dropout, prefer the LMM** (`mixedmodel`): it uses all available
observations and is valid under **Missing-At-Random (MAR)**. This matters because medical dropout is
rarely random (sicker patients leave); RM-ANOVA's complete-case analysis is only unbiased under the
stronger **MCAR**.

**But the LMM is computationally heavy** — per-channel *iterative* `fitlme` × permutations, roughly
**100–1000× slower** than the vectorized-OLS GLM presets (minutes-to-an-hour vs seconds). Cheaper
routes that usually suffice:
- **RM-ANOVA GLM** — this tool's `rmanova` models subject as a **fixed effect** (subject dummies) by
  OLS, so it's fast and tolerates *moderate* unbalance; you give up the random-effects/partial-pooling
  and the cleaner MAR theory of a true LMM.
- **Aggregate, then GLM** — if the load is *trial-level*, collapse to subject×condition **means** first,
  then run the fast `rmanova`/`mixed2way`.
Reserve the full LMM for genuinely trial-level data, real random-effects structure, or heavy
non-random dropout. Cost knobs: ≤1000 permutations, keep the `parfor` pool, analyse fewer channels.
The LMM is **script-only** (`scripts/run_lmm_analysis.m`) — not yet in the GUI.

## Known gaps (flagged by the tests)
- **Circular analyses are non-functional**: the CircStat functions (`circ_wheeler_watson_test`,
  …) are **not bundled** — they must be added to `dependencies/` for `circ_*` to run.
  `test_legacy_snpm` SKIPs the circular check and says so.
- **Event-bridge group recovery (~25%)**: `test_event_group` is a known XFAIL — its
  `planted_truth` channel indices don't match the analysed 178-subset (a fixture issue; the
  core `anova1` engine is correct per `test_glm_snpm`). Revisit with the TurtleWave importer.
