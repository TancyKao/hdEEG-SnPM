# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Orchestration model (read first)

This repository is maintained by a small squad of specialised subagents defined in
`.claude/agents/`. **The top-level session orchestrates; it does not do the work itself.**
Delegate every substantive task — reading beyond a quick orientation, changing statistics or
GUI code, writing tests or docs, reviewing — to the owning subagent via the Agent tool, then
integrate and report. Do not personally edit the stat engines, the GUI, or the docs when an
agent owns that surface.

Route by ownership:

- **`biostatistician`** — the statistical-methods authority and FIRST leg of method
  development. Decides and justifies *how* to test an effect: the permutation scheme matched
  to the effect, the design/contrast, exchangeability, TFCE-vs-cluster correction, mixed-model
  assumptions. Invoke when the method or its validity is the open question ("is this
  permutation scheme correct?", "how do we test this interaction?", "does this contrast break
  exchangeability?"). Delivers a reviewed **Method Spec** that `stats-engine-engineer` builds;
  prototypes only under `_scratch/`, writes no production MATLAB.
- **`stats-engine-engineer`** — the numerical core: `core_snpm_analysis/glm/lmm.m`, every
  `dependencies/` stat and plot engine, the montage/neighbour layer, result generation
  (Excel/HTML/topoplots), and the `scripts/` headless runners. Anything touching statistical
  correctness or the numerical output of an analysis. Receives Method Specs and implements +
  numerically validates them against MATLAB builtins; does not originate the method.
- **`gui-engineer`** — `SnPMAnalysisGui.m`, the App Designer front-end: analysis-driven
  inputs, role pickers, the Files/Spectral-folder toggle, header auto-detection, readiness
  validation, headless-construct safety, and the ergonomics of driving an analysis for
  non-statistician researchers. Implements how the analysis is *specified*, not the math.
- **`docs-test-engineer`** — the `test_all.m` harness and `test_*.m` suites (statistic-vs-
  builtin, planted-cluster recovery + negative control, output-file assertions), the synthetic
  generators, and `docs/` (especially `ANALYSIS_CATALOG.md`). Invoke after any analysis/GUI/
  engine change to extend tests and sync the catalog.
- **`quality-reviewer`** — the independent quality gate. Reviews every substantive change for
  inferential validity first, then the project's hard constraints, and returns PASS/BLOCK. Not
  the same as an engineer's own self-check; this is a separate, adversarial pass.

## Quality gate (definition of done)

No substantive change is "done" until `quality-reviewer` has reviewed it and returned PASS.
The flow: owning engineer implements → `quality-reviewer` reviews → if BLOCK, the orchestrator
routes each finding back to the owning agent → re-review → PASS. The reviewer is independent
and does not fix its own findings.

For method-driven work there is a front stage: `biostatistician` produces a Method Spec →
`stats-engine-engineer` reviews it for soundness and feasibility → the orchestrator gives
go/no-go, looping the user for live statistical-validity calls → only then does implementation
begin, after which the normal code quality gate applies. A Method Spec is reviewed *before* it
is built on, never handed straight to implementation.

Orchestrator responsibilities stay at the top level: understanding the request, choosing and
sequencing agents (running them in parallel when independent), resolving cross-cutting
decisions, running changes through the gate, and giving the user the final integrated answer.
Only trivial single-surface edits (a typo, a one-line doc fix) may skip the gate. When a task
spans surfaces (a new analysis needs the method decided, the engine built, a GUI input, tests,
and a catalog update), fan out to each owning agent rather than doing any part yourself.

## What this is

MATLAB tool for **Statistical non-Parametric Mapping (SnPM)** of high-density EEG (hd-EEG) topographies. It runs permutation-based, threshold-free cluster enhancement (TFCE) and cluster-based statistics across scalp channels, with channel-adjacency (neighbor) structure. Supports group comparisons, correlations (with optional covariate/partial-correlation control), and circular-statistics variants. Produces topoplots, an Excel table, a `.mat` results struct, and one combined HTML report (TFCE + cluster in a single toggle; stat-aware: t-map for t-tests, F-map + post-hoc for ANOVA).

Theory: Maris & Oostenveld (2007) and Nichols & Holmes (2001) — PDFs in `Ref/`.

**Working approach:** Treat this project as a sleep-EEG **statistician**, not just a coder — prioritize inference validity (permutation scheme matched to the effect, family-wise-error / cluster & TFCE correction, exchangeability, mixed-model assumptions) over engineering polish, and flag when a design or contrast violates them.

## Running

There is no build step. Requires MATLAB (uses `matlab.apps.AppBase` App Designer, `readtable`, stats toolbox functions like `tinv`/`ttest`/`corr`).

- **GUI:** run `SnPMAnalysisGui` from the repo root. It collects file/sheet/output paths and parameters, then calls `core_snpm_analysis`.
- **Headless / NCI scripts (`scripts/`):** edit the `CONFIG` block of a template and run `matlab -batch "run('scripts/<name>.m')"` (no display, exits on completion). `run_glm_analysis.m` (anova1/ancova/regression/rmanova/mixed2way + legacy; joins a subject×channel matrix + `subjects.csv`), `run_lmm_analysis.m` (`mixedmodel`), `run_spectral_report.m` (`export_report`), `run_event_report.m` (`export_event_report`). See `scripts/README_scripts.md` (incl. the `db_to_group_table.py` event-prep step). Verified: `run_glm_analysis.m` runs end-to-end via `matlab -batch` and recovers the planted synthetic cluster.
- **Synthetic test data:** `gen_synthetic_testdata.m` writes `test_data/synthetic_gui/` — one file (GLM presets) or file-pair (legacy/circular) per analysis plus `lmm_long.csv`, each with a planted 17-channel neighbour cluster near E129 (E1…E256 columns → recording system **"EGI 256 (HydroCel)"**). See that folder's `README.md` for the file→analysis→GUI-settings map. Verified: all 10 GUI analyses recover the cluster via TFCE (cluster-extent also for most).
- **Programmatic / smoke test:** edit and run `test_core_snpm.m` — it builds a `params` struct and calls `core_snpm_analysis(params)` directly. This is the fastest way to verify changes without the GUI. Note: paths in `test_core_snpm.m` are hard-coded to the author's Dropbox; repoint `data1_file`/`data2_file`/`output_path` (e.g. to `test_data/`) before running.
- **Verify everything:** `matlab -batch "test_all"` — one PASS/FAIL/SKIP summary ending in `VERIFY: ALL GREEN`: checkcode + headless GUI construct + `test_glm_snpm`/`test_lmm_snpm`/`test_legacy_snpm` (stat-correctness vs MATLAB builtins) + `test_event_group` (known XFAIL) + per-analysis output-file assertions + a synthetic spectral-folder smoke (`gen_synthetic_spectral_folder.m`) + edge/failure-mode tests. Excludes `test_core_snpm.m` (hardcoded paths). Per-analysis import/output/fixture reference: `docs/ANALYSIS_CATALOG.md`. **Known gaps it flags:** circular analyses need the CircStat toolbox in `dependencies/` (currently missing → SKIP); the event-bridge `synthetic_events` `planted_truth` has a channel-index mismatch (XFAIL).
- **Definition of done:** run `core_snpm_analysis` end-to-end and confirm it writes `<base>_<timestamp>.mat`, `.xlsx`, and a combined `<base>_<timestamp>_report.html` to `output_path` without throwing. (The single report is built by `dependencies/generateAnalysisReport.m`, called once per analysis — TFCE+cluster toggle, stat-aware. `export_report.m`'s faceted *spectral* reports are a separate path and unchanged.) Use a small `permutations` value (e.g. 100) for quick iteration; production runs use 10000+.

**Paths must be set from the repo root** — `core_snpm_analysis` calls `addpath(genpath(pwd))` and expects `eeglab2022.1/` (bundled) and `dependencies/` to resolve from the current directory. Override via `params.snpm_path` / `params.eeglab_path` if running elsewhere.

## Architecture

The pipeline is a single orchestrator calling a set of independent stat/plot functions in `dependencies/`. **`SnPMAnalysisGui.m` is only a thin front-end** — all real logic lives in `core_snpm_analysis.m`.

`core_snpm_analysis(params)` flow:
1. **Load** data1/data2 via `readtable` (Excel sheet or CSV).
2. **Correlation-only preprocessing:** match subjects across the two tables by a column containing "subject"; uses *pairwise (per-channel) deletion* of NaNs, NOT listwise — each channel uses all available subject pairs. Requires ≥3 matched subjects. Strips subject columns afterward.
3. **Channel setup** — **recording-system, label-based** (shared helper `dependencies/snpm_setup_channels.m`, driven by `params.channels`). The single source of truth is `dependencies/snpm_montage_registry.m`, which maps a system key to `{labels, chanloc_file, neighbor_file}`. Selection always matches the data file's column names against the system's channel labels (`hdeeg_scalpchannels`), then subsets/reorders data + chanlocs + neighbours into one order — count-agnostic (TurtleWave can export all 256). Systems:
   - `'egi'` — EGI 256 / HydroCel; analyzes the **egi257 178-scalp** set, reusing `EEG178chanlocs.mat` + `NeighborMatrix_178.mat`. Legacy aliases `'164 channels'`/`'178 channels'` both resolve here. **Note:** this replaced the old positional 164-good masking (radius < 0.60 ∩ cluster ≠ 0), so EGI cluster results differ from pre-2026-06 runs.
   - `'compu'` — Compumedics 257 / Neuvo; analyzes the **compu257 249-scalp** set. Assets `compu257_chanlocs.mat` + `NeighborMatrix_compu257.mat` are built from `Compumedics-257.sfp` by `dependencies/build_compu_montage.m`.
   - Add a system by editing `snpm_montage_registry.m` only. The old per-pipeline channel blocks and the positional path are retired; `TopoplotSignificant_single`/`func_genSnpmTable` take a trailing `select_mode` ('label' from all core callers).
4. **Transform** per `params.datatype`: `absolute` (none) / `logscale` (log10) / `normalize` (z-score per subject).
5. **Covariates (correlation only):** if `use_covariates`, residualize both data matrices per-channel against `[ones, covariates]` (see `dependencies/partial_correlation.m`) before permutation — this is partial correlation by residualization.
6. **Stats:** `global_stat_test` (omnibus), then `snpm_single_threshold_with_TFCE` (TFCE-corrected p per channel) and `snpm_cluster_analysis` (cluster-corrected). `circ_*` comparisons route to the `_circ` variants instead.
7. **Outputs:** topoplots (`plot_topoInd`, `TopoplotSignificant_single`), `.mat` struct, Excel (`func_genSnpmTable`), HTML (`generateAnalysisReport`, called once per report_type).

### Key dependency functions (`dependencies/`)
- `snpm_single_threshold_with_TFCE.m` — permutation TFCE; the `compstring = [comparison tail]` switch is the canonical list of supported analyses.
- `snpm_cluster_analysis.m` — cluster-extent permutation test; builds adjacency via `make_neighbors_sparse` and finds clusters with `snpm_find_clusters_graphalgs`.
- `ClusterEnhancement.m` — the TFCE integral (E=0.5, H=2 hard-coded as defaults in `core_snpm_analysis`).
- `global_stat_test.m` — omnibus test on channel-averaged data.
- `func_genSnpmTable.m`, `generateAnalysisReport.m` — Excel and HTML output.
- `*_circ.m` — circular-statistics counterparts for phase/angle data.

### Supported comparisons (`params.comparison`)
`pairedT`, `unpairedT`, `onesampleT`, `correlationS` (Spearman), `correlationP` (Pearson), `circ_wheeler_watson_Test`, `circ_WatsonsU2Test`. `params.tail` ∈ `both`/`left`/`right`. Comparison + tail are concatenated into a `compstring` and dispatched via `switch` inside the stat functions — **adding a new comparison means adding a case in every stat function** (`global_stat_test`, `snpm_single_threshold_with_TFCE[_circ]`, `snpm_cluster_analysis[_circ]`).

### Mixed linear model path (`params.comparison = 'mixedmodel'`)
Per-channel **LMM** for repeated-measures / trial-level designs (e.g. sleep awakenings), replicating Stephan et al. 2021 (`Ref/Stephan et al. - 2021 ...pdf`). It is a **separate pipeline**, not another `compstring` case: `core_snpm_analysis.m` routes `'mixedmodel'` to **`core_snpm_lmm.m`** via an early guard, so the wide-format t-test flow is untouched.
- **Input is long-format** (one row per trial), a single file: `Subject`, group, DV, time, + one power column per channel. `params.lmm_meta_cols` lists the non-channel columns; everything else is treated as channels. Channel count must equal `size(neighbors,1)` (256 for the 164-montage) — guarded with a clear error.
- Engine files in `dependencies/`: `snpm_lmm_fit.m` (per-channel `fitlme`, `parfor` over channels), `snpm_lmm_TFCE.m` and `snpm_lmm_cluster.m` (mirror the output structs of the t-test engines, so `func_genSnpmTable`/`generateAnalysisReport`/`TopoplotSignificant_single` are reused unchanged), `snpm_lmm_permute_meta.m` (permutation schemes).
- **Effect of interest** is parameterized: `params.lmm_effect_type` = `'continuous'` (signed t, e.g. power→DV; supports k≥2 groups as covariates) or `'factor'` (F via `anova(lme)`, e.g. omnibus between-group / interaction — the proper generalization past the 2-group `ttest2` ceiling).
- **Permutation scheme is coupled to the effect** (correctness-critical): `within_subject` (shuffle DV within each subject) for continuous/within effects; `group_label` (relabel whole subjects) for between-group effects. Auto-derived from `effect_type`, overridable via `params.lmm_perm`. The within-subject shuffle gives no null for a between-group mean effect.
- Reuses `ClusterEnhancement`, `make_neighbors_sparse`, `snpm_find_clusters_graphalgs` as-is. Cluster statistic is **mass** (`sum(Wald)/n_channels`, paper's recipe), unlike the t-test path's cluster-**size**.
- Verify with `test_lmm_snpm.m` (`matlab -batch "test_lmm_snpm"`): planted-cluster recovery, factor F-path + negative control, full-pipeline plumbing, and a `pairedT` no-regression check.
- **Not in the GUI** (`SnPMAnalysisGui.m`) — needs formula inputs; run via `scripts/run_lmm_analysis.m` (headless/NCI).

### GLM group-analysis presets (`anova1`/`ancova`/`regression`/`rmanova`/`mixed2way`)
Additive General-Linear-Model engine for sleep designs the legacy t-tests can't express — **>2 groups, ANCOVA, regression with nuisance covariates, repeated-measures ANOVA across >2 conditions, and group×condition interaction**. Like `mixedmodel`, these route from `core_snpm_analysis.m` (early guard) to a separate pipeline **`core_snpm_glm.m`**; legacy paths untouched. Designed for non-statisticians: pick a preset, map columns (`group_col`/`condition_col`/`subject_col`/`predictor_col`/`covariate_cols`), and the engine picks the statistic, contrast, and permutation automatically.
- Engine files in `dependencies/`: `snpm_glm_stat.m` (vectorized OLS across all channels at once — t for a 1-row contrast, F for multi-row; far faster than per-channel `fitlme`, no `parpool`), `snpm_glm_design.m` (preset → `X`/contrast `C`/nuisance/exchangeability blocks/perm type + auto pairwise post-hoc; `get_covariates` dummy-codes categorical/string covariates like `sex` so they can sit alongside numeric ones), `snpm_glm_permute.m` (**Freedman–Lane** residual permutation with exchangeability blocks: `free` for between-subject, `within` for repeated measures), `snpm_perm_correction.m` (shared driver: TFCE max-null + cluster-mass in one loop; emits the standard `T`/`p`/`Clusters` fields).
- Statistic is **t** (1-df: regression, 2-group, pairwise post-hoc) or **F** (omnibus: ≥3 groups, RM-ANOVA, interaction). Omnibus presets auto-run pairwise post-hoc maps (`results_struct.posthoc`).
- **`mixed2way` is a multi-effect report, not a single test.** It routes (early guard in `core_snpm_glm.m`) to the local `run_mixed2way`, which computes **all three** two-way mixed-ANOVA effects — each with its own error term and permutation scheme — and renders them in one combined HTML report: the **group main effect** (between-subjects, computed on subject means via the `anova1` design with `free` permutation — *not* estimable in the within design, so it is a separate computation), the **condition main effect** (within-subject; reduced model `X=[1,Ds,Dc]` **without** the interaction columns so the `Dc` contrast is the pooled/marginal main effect, not the reference-group effect), and the **interaction** (within-subject, full `X=[1,Ds,Dc,Dx]`). It also emits descriptive figures (`plot_cellmean_grid` group×condition cell-mean topos + `plot_interaction_lines` mean±SE interaction plot over the interaction cluster), a parametric **global** whole-head mixed ANOVA per effect (channel-averaged, `compute_global_effect`), and interaction **simple effects** (condition within each group, reusing the `rmanova` design). `params.effect` is **ignored** for mixed2way (the report always covers all three); `results_struct.effects(k)` is a normalized per-effect struct array (`normalize_effect`), with top-level fields kept pointing at the interaction for back-compat (`func_genSnpmTable`/`.mat`). Report rendering: `generateAnalysisReport` has a multi-effect dispatcher (`isfield(rs,'effects')&&numel>1`) → `write_effect_one`/`write_significance_grouped`/`write_descriptive`/`write_global_multi`; single-effect presets fall through unchanged. Reusable `dependencies/save_single_topo.m` was promoted from the core local for the grid helper.
- Method matches the Ref papers: per-channel t/F + neighbour cluster/TFCE + max-statistic permutation null (Maris & Oostenveld 2007; Nichols & Holmes 2001), with Freedman–Lane for nuisance (Winkler et al. 2014).
- Same channel-count invariant and output reuse as the LMM path. Verify with `test_glm_snpm.m` (stat-correctness vs `fitlm`/`anova1`/`anovan`; planted-cluster recovery per preset; 3-group post-hoc; Freedman–Lane confound control; `pairedT` no-regression; **`t8_mixed2way`**: three planted clusters — group/condition/crossover-interaction — each recovered under its own permutation scheme with no cross-leakage, global significance, simple effects, and full-pipeline output/PNG assertions).
- **GUI** (`SnPMAnalysisGui.m`, title "hdEEG-SnPM Toolbox"): "smart", analysis-driven inputs. `ComparisonDropDown` shows **plain-language labels grouped under inert `--- ... ---` headers** (`Items` = labels, `ItemsData` = internal `comparison` keys; header rows `__h*__` revert in `ComparisonDropDownValueChanged`). On analysis change, `applyAnalysisLayout` **relabels the file pickers and hides Data 2** for the single-file GLM presets (Data File vs Condition A/B, Group A/B, Measure 1/2, Angles A/B) and sets a one-line `HintLabel` (label/hint text in the `analysis_labels` local fn). On file load, `detectAndFill` reads headers via `detectImportOptions`, classifies `^E\d+$`/`Cz` as channels vs metadata, auto-fills `MetaColsField` + a `DetectCaptionLabel` ("Detected: N channels, M metadata columns"), and **populates the role pickers from the metadata columns** — Group/Predictor/Condition/Subject are `uidropdown`s, Covariate cols a multi-select `uilistbox` (`fillRoleDropdown`/`firstMatch` helpers; `updateRoleFields` shows only the roles the preset needs via `app.RoleRows`). `checkReadyToRun` is analysis-aware: single Data File for GLM, channel-count must match the montage (256/178), and every visible role must be set. `mixedmodel` stays script-only (needs formula inputs).
- `aggregated-per-channel` event data (e.g. TurtleWave spindle/SO summaries) feeds these presets; `event-level` rows feed `mixedmodel`. A TurtleWave CSV/SQLite importer is a separate future task.

### 178-channel montage + spectral-data bridge
For 178-channel hd-EEG (e.g. the example study), `dependencies/build_178_neighbors.m` builds `NeighborMatrix_178.mat` from chanlocs X/Y/Z (the tool only shipped a 256 matrix); the 178 path in `core_snpm_glm.m`/`core_snpm_lmm.m` loads and remaps it. `load_spectral_dataset.m` assembles an analysis table from a BIDS power-spectrum tree (`sub-*_condition-*_desc-{stage}_powerspect.mat`) — pick band + type (`absolute`/`normalized`); returns `meta_cols` to pass straight to `core_snpm_analysis`. It reads per-channel band power from **`EEG.features`** (falls back to `EEG.bands`) and matches the `'absolute power'`/`'normalized power'` type strings; an optional `opts.subjects_csv` left-joins subject metadata (group/covariates/predictor) on a normalized `Subject` key. `plot_global_spectrum.m` draws the global (channel-averaged) PSD comparison per group. **Note:** condition-a/-b are within-subject → use `rmanova`/`pairedT`, not `anova1`.

**Spectral folder as a data source (multi-folder, folder = one design factor level).** EEG_processor emits one BIDS spectral-power folder per cohort/condition (`<folder>/sub-*_..._desc-<stage>_powerspect.mat`, band power in `EEG.features`). The model: **each folder is one level** of the design factor — you list N folders + labels, pick a band + power type, and the analysis compares them per stage. `load_spectral_dataset.m` takes `opts.folders`(+`opts.labels`) and emits a long table with a generic **`level`** column = the folder label (single-root `condition-*` scan kept only as a fallback when `opts.folders` is absent). `dependencies/spectral_to_snpm_params.m` is the shared helper: loads one (band,type,stage) cell, slices by `level`, writes temp CSV(s), returns `params` for `core_snpm_analysis`.
- **Supported (folder-as-factor):** `unpairedT` (2 folders, between), `pairedT`/`onesampleT` (2 folders, within — subjects matched by `sub-XX`), `anova1` (≥2 folders, between, `group_col='level'`), `rmanova` (≥2 folders, within, `condition_col='level'`). **CSV-only** (need per-subject covariates/predictor or trial-level data): `ancova`, `regression`, `correlationP/S`, `mixed2way`, `mixedmodel`.
- **Headless:** `scripts/run_spectral_analysis.m` — CONFIG = `FOLDERS`+`LABELS` cell arrays + comparison + band(s)/stage(s); band×stage **sweep** → `SWEEP_grid_*.csv`.
- **GUI:** `SnPMAnalysisGui.m` **Data source: Files | Spectral folder** toggle (`DataSourceDropDown`→`applyDataSourceLayout`). Spectral mode shows a **folder `uitable`** (path + editable label) with **Add folder…/Remove selected** buttons, multi-select **Band**/**Stage** listboxes, **Power type** dropdown (nested `SpecGrid`); restricts comparisons to the 5 above; the first added folder fills bands/stages/channels via `scan_spectral_folder`; `runSpectralSweep` builds `FOLDERS`/`LABELS` from the table and loops bands×stages. `checkReadyToRun` validates folder count per analysis (2 for paired/2-group, ≥2 for anova1/rmanova) and unique non-empty labels. The Roles panel is hidden in spectral mode (factor is automatic).

**PSG + KDT file robustness** (`load_spectral_dataset.m` / `scan_spectral_folder`): files are keyed on `sub-XXX` (subject) and `desc-XXX` (stage/segment), ignoring the variable middle (`task-psg`/`task-kdtcomb`, `run`, `condition`) and suffix. Glob is `*powerspect*.mat` (matches PSG `_powerspect.mat` and KDT `_powerspect_final.mat`); stage regex is `desc-([^_.]+)` (PSG `n2`, KDT `eyesopen…`). **Run is ignored** — rows sharing `(Subject, level, stage)` are averaged (no-op for PSG's one-file-per-stage; averages KDT `run-1`/`run-4`). `.mat` files without an `EEG` struct are skipped with a warning (e.g. the the example study KDT `theta_power` folder has 2 files saved as `kdt_comb1`/`kdt_comb4` → dropped). `bandsource` prefers `EEG.bands` then `EEG.features`; both `'absolute'`/`'absolute power'` type spellings match.

### Required data files (loaded by name, must be on path)
`NeighborMatrix_256`, `egi256_chanlocsCluster`, `EEG178chanlocs` — channel locations and the channel-adjacency (neighbors) matrix. The neighbors matrix is channels × max-neighbors, NaN-padded.

## Conventions & gotchas
- `.asv` files are MATLAB autosave backups — ignore them; edit the `.m`.
- Subject ID matching is string-based; numeric IDs are normalized to `sub%03d` format.
- Channel selection is now one shared label-based path (`dependencies/snpm_setup_channels.m`) for every recording system and pipeline — there is no longer a per-pipeline `'164'`/`'178'` branch to keep in sync. Add/adjust systems in `dependencies/snpm_montage_registry.m`. Channel labels come from `dependencies/hdeeg_scalpchannels.m` (`egi257`/`compu257`).
- New comparison types are gated by the GUI dropdown in `SnPMAnalysisGui.m` (`ComparisonDropDown.Items`) AND the stat-function switches; keep them in sync. **Exception:** `'mixedmodel'` (`core_snpm_lmm.m`) and the GLM presets `anova1`/`ancova`/`regression`/`rmanova`/`mixed2way` (`core_snpm_glm.m`) are separate pipelines, not `compstring` cases — don't add them to the stat-function switches.
- Docs live in `docs/` (`PLAN.md`, `RECAP.md`, `MODIFICATIONS_SUMMARY.md`, `DESIGN_PROMPT_*.md`); `docs/MODIFICATIONS_SUMMARY.md` documents the v2.0 addition of Pearson correlation + covariate control (line numbers there may drift from the current source). Report HTML templates live in `templates/` (`anova_report_template.html`, `t_report_template.html`, `sleep_eeg_report.html`, `hdEEG_LMM.html`); `export_glm_report.m` resolves them via `fileparts(mfilename)/templates`.
- `eeglab2022.1/` is a vendored copy of EEGLAB — treat as a third-party dependency, don't modify.
