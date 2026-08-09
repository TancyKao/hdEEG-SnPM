# SnPM Sleep EEG — Statistical Analysis Expansion Plan

Consolidated plan + status for extending the per-channel SnPM tool to cover the full range of
sleep-research group statistical designs. Two efforts, both **implemented and verified**, plus
deferred items.

- **Effort 1 — Mixed Linear Model (LMM) path** — ✅ done (`mixedmodel`)
- **Effort 2 — GLM group-analysis presets** — ✅ done (`anova1`/`ancova`/`regression`/`rmanova`/`mixed2way`)
- **Deferred** — GUI wiring; TurtleWave event-table importer

---

## Context & problem

The original tool tests each EEG channel with a **t-test / correlation** where the exchangeable unit is
the **subject** (one value per subject per channel, wide format), dispatched through a duplicated
`comparison+tail` switch in three engines (`global_stat_test.m`, `snpm_single_threshold_with_TFCE.m`,
`snpm_cluster_analysis.m`) plus two `_circ` variants. It could not express what real sleep studies need:
repeated measures (many awakenings/events per subject), >2 groups, group×condition interactions, or
regression with nuisance covariates.

**User decisions:** purely **additive** (do not replace/refactor the legacy tests); **non-statistician
friendly** (named presets + column mapping, no formulas, statistic/contrast/permutation chosen
automatically); runs on **NCI/HPC**; statistics stay in **MATLAB** (EEG power is `.mat`; the entire
correction/topography/report stack is MATLAB), with Python event detection (TurtleWave) upstream and the
**file** (CSV/SQLite) as the bridge.

## Architecture

Three tiers feed the **same** existing correction + output chain
(`ClusterEnhancement`, `make_neighbors_sparse`, `snpm_find_clusters_graphalgs`, the max-statistic →
corrected-p pattern, and `func_genSnpmTable`/`generateAnalysisReport`/`TopoplotSignificant_single`,
all keyed on `T.real_T`/`p.real`/`p.correctedTFCE`/`Clusters(i).p`/`.channels`):

- **Legacy** (unchanged): `pairedT`, `unpairedT`, `onesampleT`, `correlationP/S`, `circ_*`.
- **Tier A — GLM + permutation** (fast, vectorized; fixed-effects designs).
- **Tier B — LMM** (per-channel `fitlme`; trial/event-level random effects).

Routing: `core_snpm_analysis.m` has early guards that send `mixedmodel` → `core_snpm_lmm.m` and the GLM
presets → `core_snpm_glm.m`; the legacy switch is the untouched fall-through. Channel-count invariant:
`size(power,2) == size(neighbors,1)` (256 for the 164-montage), guarded with a clear error.

### Legacy per-channel statistics (unchanged)

| comparison | per-channel test | MATLAB | permutation |
|---|---|---|---|
| `pairedT` | paired t (x vs y) | `ttest` | sign-flip within subject (exhaustive 2^n or capped) |
| `unpairedT` | two-sample t | `ttest2` | group-label shuffle (`randperm`, up to `nchoosek`) |
| `onesampleT` | one-sample t vs 0 | `ttest` | sign-flip |
| `correlationP` / `correlationS` | Pearson / Spearman r | `corr` | subject relabel |
| `circ_*` | circular tests | `_circ` engines | relabel |

All legacy tests: exchangeable unit = **subject**; FWE via max-statistic TFCE + cluster size; tail ∈ both/left/right.

### Design → GLM building blocks (Tier A)

Every GLM preset is the same computation — fit `Y = Xβ + ε`, test contrast `C` — differing only in how
`X`, `C`, and the exchangeability structure are built:

| Design | X columns | Contrast | Permutation |
|---|---|---|---|
| one-sample / paired t | intercept (on diffs) | intercept | sign-flip within subject |
| unpaired t | group dummy | group | shuffle across subjects |
| one-way ANOVA (k groups) | k−1 dummies | F | shuffle across subjects |
| regression / correlation | predictor (+ nuisance) | slope | Freedman–Lane |
| ANCOVA (group + age/sex/AHI) | group + nuisance | group | Freedman–Lane |
| two-way mixed ANOVA (group × condition) | group, condition, G×C | F on interaction | within / whole per term |
| RM-ANOVA (N2/N3/REM) | subject + condition | condition | within-subject |
| serial awakenings / events | + random(1\|subj) | term | LMM tier (Effort 1) |

## Research question each analysis answers (SWA = slow-wave activity, running example)

| Analysis | Question it answers | Sleep example |
|---|---|---|
| one-sample t *(legacy)* | Is SWA per channel ≠ zero / a fixed reference? | Is the overnight SWA change ≠ 0? |
| paired t *(legacy)* | In the same people, does SWA differ between **two** conditions? | Baseline vs recovery; early vs late night; pre/post-CPAP |
| unpaired t *(legacy)* | Does SWA differ between **two** groups? | OSA vs controls |
| one-way ANOVA `anova1` | Does SWA differ across **3+** groups? then which pair | Control vs mild vs severe OSA; 3 insomnia phenotypes |
| ANCOVA `ancova` | Do groups differ **after removing** age/sex/AHI? | Patient–control gap once age partialled |
| regression `regression` | Does SWA **scale with** a continuous variable (± covariates)? | SWA vs AHI; spindle density vs memory gain; power vs sleep depth |
| RM-ANOVA `rmanova` | Within subjects, does SWA differ across **3+** conditions/stages? | N2 vs N3 vs REM; across 4 NREM cycles |
| two-way mixed ANOVA `mixed2way` | Does the **condition effect differ between groups**? (group × condition) | Is overnight SWA dissipation slope different in insomnia? Treatment×group |
| correlation P/S *(legacy)* | Monotonic power–behaviour association | (subset of regression) |
| mixed model `mixedmodel` | Event/trial-level effect with **many observations per subject** | Per-awakening fast power vs sleep depth (Stephan); per-spindle amplitude vs memory |
| circular *(legacy)* | Phase questions | Spindle phase relative to SO; coupling-phase group difference |

---

## Effort 1 — Mixed Linear Model (`mixedmodel`)  ✅ done

Replicates **Stephan et al. 2021** (`Ref/Stephan ...pdf`): per-channel `fitlme`, effect of interest =
signed **t** (continuous predictor, e.g. power→sleep-depth; supports k≥2 groups as covariates) or **F**
(`anova`, omnibus between-group / interaction). Permutation scheme is **coupled to the effect**:
`within_subject` (shuffle DV within subject) for continuous/within effects; `group_label` (relabel whole
subjects) for between-group effects. Correction = TFCE (added) + cluster **mass** (`sum(Wald)/n`, paper's
recipe).

**Effect of interest (per channel):**

| Effect of interest | df | Statistic | Sign? | Tails |
|---|---|---|---|---|
| power (continuous predictor) | 1 | t, Wald = t² | signed | both/left/right |
| group, 2 levels | 1 | t (group contrast) | signed | both/left/right |
| group, ≥3 levels (omnibus) | k−1 | F via `anova(lme)` | non-negative | both only |
| power × group interaction | k−1 | F | non-negative | both only |

**Files:** `dependencies/snpm_lmm_fit.m`, `snpm_lmm_TFCE.m`, `snpm_lmm_cluster.m`, `snpm_lmm_permute_meta.m`;
`core_snpm_lmm.m`. Input = long format (one row per trial; `lmm_meta_cols` separates
meta from channel columns). Needs Statistics & ML Toolbox; `parfor` over channels (start `parpool` on NCI).

**Verified** (all tiers): planted neighbour cluster recovered (TFCE + cluster); factor
F-path detects a real group effect with a clean negative control; full CSV→pipeline plumbing writes
`.mat`/`.xlsx`/HTML; `pairedT` no-regression.

## Effort 2 — GLM group-analysis presets  ✅ done

One **vectorized GLM** (`β = X\Y` across all channels at once — far faster than per-channel `fitlme`, no
`parpool`) behind friendly presets. Each preset auto-builds design `X`, contrast `C` (1 row → t, >1 → F),
nuisance columns, exchangeability blocks, and permutation type.

| Preset | Sleep use case | Statistic | Permutation |
|---|---|---|---|
| `anova1` | ≥2 groups (control / treatment-A / treatment-B) | omnibus F + auto pairwise post-hoc | shuffle subjects across groups (FL) |
| `ancova` | group difference controlling age/sex/AHI | F/t on group, nuisance partialled | Freedman–Lane (permute, keep nuisance) |
| `regression` | power vs continuous (AHI, sleep-depth, cognition) ± covariates | t on slope | Freedman–Lane |
| `rmanova` | ≥2 within conditions (N2/N3/REM, baseline/recovery) | omnibus F + post-hoc | permute condition labels **within subject** |
| `mixed2way` | two-way mixed ANOVA: group × condition (does patient–control gap depend on stage) | interaction (default) or condition F | within-subject (per-term EB) |

A `transform` option (`absolute`/`log`/`rank`) applies across presets; `rank` yields nonparametric equivalents.

**Files:** `dependencies/snpm_glm_stat.m` (vectorized t/F), `snpm_glm_design.m` (preset→design + post-hoc),
`snpm_glm_permute.m` (Freedman–Lane residual permutation with exchangeability blocks), `snpm_perm_correction.m`
(shared driver: TFCE max-null + cluster-mass in one loop, emits standard `T`/`p`/`Clusters`);
`core_snpm_glm.m`. User maps columns via `group_col`/`condition_col`/`subject_col`/
`predictor_col`/`covariate_cols`; `meta_cols` separates meta from channel columns.

**Verified** (7 tiers): GLM t/F match `fitlm`/`anova1`/`anovan` to 1e-6; planted-cluster
recovery per preset with clean negative controls; 3-group omnibus + post-hoc identifies the differing pair;
Freedman–Lane keeps a true predictor effect while suppressing an age-confound cluster; `pairedT`
no-regression. Also fixed a latent `titlename` typo in `TopoplotSignificant_single.m` that was silently
killing topography figures for non-"VS" titles (repairs both new paths).

## Method grounding (`Ref/`)

- **Maris & Oostenveld 2007** — per-sensor statistic (**t for two conditions, F for >2**), neighbour
  clustering, cluster-level statistic, max-cluster permutation null for FWER; scheme set by design.
- **Nichols & Holmes 2001** — single-threshold **maximum-statistic** permutation test; **sign-flipping**
  for one-sample/within.
- **Winkler et al. 2014** (PALM) — **Freedman–Lane** permutation for the GLM with nuisance regressors.
- **Stephan et al. 2021** — the event/trial-level LMM tier.

---

## Event data (TurtleWave) — both granularities

- **Aggregated per channel** (mean spindle amplitude / density / count per subject×channel) → wide table →
  GLM presets.
- **Event-level** (one row per detected event) → long table → `mixedmodel`.

No new engine needed for either. The **TurtleWave→analysis-table importer is deferred** to its own task
(pending a sample CSV/SQLite); MATLAB reads CSV via `readtable` and SQLite via the `sqlite` interface.

## Real dataset format (local-sleep study, hd-EEG spectral power)

Path: `…/<study>/01_data/derivatives/eeg/hdeeg_analysis_all_sub/{condition-a,condition-b}/`.
One `.mat` per **subject × condition × stage × run**, BIDS-named
`sub-XX_condition-{a,b}_task-psg_run-N_desc-{n1,n2,n3,rem}_powerspect.mat`
(~27 subjects/condition; 2 conditions × 4 stages).

Each file holds an `EEG` struct:
- `.data` — 178 channels × 2049 freq bins (PSD, 0–250 Hz, freqstep 0.122 Hz)
- `.bands` — **14 bands precomputed per channel**: each `.label` (low-delta, delta, theta, alpha, sigma,
  beta, gamma, …), `.type` (`absolute`/`normalized`), `.freqrange` `[lo hi]`, `.data` = **178×1 band power**
- `.chanlocs` — 178 channels (E1…E178) with X/Y/Z; `.subject`, `.run` set; `.condition`/`.group` empty
  (parse from filename/folder)

Implications for the new engines — **both built & validated on the real data**:
1. **178-channel neighbour matrix** — `dependencies/build_178_neighbors.m` builds `NeighborMatrix_178.mat`
   from chanlocs X/Y/Z (distance threshold; mean degree ~6, connected). The 178 path in `core_snpm_glm.m`
   / `core_snpm_lmm.m` now loads + remaps it (`remap_neighbors`). Also fixes the legacy 178 path.
2. **Spectral assembly loader** — `load_spectral_dataset.m`: walks the tree, parses subject/condition/stage,
   pulls a chosen band×type 178-vector from `EEG.bands`, assembles a table (`meta_cols` =
   Subject/condition/stage/run + 178 channel columns) ready for any preset. Band power is precomputed.
3. **Global spectrum comparison** — `plot_global_spectrum.m`: mean PSD across channels per subject,
   averaged within group (condition/stage), overlaid log-power vs frequency (the reference figure).

Validated end-to-end: `anova1` (condition a vs b @N2), `rmanova` (N2/N3/REM), global spectrum — all on the
178 montage, outputs written.

**Usage caveat:** condition-a / condition-b are within-subject (same subjects). Compare them with
`rmanova` (condition as the within factor) or the legacy `pairedT`, **not** `anova1` (which assumes
independent groups). Use `anova1`/`ancova` for between-subject groups (e.g. patient vs control).

## Deferred / open items

1. **GUI** — presets are headless/NCI only; not yet added to `SnPMAnalysisGui.m`
   (`ComparisonDropDown.Items` + per-preset column-mapping panel).
2. **TurtleWave importer** — separate task once a sample event file is available; produce both an
   aggregated-per-channel table and an event-level table.
3. **Loader defaults** — wire `meta_cols` / column-mapping defaults to the user's real file column names
   once a sample wide/long file is shared.

## How to run (headless)

```matlab
parpool;                      % LMM only (GLM needs no pool)
[results, text] = core_snpm_analysis(params);
```
`params.comparison` selects the analysis; see `core_snpm_lmm.m` / `core_snpm_glm.m` headers for the full
parameter list per preset.
