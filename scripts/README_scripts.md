# Headless / NCI scripts — run analyses without the GUI

These templates run the same analyses as **hdEEG-SnPM Toolbox** (`SnPMAnalysisGui.m`)
from the command line, so jobs can be submitted to NCI/HPC. Edit the `CONFIG`
block at the top of each script, then run headless:

```bash
matlab -batch "run('scripts/run_glm_analysis.m')"
```

`matlab -batch` needs no display, exits on completion, and returns a non-zero
status on error — suitable for a PBS/Slurm job. All scripts call
`addpath(genpath(ROOT))` themselves, so only `ROOT` must be correct.

## Which script

| Script | Analysis | Engine |
|---|---|---|
| `run_glm_analysis.m`   | anova1 / ancova / regression / rmanova / mixed2way, or legacy pairedT / unpairedT / onesampleT / correlation | `core_snpm_analysis` → `core_snpm_glm` |
| `run_lmm_analysis.m`   | per-channel linear mixed model (Stephan 2021) | `core_snpm_analysis` → `core_snpm_lmm` |
| `run_spectral_report.m`| two-condition spectral report (periodograms + topo + HTML) | `export_report` |
| `run_spectral_analysis.m`| multi-folder spectral SnPM (folder = factor level; unpairedT/pairedT/onesampleT/anova1/rmanova; band×stage sweep) | `spectral_to_snpm_params` → `core_snpm_analysis` |
| `run_event_report.m`   | event-level (spindle/SO) report (group t-map + HTML) | `export_event_report` |
| `export_parity_fixtures.m` | not an analysis: writes MATLAB reference results plus the exact permutation indices, for checking another implementation against MATLAB (see below) | `snpm_perm_tap` + `core_snpm_analysis` |

There is **no `run_circular_analysis.m` template**. The circular (phase/angle) tier is driven by
building a `params` struct directly — see
[Run a circular (phase) analysis](../docs/how-to/run-a-circular-phase-analysis.md), which carries
the full struct, the required `circ_units` / `circ_convention` / precision-file inputs, and the
detector zero-convention table. Input files for all four circular configurations ship in
[`example_data/`](../example_data/README.md) — the two `circ_phase_group{A,B}.csv` angle files
with their `circ_phase_counts{A,B}.csv` precision files (Hotelling and Watson U²),
`circ_conc_group{A,B}.csv` + counts (a concentration difference at the same mean direction), and
`circ_anglinear_angles.csv` + `circ_anglinear_measure.csv` (circular–linear).

## Data prep

To try a template before you have your own data, point its `CONFIG` at
[`example_data/`](../example_data/README.md) at the repo root — one ready-made CSV (or file pair)
per analysis, each with a 17-channel cluster planted near E129 on the EGI 256 montage, so you can
confirm the run recovers an effect it was built to find. `run_glm_analysis.m` covers
`glm_anova1.csv`, `glm_ancova.csv`, `glm_ancova3.csv`, `glm_regression.csv`, `glm_rmanova.csv`,
`glm_mixed2way.csv` and the legacy pairs; `run_lmm_analysis.m` covers `lmm_long.csv`.

- **Spectral** (`run_spectral_report.m`): point `DATA` at the folder of per-subject
  spectral `.mat` files; the script scans them itself.
- **GLM / event** (`run_glm_analysis.m`, `run_event_report.m`): need a subject×channel
  matrix + a `subjects.csv` metadata table. For TurtleWave events, generate them
  with the Python bridge first (no MATLAB, stdlib only):

  ```bash
  # DATA_dir and OUT_dir are POSITIONAL. --group labels the cohort in subjects.csv;
  # pass it explicitly whenever you have more than one group.
  python3 db_to_group_table.py <DATA_dir> <OUT_dir> --group <label>
  # optional: --stages NREM2,NREM3 (default) | NREM3 ...   stages analysed
  #           --bands spindle_9-12,spindle_12-15,sw_0.5-1.25  subset of event bands
  #           --allow-hypnogram-denominator   accept pre-4.3 databases (biased density)
  # writes: subjects.csv, channels.csv, recordings.csv, eventStat_<param>_<event>.csv
  ```

  Density is events per minute of TurtleWave's stored artefact-free analysed time
  (TurtleWave 4.3+), over the stages the detection run actually searched, out of
  `--stages`; a band whose subjects would be divided over different stages is
  refused. An older database stops the run with an error unless
  `--allow-hypnogram-denominator` is given, which divides it by hypnogram minutes
  instead, with a warning per recording; a run mixing the two denominators is
  refused even then. `recordings.csv` (one row per subject and band) records which
  denominator and stages each subject got; do not merge tables from separate runs
  across denominator types or stage sets. A band never detected in a subject is
  written as a blank row with a warning. The bridge reads
  one detection run per band and stops with an error if two runs cover the same
  stage, if a band's events all lie outside `--stages`, or if a band with no event
  has more than one method in TurtleWave's processing record; it sets interpolated
  channels to missing. Rules and options:
  `docs/how-to/analyse-turtlewave-event-data.md`. `run_event_report.m` expects all
  three bands, so do not combine it with `--bands`.

  `run_glm_analysis.m` joins `matrix_csv` + `subjects.csv` on `Subject` and runs
  one matrix; `run_event_report.m` loops all parameters × event bands itself.

## Spectral folder as a data source (`run_spectral_analysis.m`)

EEG_processor emits one BIDS spectral-power folder per cohort/condition
(`<folder>/sub-*_..._desc-<stage>_powerspect.mat`, band power in `EEG.features`).
The model: **each folder is one level of the design factor.** List the folders +
labels, fix a `band` + power `type`, and the analysis compares them per stage; the
178 channels are the response map. Same subject across two folders (matched by the
`sub-XX` token) = within-subject; distinct subjects per folder = between-group.

CONFIG: `FOLDERS = {pathA, pathB, ...}` with parallel `LABELS = {'a','b',...}`.
`spectral_to_snpm_params.m` slices each cell (`level` = folder label) → temp CSV(s) →
`core_snpm_analysis` (full `.mat`/`.xlsx`/HTML per cell). A `sweep_bands × sweep_stages`
grid loops the cells and writes `SWEEP_grid_*.csv` (`nUncorr`/`nTFCE`/`nCluster`/`minClusterP`).

| comparison | #folders | factor (= folder label) |
|---|---|---|
| `unpairedT` | 2 | two independent groups (between) — **subject counts may differ between the two folders** |
| `pairedT` / `onesampleT` | 2 | two conditions, same subjects (within) |
| `anova1` | ≥2 | groups (between) |
| `rmanova` | ≥2 | conditions, same subjects (within) |

Two cohort folders with unequal subject counts on the `unpairedT` route now run end to end; before
the 2026-08 fix that combination threw a dimension error before reaching a statistic.

**CSV-only** (need per-subject covariates/predictor or trial-level data, so not this
folder source): `ancova`, `regression`, `correlationP/S`, `mixed2way`, `mixedmodel` —
use `run_glm_analysis.m` / `run_lmm_analysis.m`.

**Faceted dashboard:** set `MAKE_DASHBOARD = true` (default on) and `run_spectral_analysis.m`
also writes ONE self-contained dashboard HTML — `OUT/dashboard/sleep_eeg_report_filled.html` —
with **Absolute + Relative power × all stages × Uncorrected/TFCE/Cluster** toggles (via
`export_report`, which now accepts `comparison`/`folders`/`labels`). Supported for the **2-level
contrasts only** (pairedT/onesampleT/unpairedT); anova1/rmanova still get the per-cell reports +
grid CSV (faceted dashboard for omnibus designs is a follow-up).

Notes: `type` is the band **power** type (`absolute`/`normalized`); `datatype` is a
further transform (auto: absolute→`logscale`, normalized→none). `level_A`/`level_B`
optionally pick/order the two levels for the 2-folder tests. Filenames must tokenize
cleanly — a stray `..._condition-b-psg_...` (missing `_task-` separator) mis-parses and
drops that subject.

## Reference fixtures for another implementation (`export_parity_fixtures.m`)

This script does not answer a research question. It produces reference files for checking a
separate implementation of the same analyses (for example a Python one) against MATLAB. Two
programs that use the same seed still draw different permutations, so the script saves the
permutation indices MATLAB actually drew, and the other program replays them.

```bash
matlab -batch "run('scripts/export_parity_fixtures.m')"
```

For each case the script seeds the generator (`rng(SEED,'twister')`), switches on the recording
hook `dependencies/snpm_perm_tap.m`, runs `core_snpm_analysis`, and saves one
`<case_name>.mat` (format `-v7`, readable with `scipy.io.loadmat(path, simplify_cells=True)`)
holding the inputs, the montage, the permutation indices and the full results; for the t-test and
correlation analyses it also holds the matrices handed to the statistic engines, the
cluster-size null and the channel set of the global test. It then writes a
`README.md` in the output folder with one row per fixture and a description of every field.
Cases cover `pairedT`, `unpairedT`, `onesampleT`, `correlationP` and `correlationS` (with and
without covariates), `anova1`, `ancova`, `regression`, `rmanova` and `mixed2way`, using the
files in [`example_data/`](../example_data/README.md) plus a few inputs the script builds from
them. The circular and mixed-model tiers are not exported yet.

| CONFIG field | Default | Meaning |
|---|---|---|
| `OUT_DIR` | `test_fixtures/parity/` under the repo root | where the fixtures and their `README.md` go (this folder is not tracked by git) |
| `SEED` | `20261013` | generator seed set immediately before every case |
| `PERMUTATIONS` | `500` | permutations requested per permutation loop |
| `INPUT_SEED` | `4242` | seed for the inputs the script builds itself |
| `ONLY_CASES` | `{}` | `{}` runs every case; a list of case names re-exports only those (the `README.md` is then left as it was) |
| `RUN_DIR` | a folder under `tempdir` | where each run writes its own `.mat`/`.xlsx`/HTML/PNG output; keep it outside the repo, because `core_snpm_analysis` adds everything under the current folder to the MATLAB path |
| `RESUME` | `false` | `true` keeps any fixture that already exists with the same `SEED` and `PERMUTATIONS` and runs only the missing cases, for use after an interrupted export. It checks nothing else: after editing the script or the toolbox, export with `false` |
| `KEEP_RUN_OUTPUTS` | `false` | `false` deletes each run's own outputs once its fixture is saved |

For the General Linear Model presets the pipeline does not keep its design matrices, so the
script rebuilds every permutation loop (design, contrast, rows) with the same functions the
pipeline calls, replays it from the recorded indices, and refuses to save a fixture whose replay
does not match the pipeline's own statistic, uncorrected p, enhanced map and sorted enhanced-map
null exactly. The cluster-mass null is not recomputed for the fixture: it is recorded inside the
pipeline where it is computed, and the replay's own value is used only as a cross-check that must
equal the recorded one. GLM fixtures carry `glm.scheme = '2026-10-14'`, which marks the method in
which each pairwise post-hoc comparison is its own two-level test on the rows of its two levels
and the two-way mixed design uses one channel set for every effect. The recording hook is off
unless this script switches it on, and it changes no result either way.

## Column roles (GLM)

Set only the roles the chosen `comparison` needs (leave `''`):

| comparison | needs |
|---|---|
| `anova1`     | `group_col` (≥2 levels) |
| `ancova`     | `group_col` + `covariate_cols` |
| `regression` | `predictor_col` (+ optional `covariate_cols`) |
| `rmanova`    | `subject_col` + `condition_col` |
| `mixed2way`  | `subject_col` + `group_col` + `condition_col` |

`mixed2way` needs a balanced table: exactly one row per condition for every subject, and one
group per subject. Otherwise the run stops with `core_snpm_glm:mixed2wayUnbalanced`, naming the
subjects; remove or complete them, or use `run_lmm_analysis.m`. With three or more levels,
`anova1`, `ancova`, `rmanova` and `mixed2way` also write pairwise post-hoc maps, each a two-level
test on the rows of its two levels. The rules are in
[`docs/reference/ANALYSIS_CATALOG.md`](../docs/reference/ANALYSIS_CATALOG.md).

## NCI notes

- One analysis per job; fan out parameters/bands as an array job (each a separate
  `matlab -batch` call with a different `CONFIG`), or copy a script per cell.
- `permutations=1000` is the default; raise for final inference. Permutation p is
  `(b+1)/(N+1)` across all tiers, so the smallest reportable p is `~1/(N+1)` — a real effect
  can't clear α unless `permutations` is large enough (≥5000–10000). The legacy
  `unpairedT`/correlation paths switch to an **exact** test when `permutations ≥` the number of
  distinct labelings (`nchoosek(nSubj,nGrp)` / `nSubj!`); very small samples then hit a coarse
  p-grid (e.g. unpaired nGrp=3, nSubj=6 → 20 perms → min p ≈ 0.048).
- `channels` now selects a **recording system** (label-based, count-agnostic): `'egi'` =
  EGI 256 / HydroCel (egi257 178-scalp; legacy `'164 channels'`/`'178 channels'` still resolve
  here) and `'compu'` = Compumedics 257 / Neuvo (compu257 249-scalp). The data file's column
  names are matched against the system's channel labels (`hdeeg_scalpchannels`), so columns must
  use real labels (`E1..E256`/`Cz` for EGI, `Fp1`,`Fpz`,… for Compumedics). Add a system in
  `dependencies/snpm_montage_registry.m`.
- Outputs (Excel + topoplots + HTML) land in `output_dir` / `OUT`.
