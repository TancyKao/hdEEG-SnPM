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

## Data prep

- **Spectral** (`run_spectral_report.m`): point `DATA` at the folder of per-subject
  spectral `.mat` files; the script scans them itself.
- **GLM / event** (`run_glm_analysis.m`, `run_event_report.m`): need a subject×channel
  matrix + a `subjects.csv` metadata table. For TurtleWave events, generate them
  with the Python bridge first (no MATLAB, stdlib only):

  ```bash
  python3 db_to_group_table.py --root <cohort_dir> --out <DATA>
  # writes: subjects.csv, channels.csv, eventStat_<param>_<event>.csv
  ```

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
| `unpairedT` | 2 | two independent groups (between) |
| `pairedT` / `onesampleT` | 2 | two conditions, same subjects (within) |
| `anova1` | ≥2 | groups (between) |
| `rmanova` | ≥2 | conditions, same subjects (within) |

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

## Column roles (GLM)

Set only the roles the chosen `comparison` needs (leave `''`):

| comparison | needs |
|---|---|
| `anova1`     | `group_col` (≥2 levels) |
| `ancova`     | `group_col` + `covariate_cols` |
| `regression` | `predictor_col` (+ optional `covariate_cols`) |
| `rmanova`    | `subject_col` + `condition_col` |
| `mixed2way`  | `subject_col` + `group_col` + `condition_col` |

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
