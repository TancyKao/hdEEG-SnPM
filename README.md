# hdEEG-SnPM Toolbox

MATLAB tool for **Statistical non-Parametric Mapping (SnPM)** of high-density EEG
topographies — permutation-based TFCE and cluster statistics across scalp channels with
channel-adjacency (neighbour) structure. Group comparisons, correlations, GLM presets
(ANOVA/ANCOVA/regression/RM-ANOVA/two-way mixed), per-channel LMM, and circular variants.
Produces topoplots, an Excel table, a `.mat` results struct, and HTML reports.

Theory: Maris & Oostenveld (2007), Nichols & Holmes (2001) — PDFs in `Ref/`.

## Run it

No build step; requires MATLAB (App Designer + Statistics Toolbox). **Run from the repo
root** (`core_snpm_analysis` does `addpath(genpath(pwd))` and expects `eeglab2022.1/` and
`dependencies/` to resolve from the current directory).

- **GUI:** `cd` to the repo root, then run `SnPMAnalysisGui`. Pick an analysis, load
  file(s) or a spectral folder, set parameters, Run.
- **Headless (NCI/HPC):** edit the `CONFIG` block of a `scripts/run_*.m` template and run
  `matlab -batch "run('scripts/<name>.m')"`. See `scripts/README_scripts.md`.
- **Verify everything:** `matlab -batch "test_all"` → one PASS/FAIL/SKIP summary
  (checkcode + GUI construct + all test suites + per-analysis outputs + spectral smoke +
  edge cases). Individual suites: `test_glm_snpm`, `test_lmm_snpm`, `test_legacy_snpm`.
  Per-analysis import/output reference: `docs/ANALYSIS_CATALOG.md`.

## Layout

| Path | What |
|---|---|
| `SnPMAnalysisGui.m` | the App Designer GUI (thin front-end) |
| `core_snpm_analysis.m` | orchestrator; routes to `core_snpm_glm.m` / `core_snpm_lmm.m` |
| `export_*.m`, `plot_global_spectrum.m`, `load_spectral_dataset.m` | reports + I/O |
| `gen_synthetic_*.m`, `test_*.m` | synthetic-data generators + test suites |
| `dependencies/` | stat/plot/channel engines + montage assets (`.mat`/`.sfp`/`.csv`) |
| `templates/` | report HTML templates (anova / t / spectral / LMM) |
| `scripts/` | headless runners + `README_scripts.md` + Python bridges |
| `docs/` | `PLAN.md` (design + analysis tables), `RECAP.md`, `MODIFICATIONS_SUMMARY.md`, `DESIGN_PROMPT_*.md` |
| `test_data/` | inputs + synthetic fixtures (run outputs are git-ignored / regenerable) |
| `Ref/` | method papers |
| `eeglab2022.1/` | vendored EEGLAB (third-party; don't modify) |

See `CLAUDE.md` for the full architecture and conventions, and `docs/PLAN.md` for the
analysis design and which test answers which research question.
