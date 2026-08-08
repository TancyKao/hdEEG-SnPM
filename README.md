# hdEEG-SnPM Toolbox

**Version 1.0.0** (tag `v1.0.0`) · MIT licence · [Changelog](CHANGELOG.md) ·
[Documentation site](https://tancykao.github.io/hdEEG-SnPM/)

MATLAB tool for **Statistical non-Parametric Mapping (SnPM)** of high-density EEG
topographies — permutation-based TFCE and cluster statistics across scalp channels with
channel-adjacency (neighbour) structure. Group comparisons, correlations, GLM presets
(ANOVA/ANCOVA/regression/RM-ANOVA/two-way mixed), per-channel LMM, and a circular (phase/angle)
tier — Hotelling T², Watson's U², circular–linear correlation.
Produces topoplots, an Excel table, a `.mat` results struct, and HTML reports.

It was developed for sleep hd-EEG research (spectral power, spindle/slow-wave measures,
condition and group contrasts across the scalp), but nothing in the statistics is
sleep-specific: any subject × channel measure with a montage the toolbox knows will run.
The documentation lives in **[`docs/`](docs/README.md)** — start there.

**What it does, in one line:** you give it a table of one number per subject per channel, tell
it which columns are the design (group, condition, subject, covariates), and it returns which
scalp regions differ, corrected for testing every channel at once by a permutation null rather
than a parametric assumption.

## Where this fits — the other half of the pipeline

hdEEG-SnPM is the **group-statistics** half of a two-program sleep hd-EEG pipeline. The
**detection** half is [**TurtleWave-hdEEG**](https://github.com/TancyKao/TurtleWave-hdEEG)
(Python ≥ 3.10, MIT): a Wonambi-based detector for spindles, slow waves and K-complexes with
phase–amplitude coupling, parallel/HPC batch processing, and two PyQt5 interfaces
(`turtlewave_gui` for detection, `eeg_review_gui` for review). It answers *where are the events
in this recording*; hdEEG-SnPM answers *do these groups differ across the scalp, corrected for
testing every channel*.

The join between them is a file, not a library: `db_to_group_table.py` in this repository reads
TurtleWave's Wonambi SQLite database and staged hypnogram and writes the wide subjects × channel
tables plus a `subjects.csv` that the analyses consume. Aggregated per-channel measures (spindle
density, amplitude, duration) feed the t-tests, correlations and GLM presets; event-level rows
feed the per-channel mixed model; coupling **phase** feeds the circular tier.

**TurtleWave is a supported source, not a dependency.** This toolbox takes a subjects × channels
table from wherever you like — a spectral folder, another detector, a spreadsheet — and needs no
Python installed.

> **If you analyse coupling phase from TurtleWave, read this first.** Exports predating
> **v4.0** report preferred phase **exactly 180° inverted** (fixed upstream in
> [commit `d341690`](https://github.com/TancyKao/TurtleWave-hdEEG/commit/d3416909571b7f6fa39c624b72136acf8fc566ef),
> 16 July 2026). An inverted angle looks entirely plausible — tight distribution, significant
> cluster — and says events couple to the down-state when they couple to the up-state. Either
> re-export from v4.0+, or declare `circ_convention = 'turtlewave_pre_v4'`. Details:
> [Take TurtleWave detections into a group analysis](docs/how-to/analyse-turtlewave-event-data.md).

Theory:

- Maris E, Oostenveld R (2007). *Nonparametric statistical testing of EEG- and MEG-data.*
  Journal of Neuroscience Methods, 164(1), 177–190. doi:10.1016/j.jneumeth.2007.03.024
- Nichols TE, Holmes AP (2001). *Nonparametric permutation tests for functional neuroimaging:
  a primer with examples.* Human Brain Mapping, 15(1), 1–25. doi:10.1002/hbm.1058

Freedman–Lane nuisance permutation follows Winkler AM et al. (2014), *Permutation inference for
the general linear model*, NeuroImage 92, 381–397 (doi:10.1016/j.neuroimage.2014.01.060); the
per-channel mixed-model path follows Stephan AM et al. (2021), *Conscious experiences and
high-density EEG patterns predicting subjective sleep depth*, Current Biology 31(24), 5487–5500.

## Requirements

- **MATLAB** — developed and tested on R2025a. The GUI uses App Designer
  (`matlab.apps.AppBase`), so a MATLAB new enough to run App Designer apps is required; the
  headless scripts work without a display.
- **Statistics and Machine Learning Toolbox** — required, not optional (`fitlm`, `fitlme`,
  `anovan`, `tinv`, `corr`, …).
- **EEGLAB — not included in this repository.** See below.

### EEGLAB is not bundled — install it yourself

The toolbox calls EEGLAB's `topoplot` to render every topographic map. **EEGLAB is not in this
repository** (it is ~116 MB of third-party code and is git-ignored), so a fresh clone has code
that calls EEGLAB functions but no EEGLAB. Statistics will compute, but any plotting step will
fail with an undefined-function error until you install it.

Download EEGLAB from <https://sccn.ucsd.edu/eeglab/> and either:

- unpack it at the repo root as `eeglab2022.1/` (the path `core_snpm_analysis` looks for by
  default), or
- point at your own copy with `params.eeglab_path = '/path/to/eeglab'`.

Version 2022.1 is what this was developed against; other recent versions should work, since
only `topoplot` and its helpers are used.

## Run it

No build step. **Run from the repo root** — `core_snpm_analysis` does `addpath(genpath(pwd))`
and resolves `eeglab2022.1/` and `dependencies/` relative to the current directory. Override
with `params.snpm_path` / `params.eeglab_path` if you must run from elsewhere.

- **GUI:** `cd` to the repo root, then run `SnPMAnalysisGui`. Pick an analysis, load
  file(s) or a spectral folder, set parameters, Run.
- **Headless (HPC/cluster):** edit the `CONFIG` block of a `scripts/run_*.m` template and run
  `matlab -batch "run('scripts/<name>.m')"`. See `scripts/README_scripts.md`.
- **Verify everything:** `matlab -batch "test_all"` → one PASS/FAIL/SKIP summary
  (checkcode + GUI construct + all test suites + per-analysis outputs + spectral smoke +
  edge cases). Individual suites: `test_glm_snpm`, `test_lmm_snpm`, `test_legacy_snpm`,
  `test_circ_snpm`, `test_circ_stats`, `test_cluster_enhancement_identity`, `test_source_snpm`.
  Per-analysis import/output reference: `docs/reference/ANALYSIS_CATALOG.md`.

### First run — generate the example data

No real data ships with this repository. Generate the synthetic fixtures instead — each one has
a known 17-channel cluster planted near channel E129, so you can check the tool recovers an
effect it was told to find:

```matlab
gen_synthetic_testdata      % writes test_data/synthetic_gui/ + a README mapping
                            % each file -> analysis -> GUI settings
```

Then follow [the getting-started tutorial](docs/tutorials/getting-started.md).

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
| `docs/` | Diátaxis docs, also published as a website (`README.md` index): `tutorials/`, `how-to/`, `reference/` (`ANALYSIS_CATALOG.md`), `explanation/`, `archive/` (frozen history: `PLAN.md`, `MODIFICATIONS_SUMMARY.md`, `DESIGN_PROMPT_*.md`) |
| `test_data/` | *not in the repository* — generated inputs and run outputs (git-ignored, regenerable via `gen_synthetic_*.m`) |
| `eeglab2022.1/` | *not in the repository* — install EEGLAB here yourself (see above) |

See `CLAUDE.md` for the full architecture and conventions. Start the docs at
[`docs/README.md`](docs/README.md); the analysis design and research-question rationale are in
`docs/archive/PLAN.md` (history) and `docs/explanation/choosing-an-analysis.md`.
Release-by-release changes are in [`CHANGELOG.md`](CHANGELOG.md).

## License

This toolbox is released under the **MIT License** — see [`LICENSE`](LICENSE).
Copyright (c) 2026 Tancy Kao.

Two third-party components are governed by their own terms, not by the MIT licence above:

- **CircStat** (Berens P, 2009) — a partial, unmodified copy is vendored under
  `dependencies/circstat/` and is distributed under the **BSD 2-clause licence**; its
  `license.txt` is included verbatim in that folder, as the licence requires. See
  `dependencies/circstat/README_VENDORED.md` for what is vendored and why (including which
  functions are validation-only and must not be used for inference). If you publish work using
  it, the upstream readme asks you to cite: Berens P (2009), *CircStat: A MATLAB Toolbox for
  Circular Statistics*, Journal of Statistical Software 31(10).
- **EEGLAB** — **not included in this repository** and not covered by this licence. You must
  obtain it separately (see above) and comply with the licence of the copy you install. The
  EEGLAB 2022.1 core (`eeglab.m` and `functions/`) is BSD 2-clause per its own
  `eeglablicense.txt`, but EEGLAB plugins are released under other licences, some of them GPL —
  check the terms of what you actually install before redistributing anything that includes it.
