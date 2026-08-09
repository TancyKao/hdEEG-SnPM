# Documentation

**Version 1.0.0** · MATLAB toolbox for permutation-based (TFCE and cluster) statistics on
high-density EEG topographies · [source on GitHub](https://github.com/TancyKao/hdEEG-SnPM)

This documentation is organised with the [Diátaxis](https://diataxis.fr) framework: four kinds
of material, each serving a different need. Pick by what you're trying to do right now.

| I want to… | Go to | What it is |
|---|---|---|
| **Learn the tool** by doing a first run | **[Tutorials](tutorials/README.md)** | A guided lesson on the shipped example data, end to end. |
| **Get a specific task done** | **[How-to guides](how-to/README.md)** | Short, direct recipes for a known goal. |
| **Look up a fact** (inputs, outputs, keys) | **[Reference](reference/ANALYSIS_CATALOG.md)** | Neutral, complete description of the machinery. |
| **Understand why** it works this way | **[Explanation](explanation/README.md)** | Design rationale and statistical background. |

New here? Start with **[Getting started: your first SnPM analysis](tutorials/getting-started.md)**.
It runs on the ready-made data in
[**`example_data/`**](https://github.com/TancyKao/hdEEG-SnPM/blob/main/example_data/README.md) —
one file (or file pair) per analysis, each with a 17-channel cluster planted near E129, already
in the repository.

**Arrived from TurtleWave?** hdEEG-SnPM is the group-statistics half of a two-program pipeline
whose detection half is [TurtleWave-hdEEG](https://github.com/TancyKao/TurtleWave-hdEEG). Go
to [Take TurtleWave detections into a group analysis](how-to/analyse-turtlewave-event-data.md)
— and if you are analysing coupling phase, read the pre-v4.0 inversion warning there before you
run anything. TurtleWave is a supported source, not a requirement: this toolbox takes a
subjects × channels table from anywhere and needs no Python.

## The four quadrants

### Tutorials — learning-oriented
Take you by the hand through a first, working analysis so you build confidence.
- [Getting started: your first SnPM analysis](tutorials/getting-started.md)

### How-to guides — goal-oriented
Assume you know what you want; get you there.
- [Run an analysis in the GUI](how-to/run-an-analysis-in-the-gui.md)
- [Take TurtleWave detections into a group analysis](how-to/analyse-turtlewave-event-data.md)
- [Run a circular (phase) analysis](how-to/run-a-circular-phase-analysis.md)
- [Run a source-level (2447-voxel) analysis](how-to/run-a-source-level-analysis.md)
- Run analyses headless (HPC / compute cluster) → [`scripts/README_scripts.md`](https://github.com/TancyKao/hdEEG-SnPM/blob/main/scripts/README_scripts.md)

### Reference — information-oriented
The source of truth; consult, don't read cover to cover.
- [Analysis catalog](reference/ANALYSIS_CATALOG.md) — the backbone (one row per analysis)
- [Reference index](reference/README.md) — comparison keys + colocated example-data/script references
- Example data → [`example_data/README.md`](https://github.com/TancyKao/hdEEG-SnPM/blob/main/example_data/README.md) (in the repository; one file or file pair per analysis)

### Explanation — understanding-oriented
Background reading; the *why* behind the design.
- [Detection and statistics are two tools](explanation/detection-and-statistics-two-tools.md)
- [Choosing an analysis: design, missing data, and cost](explanation/choosing-an-analysis.md)
- [Two-way mixed ANOVA — how the interaction is computed](explanation/mixed2way-interaction.md)
- [About circular statistics for phase](explanation/circular-statistics-for-phase.md)
- [Missing data and the complete-column rule](explanation/missing-data-and-excluded-channels.md)
- [Interpreting cluster results](explanation/interpreting-cluster-results.md)
- [Why this tool exists alongside a PALM-based pipeline](explanation/eeg-processor-vs-snpm.md)

## Archive

Planning notes, the v2.0 changelog, and design-assistant prompts are kept in
[`archive/`](archive/README.md) for provenance. They are **frozen history**, not maintained
user documentation — don't follow them as instructions. Releases are tagged in the repository;
see the [release history](https://github.com/TancyKao/hdEEG-SnPM/releases).

## Two READMEs live outside `docs/`

By design, two reference/how-to documents sit next to the files they describe so they stay in
sync with them. Neither is part of this published site — only `docs/` is published — so the
links below go to the files on GitHub:

- [`scripts/README_scripts.md`](https://github.com/TancyKao/hdEEG-SnPM/blob/main/scripts/README_scripts.md)
  — the headless/HPC runners (how-to + reference).
- [`example_data/README.md`](https://github.com/TancyKao/hdEEG-SnPM/blob/main/example_data/README.md)
  — the example data (reference): which file drives which analysis, the roles to set, and what
  each one plants.

The reference section below is the source of truth for how each analysis is wired.
