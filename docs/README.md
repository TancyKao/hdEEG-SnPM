# hdEEG-SnPM Toolbox — documentation

This documentation is organised with the [Diátaxis](https://diataxis.fr) framework: four kinds
of material, each serving a different need. Pick by what you're trying to do right now.

| I want to… | Go to | What it is |
|---|---|---|
| **Learn the tool** by doing a first run | **[Tutorials](tutorials/)** | A guided lesson on synthetic data, end to end. |
| **Get a specific task done** | **[How-to guides](how-to/README.md)** | Short, direct recipes for a known goal. |
| **Look up a fact** (inputs, outputs, keys) | **[Reference](reference/ANALYSIS_CATALOG.md)** | Neutral, complete description of the machinery. |
| **Understand why** it works this way | **[Explanation](explanation/README.md)** | Design rationale and statistical background. |

New here? Start with **[Getting started: your first SnPM analysis](tutorials/getting-started.md)**.

## The four quadrants

### Tutorials — learning-oriented
Take you by the hand through a first, working analysis so you build confidence.
- [Getting started: your first SnPM analysis](tutorials/getting-started.md)

### How-to guides — goal-oriented
Assume you know what you want; get you there.
- [Run an analysis in the GUI](how-to/run-an-analysis-in-the-gui.md)
- Run analyses headless (NCI/HPC) → [`scripts/README_scripts.md`](../scripts/README_scripts.md)

### Reference — information-oriented
The source of truth; consult, don't read cover to cover.
- [Analysis catalog](reference/ANALYSIS_CATALOG.md) — the backbone (one row per analysis)
- [Reference index](reference/README.md) — comparison keys + colocated fixture/script references
- Synthetic fixtures → [`test_data/synthetic_gui/README.md`](../test_data/synthetic_gui/README.md)

### Explanation — understanding-oriented
Background reading; the *why* behind the design.
- [Choosing an analysis: design, missing data, and cost](explanation/choosing-an-analysis.md)
- [Two-way mixed ANOVA — how the interaction is computed](explanation/mixed2way-interaction.md)
- [EEG_Processor group analysis vs. SnPM_2025](explanation/eeg-processor-vs-snpm.md)

## Archive

Planning notes, session recaps, the v2.0 changelog, and design-assistant prompts are kept in
[`archive/`](archive/README.md) for provenance. They are **frozen history**, not maintained
user documentation — don't follow them as instructions.

## Two READMEs live outside `docs/`

By design, two reference/how-to documents sit next to the files they describe so they stay in
sync with them:

- [`scripts/README_scripts.md`](../scripts/README_scripts.md) — the headless/NCI runners (how-to + reference).
- [`test_data/synthetic_gui/README.md`](../test_data/synthetic_gui/README.md) — the synthetic fixtures (reference).

The project architecture and orchestration model live in the root `CLAUDE.md`.
