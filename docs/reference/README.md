# Reference

*Information-oriented.* Look things up here; these pages describe the machinery neutrally and
completely. They do not teach (see [tutorials](../tutorials/)) or walk you through a task (see
[how-to guides](../how-to/README.md)).

## Pages

- **[Analysis catalog](ANALYSIS_CATALOG.md)** — the source of truth. One row per analysis:
  research question, what to import, required columns/roles, the statistic and its reference,
  its test fixture, and what it outputs. The backbone of this quadrant.

## Reference that lives next to the code it describes

Two reference documents are colocated with the assets they describe (kept there on purpose so
they stay in sync with the files):

- **Synthetic fixtures** — [`test_data/synthetic_gui/README.md`](../../test_data/synthetic_gui/README.md).
  Maps each synthetic file → analysis → GUI settings. Every file plants a 17-channel neighbour
  cluster near E129 on the EGI 256 montage.
- **Headless / NCI scripts** — [`scripts/README_scripts.md`](../../scripts/README_scripts.md).
  Which `run_*.m` template runs which analysis/engine, the column-role table per preset, the
  spectral-folder-as-factor model, and NCI notes.

## Comparison keys (`params.comparison`)

The internal analysis keys used by `core_snpm_analysis` and the GUI's `ComparisonDropDown`
(`ItemsData`):

| Key | Analysis | Engine |
|---|---|---|
| `onesampleT`, `pairedT`, `unpairedT` | one-sample / paired / unpaired t | legacy `compstring` switch |
| `correlationP`, `correlationS` | Pearson / Spearman correlation | legacy `compstring` switch |
| `circ_wheeler_watson_Test`, `circ_WatsonsU2Test` | circular tests | legacy `_circ` switch (CircStat toolbox required) |
| `anova1`, `ancova`, `regression`, `rmanova`, `mixed2way` | GLM presets | `core_snpm_glm` |
| `mixedmodel` | per-channel linear mixed model | `core_snpm_lmm` |

`params.tail` ∈ `both` / `left` / `right` for the legacy t-tests and correlations. The GLM
and LMM presets pick statistic, contrast and permutation scheme automatically from the design.
