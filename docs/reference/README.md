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
| `circ_phase_group` | phase, 2 groups — covariate-adjusted Hotelling T² | `core_snpm_circ` |
| `circ_phase_group_u2` | phase, 2 groups — Watson's U² | `core_snpm_circ` |
| `circ_corrAngLinear` | circular–linear correlation | `core_snpm_circ` |
| `anova1`, `ancova`, `regression`, `rmanova`, `mixed2way` | GLM presets | `core_snpm_glm` |
| `mixedmodel` | per-channel linear mixed model | `core_snpm_lmm` |

### `snpm_perm_correction` optional arguments

The shared TFCE + cluster-mass driver (`dependencies/snpm_perm_correction.m`) takes two optional
trailing positional arguments, in this order:

```
snpm_perm_correction(real_stat, real_p, perm_stat_fn, neighbors, E, H, alpha, ...
                     permutations, contrast_type, evaluable, dh)
```

- **`evaluable`** (10th) — `1 × nCh` logical; channels that are `false` are NaN in the observed
  map *and* in every permuted map. Omit or pass `[]` for all channels evaluable.
- **`dh`** (11th) — TFCE integration step forwarded to `ClusterEnhancement`. Omit or pass `[]`
  to let `ClusterEnhancement` apply its own default of 0.1; when `dh` is absent the fifth
  argument is genuinely not passed, so the default path is byte-identical to the pre-2026-08
  behaviour. The same `dh` is applied to the observed map and to every permuted map. A
  non-positive or non-finite `dh` errors with `snpm_perm_correction:badDh`.

`dh` is positioned **after** `evaluable`, so **to set `dh` while leaving all channels evaluable,
pass `[]` for `evaluable`**. The Watson U² path (`core_snpm_circ`) is the only caller that sets
it, at `dh = 0.005`; every t- and F-scale caller omits it. Pinned by `test_circ_snpm` T2.

`params.tail` ∈ `both` / `left` / `right` for the legacy t-tests and correlations. The GLM
and LMM presets pick statistic, contrast and permutation scheme automatically from the design.
The circular keys **lock** `tail` to `both` and `datatype` to `absolute`, and additionally
require `circ_units` and `circ_convention` (neither has a default) — see the
[how-to](../how-to/run-a-circular-phase-analysis.md).

`circ_wheeler_watson_Test` and `circ_WatsonsU2Test` are **retired**, together with the
`snpm_single_threshold_with_TFCE_circ` / `snpm_cluster_analysis_circ` engines. Use the three
`core_snpm_circ` keys above.
