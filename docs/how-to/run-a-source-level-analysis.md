# Run a source-level (2447-voxel) analysis

Run an SnPM permutation t-test on **reconstructed cortical band power** (2447 sLORETA/GeoSource
voxels) instead of scalp electrodes. The statistic core is spatial-unit-agnostic: once the data
is in the source layout, it flows through the same permutation-t / TFCE / cluster-extent engines
as the scalp t-test. This path is **headless only** (not in the GUI) — you drive it with a
`params` struct.

> **Precondition — coordinates are PLACEHOLDER.** The shipped `source2447_coords.mat` is a
> deterministic graph-Laplacian layout, **not real MNI space**. Statistics, significance flags
> and cluster membership in the output are exact; the `X/Y/Z` columns of the voxel table are
> placeholder until you replace the asset with your real GeoSource MNI export (see
> [Swap in real coordinates](#swap-in-real-coordinates)). Do not report voxel coordinates until
> you have done this.

## Before you start — the export contract

Export from GeoSource **two** things, in the SAME 2447-node order:

1. **The data.** A subjects×2447 table of cortical band power. Requirements:
   - Columns named exactly **`src0001 .. src2447`**, in canonical graph-node order (column *i* =
     neighbour-graph row *i* = voxel *i*). `snpm_assert_source` errors (`snpm:source:order`) if
     the columns are shuffled or misnamed.
   - Values are **band-power MAGNITUDES** (non-negative), **not** signed current density. The
     `absolute` group t-test rejects negative input (`snpm:source:signedMagnitude`) because
     arbitrary per-subject sign cancels a group magnitude effect. One wide CSV per condition/
     group, rows aligned across the paired files.
2. **The coordinate table** (for later): a 2447-row `label,X,Y,Z` CSV in graph-node order, to
   drop in via `build_source2447_coords` (below).

**Uniform reconstruction is a precondition the toolbox cannot check.** Every subject must share
the same head model / inverse operator / voxel grid, or voxel *i* is not the same location
across subjects and the group test is meaningless. Confirm this in your GeoSource pipeline.

## Run it

From the repo root, build a `params` struct and call `core_snpm_analysis`. Set
`channels = 'source2447'` (aliases `source` / `sources` / `loreta2447` / `voxel2447` also work):

```matlab
params = struct( ...
    'comparison',    'pairedT', ...          % or onesampleT / unpairedT
    'data1_file',    'condA_source.csv', ...  % subjects x 2447, cols src0001..src2447
    'data1_sheet',   'CSV File', ...
    'data2_file',    'condB_source.csv', ...
    'data2_sheet',   'CSV File', ...
    'channels',      'source2447', ...        % selects the source system
    'datatype',      'absolute', ...          % magnitudes -> non-negative guard active
    'tail',          'both', ...
    'permutations',  10000, ...               % 100 only for a smoke test
    'output_path',   '/path/to/out', ...
    'use_covariates', false, 'covariate_file', '');

results_struct = core_snpm_analysis(params);
```

For a trial-level (many-observations-per-subject) source design, use the LMM path instead:
`comparison = 'mixedmodel'` with `channels = 'source2447'` (see
[`scripts/run_lmm_analysis.m`](https://github.com/TancyKao/hdEEG-SnPM/blob/main/scripts/README_scripts.md)). The GLM presets
(`anova1`/`ancova`/`regression`/`rmanova`/`mixed2way`) are **not yet source-wired** — run those
on scalp montages.

## Read the output

Source runs **bypass the scalp topoplot** (a 2-D disc is meaningless for cortical voxels — the
run writes **zero** scalp PNGs) and emit a significant-voxel list instead:

- `<base>_<timestamp>_sigvoxels.csv` — one row per voxel: `idx, label, X, Y, Z, stat, p, pTFCE,
  sig_uncorr, sig_TFCE, in_sig_cluster, cluster_id, cluster_p`. The same table is also a
  `sigVoxels` sheet in the `.xlsx`.
- `<base>_<timestamp>.mat` / `.xlsx` / `_report.html` as for every analysis;
  `results_struct.is_source` is `true`.

**TFCE is the primary correction.** The cluster statistic on the source **t-test** path is
cluster **extent** (size) — the same extent-based FWE correction the scalp t-tests use, since the
source t-test reuses `snpm_cluster_analysis` unchanged. Cluster **mass** applies only to the
source **LMM** path (`core_snpm_lmm`, Stephan's mean-Wald recipe). Filter the voxel table on
`sig_TFCE == 1` or `in_sig_cluster == 1`.

## Interpret with care

- **Report REGIONAL, not per-voxel, effects.** sLORETA has low spatial resolution and spatial
  leakage; a single significant voxel is not a reliable focal claim. Describe significant
  *regions* (contiguous clusters), not individual voxels.
- **Volumetric adjacency bridges sulci and hemispheres.** The 26-connectivity source graph joins
  voxels that are close in 3-D but may sit on opposite banks of a sulcus or across the midline, so
  a cluster can span anatomically distinct areas. Inspect cluster extent against anatomy.

## Swap in real coordinates

When you have the real GeoSource MNI export (`geosource_mni.csv`, columns `label,X,Y,Z` in
graph-node order):

```matlab
T = readtable('geosource_mni.csv');
load('dependencies/NeighborMatrix_Sources_2447_Full.mat');   % node order reference
source2447_coords = struct('labels', cell(1,2447), 'X',[], 'Y',[], 'Z',[]);
for i = 1:2447
    source2447_coords(i).labels = sprintf('src%04d', i);
    source2447_coords(i).X = T.X(i);  source2447_coords(i).Y = T.Y(i);  source2447_coords(i).Z = T.Z(i);
end
provenance = struct('placeholder', false, 'space', 'MNI', 'source', 'geosource_mni.csv');
save('dependencies/source2447_coords.mat', 'source2447_coords', 'provenance');
```

Keep node order == graph node order; `snpm_assert_source` enforces `src0001..src2447` in
sequence and checks the coordinate labels agree. After the swap the `X/Y/Z` columns in
`_sigvoxels.csv` are real MNI.

## Verify

`matlab -batch "test_source_snpm"` — graph sanity + assert guards, planted-cluster recovery
(TFCE + cluster-extent), pure-null and scheme-mismatch negative controls, per-voxel t vs MATLAB
`ttest`, and the magnitude guard. See the [analysis catalog](../reference/ANALYSIS_CATALOG.md)
for the full source-space input/output contract.
