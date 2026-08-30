# Run an analysis in the GUI

Drive `SnPMAnalysisGui` (title "hdEEG-SnPM Toolbox") to run one analysis and write its report.
This assumes you already know which analysis you need — if not, see
[choosing an analysis](../explanation/choosing-an-analysis.md).

## Launch

From the repo root (paths must resolve from here):

```matlab
SnPMAnalysisGui
```

## Pick the data source

Use the **Data source** toggle:

- **Files** — one or more CSV/Excel tables you already have. Use this for every analysis.
- **Spectral folder** — one EEG_processor BIDS power-spectrum folder per cohort/condition,
  each folder a level of the design factor. Only `unpairedT` / `pairedT` / `onesampleT` /
  `anova1` / `rmanova` are available here; the factor is set automatically, so the Roles panel
  is hidden.

## Files mode

1. Choose the analysis in **Comparison**. Labels are grouped under inert `--- ... ---` headers;
   pick a real analysis, not a header.
2. Load your table(s). For the single-file GLM presets (`anova1` / `ancova` / `regression` /
   `rmanova` / `mixed2way`) only **Data File** shows — Data 2 is hidden. For the legacy
   two-file tests, load **Data 1** and **Data 2** (condition A/B, group A/B, measure 1/2, or
   angles A/B, per the analysis). For the **correlations**, Data 2 may be either a second
   subjects × channels file or just a `Subject` column plus **one** behavioural measure column
   (e.g. a reaction time) — the single value is broadcast across Data 1's channels. A file with
   several unrecognised columns is matched to Data 1 **by position** if the widths agree
   (with a caveat banner in the report) and refused otherwise, so keep a behavioural file to one
   measure column. For an EEG-vs-behaviour run keep **Data Type = `absolute`**, as the GUI
   hint says: under `normalize` each subject's EEG row is z-scored across channels and sums to
   zero, so the single-score r-map is forced to balance positive against negative channels and
   a global association cannot show; the run warns, marks the results text `CAVEAT:` and puts an
   amber **Data Type caveat** banner in the report (see the
   [catalog](../reference/ANALYSIS_CATALOG.md)). For **any** two-file comparison, Data 1 must be
   `Subject` plus channel columns only — a metadata column there shifts the channel pairing and
   the run refuses with `snpm:correlation:channelOrderMismatch` naming the column.
3. On load, the GUI auto-detects columns with the same classifier the engine uses
   (`snpm_detect_behavior_table`): a column is a channel if its name matches `^E\d+$` or `Cz`
   (case-insensitive) **or** is a label of any recording system in `snpm_montage_registry`;
   everything else is metadata. Check the **Detected: N channels, M metadata columns** caption.
4. Set the **Recording system** to match your channel labels — `E1..E256` / `Cz` → EGI 256
   (HydroCel); `Fp1..` → Compumedics 257. The channel count must match the montage.
5. Map the **column roles** the analysis needs (only the required pickers are shown):
   - `anova1` → Group
   - `ancova` → Group + Covariate cols
   - `regression` → Predictor (+ optional Covariate cols)
   - `rmanova` → Subject + Condition
   - `mixed2way` → Group + Subject + Condition
   - legacy t-tests / correlation → no roles
6. Set **Data Type** (`absolute` / `logscale` / `normalize`), **Permutations** (use ~500–1000
   for a quick check, 10000+ for final inference), and the **output** path.
7. When every required field is set, **Run** enables. Click it.

## Spectral-folder mode

1. Pick the analysis (one of the five supported).
2. **Add folder…** once per level; edit each row's label in the table. Paired / 2-group tests
   need exactly 2 folders; `anova1` / `rmanova` need ≥2. Labels must be unique and non-empty.
3. The first added folder fills the **Band** / **Stage** listboxes and channel count. Select
   one or more bands and stages — the run loops every band×stage cell.
4. Pick **Power type** (`absolute` / `normalized`), permutations, and output path, then **Run**.

## Outputs

Every run writes to the output path, per analysis (per band×stage cell in spectral mode):

- `<base>_<timestamp>.mat` — `results_struct` (`T`, `p`, `Clusters`, `correctTFCEsigch`,
  `SnPMsigch`, `chanlocs`, …).
- `<base>_<timestamp>.xlsx` — significant-channel table.
- `<base>_<timestamp>_report.html` — one self-contained report (images embedded as
  `data:image/png;base64,` URIs, so the single file can be moved, emailed or imported into
  Notion; the PNGs are still written beside it, and `params.report_embed_images = false`
  switches back to relative filenames): summary strip, model card,
  descriptive maps, one section per effect with a Cluster / TFCE / Uncorrected toggle and
  numbered clusters (post-hoc pairwise tabs for the omnibus ANOVA presets), and a Global test
  card. Its **Results CSV** button opens the sibling `<base>_<timestamp>_report.csv`.
- Topoplot PNGs.

See the [Analysis catalog](../reference/ANALYSIS_CATALOG.md) for the exact fields per analysis.

## If Run stays disabled

`checkReadyToRun` is analysis-aware. Confirm: the right number of files/folders is loaded, the
channel count matches the montage, and **every visible role** is set. `mixedmodel` is not in
the GUI — run it with [`scripts/run_lmm_analysis.m`](https://github.com/TancyKao/hdEEG-SnPM/blob/main/scripts/README_scripts.md).
