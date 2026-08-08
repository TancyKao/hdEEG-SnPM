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
   angles A/B, per the analysis).
3. On load, the GUI auto-detects columns: `^E\d+$` / `Cz` are treated as channels, everything
   else as metadata. Check the **Detected: N channels, M metadata columns** caption.
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
- `<base>_<timestamp>_report.html` — one self-contained report with a TFCE / Cluster toggle,
  stat-aware (t-map for t-tests; F-map + post-hoc for the omnibus ANOVA presets).
- Topoplot PNGs.

See the [Analysis catalog](../reference/ANALYSIS_CATALOG.md) for the exact fields per analysis.

## If Run stays disabled

`checkReadyToRun` is analysis-aware. Confirm: the right number of files/folders is loaded, the
channel count matches the montage, and **every visible role** is set. `mixedmodel` is not in
the GUI — run it with [`scripts/run_lmm_analysis.m`](https://github.com/TancyKao/hdEEG-SnPM/blob/main/scripts/README_scripts.md).
