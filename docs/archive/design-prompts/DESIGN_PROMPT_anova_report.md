# Design-assistant prompt — ANOVA/omnibus-F report: per-group topos + F-map, and per-pair post-hoc significant channels

Edit **two files** — `core_snpm_glm.m` (must *produce* the new plots + carry the per-pair data) and
`dependencies/generateAnalysisReport.m` (renders them). The data the report needs is not generated yet,
so the report rewrite alone is not enough.

## Files to give the designer
**Edit:**
- `core_snpm_glm.m` — compute/save the per-group mean topos + F-map PNGs, store group labels and the
  per-pair stat/p/cluster data.
- `dependencies/generateAnalysisReport.m` — render per-group topos + F-map, and per-pair significant-channel tables.

**Reference (read-only, don't edit):**
- `dependencies/plot_topoInd.m` — has `create_topo_plot(data, chanlocs, title, filepath, ...)`, the single-map helper to reuse.
- `dependencies/snpm_glm_design.m` — where group levels/labels come from (`ref_dummy` → level list; `pairwise_contrasts`).
- `dependencies/snpm_perm_correction.m` — what the per-contrast `Tph` / `pph` / `Cl_ph` structs contain.
- `templates/sleep_eeg_report.html` — visual reference for the figure cards / styling.
- `temp_check_csv/glm_anova1_*_report.html` — the current ANOVA report (the "before"; note the broken
  `_Data1_topo.png` / `_Data2_topo.png` image links this task fixes).
- `temp_check_csv/unpaired_*_report.html` — the t-test report that must **not** regress.

---

> **Task:** Improve the ANOVA / omnibus-F report in my MATLAB SnPM EEG tool. Two changes: (1) the
> **Topographies** section must show a mean topography **for each group** plus the **F-map**; (2) the
> **Post-hoc** section must show, for each pairwise contrast, the **actual significant channels** broken
> down into **Cluster / TFCE / Uncorrected** (not just counts). Edit two files: `core_snpm_glm.m`
> (compute and save the new plots + carry the data) and `dependencies/generateAnalysisReport.m` (render them).
>
> **Important context — the data isn't generated yet.** For the GLM/ANOVA path, `core_snpm_glm.m`
> (line ~123) currently calls only `TopoplotSignificant_single(T.real_T, uncorrsigch, correctTFCEsigch,
> SnPMsigch, ...)`, which saves just the two significance maps (`<base> Cluster.png`, `<base> TFCE.png`).
> It does **not** make per-group mean topos or a standalone F-map. As a result the report's Topographies
> section references `<base>_Data1_topo.png` / `<base>_Data2_topo.png` / `<base>_topo.png`, which **don't
> exist for ANOVA** (broken image links). This task fixes that.
>
> **Change 1 — Topographies (each group + F-map):**
> - In `core_snpm_glm.m`, after the stats are computed, build the per-channel **group means** from the
>   data matrix `power` (rows = observations/subjects, cols = channels) grouped by the design's group
>   factor. The group **labels** come from `snpm_glm_design` (`ref_dummy` returns the level list `gl`;
>   expose/store it). Save one mean-topography PNG per group, e.g. `<outputSname>_group-<label>_mean.png`,
>   using the existing single-map helper `create_topo_plot(data, chanlocs, title, filepath, ...)` in
>   `dependencies/plot_topoInd.m` (this is what the legacy path's `plot_topoInd` wraps).
> - Save the **F-map** as its own PNG, e.g. `<outputSname>_Fmap.png`, from `T.real_T` (the per-channel F
>   for the F path) via the same `create_topo_plot`.
> - Store in `results_struct` what the report needs to enumerate these: the **group labels** (cellstr) and
>   the filenames (or a stable naming convention). Suggest `results_struct.glm.group_labels` and
>   `results_struct.glm.group_mean_png` (cellstr, parallel) + `results_struct.glm.fmap_png`.
> - In `generateAnalysisReport.m`, rewrite the Topographies section so that **for F (omnibus) analyses**
>   it lays out one card per group mean (titled by label) followed by the F-map card. For non-F analyses
>   (t / r), keep the current two-means + effect-map layout unchanged. Use the spectral-style figure cards
>   already in the file.
>
> **Change 2 — Post-hoc per-pair significant channels (Cluster / TFCE / Uncorrected):**
> - The pairwise loop in `core_snpm_glm.m` (~lines 98–110) already computes, per contrast, the full
>   per-channel result but discards most of it. Store it: for each post-hoc entry keep the per-channel
>   statistic `Tph.real_T`, uncorrected p `pph.real`, TFCE-corrected p `pph.correctedTFCE`, and the
>   cluster struct `Cl_ph` (each with `.channels` and `.p`), in addition to the index sets it already
>   stores (`uncorrsigch`, `correctTFCEsigch`, `SnPMsigch`). These pairwise contrasts are 2-group **t**
>   contrasts, so the statistic is **t**.
> - In `generateAnalysisReport.m`, replace the current post-hoc **count scoreboard** (function
>   `write_posthoc`, ~line 428) with, **per contrast**, a three-way **Cluster · TFCE · Uncorrected**
>   significant-channel presentation that mirrors the omnibus `write_significance` section: each view
>   lists its significant channels with columns **Channel label · t value · p-value**, sorted by ascending
>   p; the Cluster view groups channels under their cluster and shows the cluster-level p. Reuse the
>   existing `write_flat_table` / `write_cluster_table` helpers. Keep the contrasts in the order they
>   arrive. Use channel **labels** (`results_struct.chanlocs(i).labels`), not indices.
>
> **Data already available (no new stats math needed for the report):**
> - Omnibus: `results_struct.T.real_T` (per-channel F), `results_struct.p.real` / `.correctedTFCE`,
>   `results_struct.Clusters` (`.channels`,`.p`), `results_struct.chanlocs(i).labels`,
>   `results_struct.glm.contrast_type` ('F'), `.effect_label`.
> - Post-hoc: `results_struct.posthoc` struct array — currently `.label`, `.uncorrsigch`,
>   `.correctTFCEsigch`, `.SnPMsigch`; add the per-channel `t` / `p` / clusters per the above.
> - α = 0.05; permutation p can be 0 → display as `< 1/nperm` (use the existing `fmt_p`).
>
> **Constraints (keep):**
> - Keep the signature `generateAnalysisReport(results_struct, params, base_filename, outputSname,
>   uncorrsigch, correctTFCEsigch, SnPMsigch)` and output filename `<outputSname>_report.html`;
>   self-contained HTML (inline CSS/JS, references sibling PNGs).
> - Keep the significance toggle order **Cluster · TFCE · Uncorrected** (already set).
> - Don't regress the t-test / correlation / LMM reports — only the F branch changes layout; everything
>   else stays.
> - Both edited files must pass `checkcode` with no errors, and `matlab -batch "test_all"` must end in
>   `VERIFY: ALL GREEN`.
> - Handle empty significant sets and groups with missing data gracefully; HTML-escape free text.
>
> Return the two full updated files plus a one-paragraph summary.

---

## Verify when it returns
- `matlab -batch "checkcode('core_snpm_glm.m'); checkcode('dependencies/generateAnalysisReport.m')"` — no errors.
- `matlab -batch "test_all"` — `VERIFY: ALL GREEN`.
- Regenerate an `anova1` sample (synthetic: `test_data/synthetic_gui/glm_anova1.csv`, `group_col='group'`,
  `meta_cols={'Subject','group'}`, `channels='egi'`, low perms) and open the `_report.html`: confirm one
  mean topo **per group** + the **F-map** in Topographies, and per-pair **Cluster/TFCE/Uncorrected**
  channel tables (Channel · t · p) in Post-hoc, with no broken image links.
