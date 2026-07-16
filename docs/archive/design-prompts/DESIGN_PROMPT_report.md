# Design-assistant prompt — restyle the single-analysis report (spectral look, no facets) + add uncorrected channels & p-values

Hand this to the design/coding assistant with **three files**:
- `dependencies/generateAnalysisReport.m` — the function to rewrite.
- `templates/sleep_eeg_report.html` — the **visual target** (the polished spectral dashboard style).
- a current sample, e.g. `temp_check_csv/..._report.html` — the **before** (what it looks like now).

**Scope note (for whoever runs this):** `generateAnalysisReport.m` is the **standard per-analysis
report**, produced by `core_snpm_analysis` for *every* run — CSV imports (t-tests / correlation /
GLM presets / LMM) **and** spectral-folder analyses (one report per band×stage cell). The *faceted*
spectral dashboard (`export_report.m` → `sleep_eeg_report_filled.html`, Power×Stage×Band) is a
**separate** report and out of scope — we only want this single-analysis report to *look like* that
polished style, adapted to one comparison.

---

> **Task:** Redesign the HTML produced by my MATLAB function `dependencies/generateAnalysisReport.m`
> so it matches the polished look of my spectral reports, but for a **single statistical comparison**
> (no facets).
>
> **Three files I'm giving you:**
> - `dependencies/generateAnalysisReport.m` — the function to rewrite (builds self-contained HTML via `fprintf`).
> - `templates/sleep_eeg_report.html` — the **visual target**: card-based CSS, calm palette, a segmented
>   significance toggle, a "Reading the statistics" explainer, significance chips. Copy this *look*.
> - the current sample report — so you can see what's there now.
>
> **Critical adaptation — drop all faceting.** The spectral report is a grid of
> Power(Absolute/Relative) × Sleep-stage × Band. **My report has none of that** — it's ONE comparison,
> ONE set of topographies. So adopt the CSS, layout, cards, typography, toggle, and explainer, but
> **remove every facet control**: no power selector, no stage selector, no band selector, no
> periodograms, no by-band heatmap. Just one clean page.
>
> **What the redesigned report must contain (single comparison):**
> 1. **Header** — analysis name, the statistic (**t / F / correlation r**, already derived in the code),
>    the effect label, and a one-line methods sentence (n, permutations, α = .05, TFCE + cluster correction).
> 2. **Topographies** — the existing mean topo PNGs (group/condition means) + the stat map, in spectral-style cards.
> 3. **A three-way segmented toggle: Uncorrected · TFCE · Cluster.** Each view shows that view's
>    significance map **and a table of its significant channels** with columns
>    **Channel label · statistic value · p-value**, sorted by ascending p. For the **Cluster** view,
>    group channels under their cluster and show the **cluster-level p**.
> 4. **Post-hoc pairwise table** preserved for ANOVA (F) presets.
> 5. A short **"Reading the statistics"** explainer, like the spectral report's.
>
> **All data is already in `results_struct` — no new computation:**
> - `results_struct.p.real` — uncorrected per-channel p (`uncorrsigch` arg = channels with p ≤ .05).
> - `results_struct.p.correctedTFCE` — TFCE-corrected per-channel p (`correctTFCEsigch` arg).
> - `results_struct.Clusters` — struct array, each `.channels` + cluster-level `.p`
>   (`SnPMsigch` arg = channels in clusters with p ≤ .05).
> - `results_struct.T.real_T` — the statistic per channel. `results_struct.chanlocs(i).labels` —
>   channel name (e.g. "E129"); **use labels, not indices**.
> - `results_struct.glm.contrast_type` ('t'/'F') + `.effect_label`; `results_struct.posthoc`
>   (`.label`, `.correctTFCEsigch`, `.SnPMsigch`).
> - α = 0.05; `params.comparison`, `params.permutations`, `params.tail`, `params.datatype`.
>
> Per view: Uncorrected p = `p.real(ch)`; TFCE p = `p.correctedTFCE(ch)`; Cluster p = the `.p` of the
> cluster that channel belongs to; statistic for any channel = `T.real_T(ch)`.
>
> **Constraints (keep exactly):**
> - Return the **full updated `generateAnalysisReport.m`** — same `fprintf`-built, **self-contained**
>   HTML (inline CSS + JS, references the sibling PNGs, no runtime template-file dependency).
> - Keep the signature
>   `generateAnalysisReport(results_struct, params, base_filename, outputSname, uncorrsigch, correctTFCEsigch, SnPMsigch)`
>   and output filename `<outputSname>_report.html`.
> - Work for **all analyses** (t-tests, correlation, ANOVA F + post-hoc, LMM) and handle **empty**
>   significant sets gracefully.
> - Permutation p can be 0 → display as `< 1/nperm` (from `params.permutations`), never "0.000".
> - Must pass `checkcode('dependencies/generateAnalysisReport.m')` with no errors and construct
>   headless. HTML-escape any free text.
>
> Return the full updated `.m` plus a one-paragraph summary of changes.

---

## Verify when it returns (do this before relying on it)
- `matlab -batch "checkcode('dependencies/generateAnalysisReport.m')"` — no errors.
- `matlab -batch "test_all"` — `VERIFY: ALL GREEN` (the `outputs-asserted` check confirms each analysis
  still writes `<base>_report.html`).
- Regenerate a sample (run a `unpairedT` + an `anova1` to a folder) and open the two `_report.html`
  files: confirm the **Uncorrected/TFCE/Cluster** toggle, the channel + p-value tables, and the
  F-map + post-hoc for the ANOVA. The page should have **no** stage/band/power selectors.
