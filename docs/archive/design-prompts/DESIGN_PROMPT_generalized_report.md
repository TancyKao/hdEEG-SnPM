# Design-assistant prompt — generalise the report to ALL analyses

Hand this to the design assistant **together with the current finalised report**
(`…/hdEEG_overnightPSA/sleep_eeg_report_filled.html`). It keeps that report's look/behaviour as the default
and adds a few optional, data-driven `REPORT` fields so ONE HTML renders every analysis we run. Supersedes
the earlier `DESIGN_PROMPT_other_analyses.md`.

---

> Here is my working report `sleep_eeg_report.html` (renders from one `REPORT` object + fixed-name images).
> **Keep the current visual design and behaviour exactly as the default** — do not regress any of these:
> - 3 faceting controls: a **segmented toggle** (currently Power: Absolute/Relative), **tabs** (currently
>   sleep Stage), and **rows** within each tab (currently frequency Band).
> - **Global Overview**: one periodogram image per tab — a single-panel PSD plot (gray shading marks
>   significant frequency bins; NO separate p-value panel).
> - **Regional Overview heatmap**: tabs × rows, cells coloured by significant-channel count, metric toggle
>   (Cluster / TFCE / Uncorrected).
> - **Per-tab**: a band-grid image; a **per-row detail** strip of topoplot images; a **summary table**
>   (significant rows highlighted, click a count to reveal channel labels). No per-row "badge" pill.
> - Topoplots are square; legend explains the diverging T-map (warm = A>B) and significance dots
>   (black = uncorrected, white = corrected).
>
> Now generalise it so the SAME file also renders other analyses, driven entirely by NEW OPTIONAL `REPORT`
> fields. When a field is absent, behave exactly as today. Do not break the current report.
>
> **1. Generic axes** (replace the hardcoded Power/Stage/Band semantics with named axes; default to today's):
> ```
> axes: {
>   toggle: { name:"Power",  levels:[{key:"absolute",name:"Absolute"},{key:"normalised",name:"Relative"}] } | null,
>   tabs:   { name:"Stage",  levels:[{key:"n1",name:"N1"}, ...] },
>   rows:   { name:"Band",   levels:[{key:"sigma",name:"Sigma",sub:"12–15 Hz"}, ...] }
> }
> ```
> Use `axes.*.name` for headings/labels (e.g. show "Parameter" or "Event" instead of "Power"/"Band" when
> set). If `toggle` is null, hide the segmented control. Drive image filenames from the axis keys (below).
>
> **2. `analysis` block** (optional; absent ⇒ today's two-condition t-map):
> ```
> analysis: {
>   stat: "t" | "F",            // t = signed, diverging colormap + signed legend; F = non-negative, sequential, positive-only
>   effect_label: "Condition A vs B" | "Group (3 levels)" | "slope of MoCA" | "Group × Stage",
>   legend: "warm = A > B" | "F statistic",
>   model_info: "density ~ group + (covariates)"   // optional caption under the header
> }
> ```
>
> **3. `mean_panels`** — the descriptive topo columns in each per-row detail (optional; default = condA, condB):
> ```
> mean_panels: [ {label:"control mean", token:"g_control"}, {label:"case mean", token:"g_case"} ]   // 0..k columns
> ```
> Render however many are listed (0 → none, e.g. regression; 2 → today; k → one per group/condition),
> followed by the main **stat** map. Keep today's `condA`/`condB`/`Tmap` tokens working as the default.
>
> **4. `posthoc`** — optional pairwise/main-effect maps shown as an expandable block under a row (for F /
> ≥3-level designs), each with its own stat map image + significance + channel reveal.
>
> **5. `has_periodogram`** (bool, default true) — when false, hide the Global Overview periodogram section
> (for non-spectral analyses, e.g. event-based or LMM).
>
> **6. Keep** the `values[toggle][tab][row] = [n_uncorr,n_TFCE,n_cluster,min_p]` and
> `channels[toggle][tab][row] = {u,tfce,cluster}` schema. New fields are additive.
>
> **Image filename template** (generic; today's pattern is the default special case):
> `<tab>_<toggle>_<row>_<token>.png` for mean panels, `<tab>_<toggle>_<row>_stat.png` for the main map,
> `global_periodogram_<tab>_<toggle>.png` for periodograms. (Today: tab=stage, toggle=power, row=band,
> tokens condA/condB, stat=Tmap.)
>
> Return the single updated HTML file.

---

## How each analysis fills the generic schema (the MATLAB exporter sets these)

| Analysis | toggle | tabs | rows | stat | mean_panels | periodogram |
|---|---|---|---|---|---|---|
| **Two-condition spectral** (current) | Power (abs/rel) | Stage | Band | t | condA, condB | yes |
| Spectral >2-group ANOVA | Power | Stage | Band | F | k group means | yes |
| Spectral regression/ANCOVA | Power | Stage | Band | t | (none) | yes |
| Spectral group×condition | Power | Stage | Band | F | cell means | yes |
| **Event group (spindle/SO)** | Parameter (density/amp/dur) | (single, or none) | Event-band (slow-sp/fast-sp/SO) | t (2-grp) / F (>2) | per-group means | **no** |
| Event regression (cognition) | Parameter | (single) | Event-band | t | (none) | no |
| LMM | (effect levels) | Stage or effect | Band | t or F | optional | no |

Once the template supports these, the same `export_report` (generalised) emits a matching `REPORT` + images
per analysis — no further design work per analysis.
