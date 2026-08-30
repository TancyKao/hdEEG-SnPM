# Handoff prompt — SnPM report v2: apply the master template to all analyses

Paste the prompt below into Claude Code, attaching `snpm_report_v2_master.html` (the visual spec) and pointing it at the report generator (`dependencies/generateAnalysisReport.m` in SnPM_2025).

---

## Prompt

`snpm_report_v2_master.html` is the approved v2 design for ALL hdEEG-SnPM report HTML outputs. It is built on a two-way mixed ANOVA (2×3) because that exercises every component. Refactor `generateAnalysisReport.m` so every analysis emits this template. Do NOT redesign: keep the CSS, class names, markup and inline JS exactly as in the mockup; only the content is generated.

### Universal structure (every analysis, in this order)

1. **Masthead** — eyebrow "SnPM Topographic Statistics", title, methods sentence (n permutations, tail, α, corrections), buttons: `Results CSV` (link to the sibling `_report.csv`) + `Print / PDF`.
2. **Sticky summary strip** (`.strip`) — one `.chipfx` per effect: traffic light (`.lamp.sig`/`.lamp.ns` = any corrected-significant channel), counts `Clu n · TFCE n · Unc n`, global stat + p, `href="#fx-…"` jump. Single-effect analyses get one chip + the Global chip.
3. **Model card** — formula line (spans `.v` = DV, `.r` = effect of interest) + spec rows: design, factors/levels, subjects (per group), data type, montage, permutations · α, correction. Adapt rows per analysis (e.g. covariates for ANCOVA, predictor for regression, random effects for LMM).
4. **Legend** — always all applicable scales: signed-t diverging bar and/or F hot bar, channel-dot key, ① cluster-number key.
5. **Descriptive block** — analysis-specific: cell-mean topo grid + interaction plot (mixed ANOVA), per-group means (ANOVA/ANCOVA), condition means (paired/rm), scatter (correlation/regression/LMM), circular summary (phase tests).
6. **One `.effect` section per effect** — collapsible; `open` class when the effect has any corrected-significant channel, collapsed otherwise (the mockup's Interaction shows the collapsed + empty states). Inside: effect map card, then Significant channels with the 3-way `.seg` toggle (`showView(fx,view)` — namespace `fx` per effect: 'g','c','i',…).
7. **Global test** — one `.gstat` per effect in a single `.global-card`.
8. **Footer** — self-contained note + generation timestamp.

### Cluster table rules (the key upgrade)

- Cluster panel table groups rows under `tr.grp-row` headers carrying `<span class="cnum">N</span>`, cluster p, channel count.
- For signed-t effects, insert `tr.sign-row` dividers: "Positive clusters (t > 0)" then "Negative clusters (t < 0)"; number clusters continuously across both.
- The same numbers must be drawn on the cluster-map PNG at each cluster's centroid (white disc, 1.5px accent #1f6f8b border, mono digit) — update the MATLAB topo-plot code accordingly.
- Empty result → `.empty-note` card stating the smallest p (see Interaction panels in the mockup).

### Table capping

Uncorrected tables >20 rows: first 20 visible, rest as `tr.xrow`, plus the `.morebtn` ("Show all N channels"). Print always keeps the cap (CSS handles it). Cluster/TFCE tables are never capped.

### Post-hoc (ANOVA / ANCOVA-3 / rmANOVA / mixed ANOVA)

Tabs per pairwise comparison (`.tabs` + `showTab`), count badge per tab (grey `.ct.none` when 0), each tab = TFCE map + table (Bonferroni-noted), empty tab = `.empty-note`. Print expands all tabs stacked with their `data-tab-name` labels.

### Per-analysis mapping

- **unpairedT / pairedT / onesampleT** — one t effect; descriptive = condition/group mean topos.
- **correlation** — one r effect (diverging bar labelled r); descriptive = scatter.
- **regression** — one t effect (slope); model card lists predictor + covariates; ADD the missing global test + descriptive scatter.
- **anova1 / ancova3 / rmanova** — one F effect + post-hoc tabs; descriptive = per-level mean topos; ancova model card lists covariates.
- **ancova (2 groups)** — t effect; ADD missing global test + descriptive.
- **mixed2way / 2x3** — exactly the mockup (group t, condition F, interaction F, post-hoc tabs under condition; simple effects may follow the same `.effect` pattern).
- **LMM** — keep its model card content (formula, random effects, trials·subjects); adopt strip/collapsible/cluster-numbering/capping from v2.
- **circular tests** — one effect (T²/U²/F, hot bar); descriptive = circular summary; no global test (note why in place of the card).

### Constraints

- Fully self-contained HTML: inline CSS + the mockup's small inline JS, no CDNs.
- p-format: mono, tabular-nums, no leading zero, `<.001` below.
- All images relative filenames; every image slot keeps the onerror grey "No data" placeholder.
- Keep channel rows, cluster headers, sign dividers, effect sections, strip chips and tabs as clean loop-generated blocks.

---

## Notes for you (not in the prompt)

- The mockup's CSV button targets a placeholder filename; the generator writes the real `<base>_report.csv` name.
- Centroid numbering requires a change in the MATLAB plotting code, not just the HTML.
