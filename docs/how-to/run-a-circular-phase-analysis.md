# Run a circular (phase) analysis

Compare **preferred phase** — where on the slow oscillation a subject's spindles (or other
events) sit — between two independent groups, or correlate it with a per-subject measure.

> **Check you need this first.** Coupling *strength* and event *prevalence* are ordinary linear
> numbers and go through `unpairedT` / `anova1` / `regression` today, with more power and no
> circular machinery. Only *preferred phase* needs this page. See
> [About circular statistics for phase](../explanation/circular-statistics-for-phase.md).

## Before you start

- You have one wide CSV of **angles** per group (or one angle file plus one measure file).
- You know which **zero convention** your detector uses — see
  [Pick the convention](#pick-the-convention). Getting this wrong is the single most common
  failure on this path.
- You have a per-channel **event count** file per group — `count1_file` / `count2_file`, labelled
  **"Group A / Group B Event Count File"** in the GUI. This is a required input, not an optional
  one: it is the precision covariate that stops a difference in how many events each group had
  from masquerading as a difference in phase.

## Export the right columns upstream

If you are running [TurtleWave-hdEEG](https://github.com/TancyKao/TurtleWave-hdEEG)'s
phase–amplitude coupling module across a cohort — see
[Take TurtleWave detections into a group analysis](analyse-turtlewave-event-data.md) for the
rest of that hand-off — export **per subject per channel**:

| Column | Why you need it |
|---|---|
| **preferred phase** | the analysed quantity |
| **`n_events`** | **what the engine reads today.** This is the column that goes into `count1_file` / `count2_file` (the GUI's *Event Count File*), and non-integer values are rejected. |
| `rayleigh_z` | the precision covariate the analysis is moving to. Z = nR² captures both how many events there were and how tightly they clustered. Export it now so you are ready. |
| `rayleigh_p` | the descriptive quality panel, and an optional explicitly-declared sensitivity analysis |

**None of the three gates the analysis.** There is no Rayleigh mask — cells that fail the
Rayleigh test are kept, and their imprecision is handled by the covariate rather than by
deletion. The reasoning is in
[Why there is no Rayleigh mask](../explanation/circular-statistics-for-phase.md#why-there-is-no-rayleigh-mask).

Never put `rayleigh_z` and `n_events` in the same model. They are strongly collinear and the
second adds nothing.

> **Current build.** The engine reads `count1_file` / `count2_file` as **integer event counts**
> and rejects non-integer values (`core_snpm:circCountsInvalid`), which is why the GUI calls
> them the *Event Count* files. Export `rayleigh_z` now — it is the covariate the analysis is
> moving to — but feed `n_events` into those two files until the queued engine change lands.

## File layout

Every file is **wide**: one row per subject, one column per channel, plus a `Subject` column.
Channel columns use real montage labels (`E1..E256` / `Cz` for EGI).

**Two independent groups** (`circ_phase_group`, `circ_phase_group_u2`) — four files:

| File | GUI picker | Contents |
|---|---|---|
| `data1_file` | Group A Angles File | group A angles |
| `data2_file` | Group B Angles File | group B angles |
| `count1_file` | Group A Event Count File | group A per-channel event counts |
| `count2_file` | Group B Event Count File | group B per-channel event counts |

The two event-count files must have **exactly the same channel columns** (order free) and
**exactly the same subject ids** as their angle file. Mismatches are hard errors that name the
offending entries.

**Circular–linear** (`circ_corrAngLinear`) — two files, no event-count file: `data1_file` =
angles (GUI: *Angles File*), `data2_file` = a `Subject` column plus one numeric measure column
(GUI: *Measure File (linear)*; name the column in `measure_col` if the file has more than one).

### Subject ids: disjoint or matched

| Analysis | Subject ids |
|---|---|
| `circ_phase_group`, `circ_phase_group_u2` | **disjoint** across the two angle files — these are independent groups |
| `circ_corrAngLinear` | **matched** between the angle file and the measure file — paired per subject |

If your two angle files contain the same subjects, you have a within-subject design and this is
the wrong test. Use the signed linearised measure with `pairedT` or `rmanova` instead.

## Angles must be in radians (or declare degrees)

Set `circ_units` explicitly — in the GUI this is **Angle units** in the
**"Phase / angle settings (required)"** panel, which appears as soon as you pick a phase
analysis. There is **no default** — reading degrees as radians scrambles every angle and nothing
downstream can detect it, so the engine errors (`core_snpm:circUnitsRequired`) rather than
guessing.

- Already in radians on `[0, 2π)` or `(−π, π]` → `circ_units = 'rad'`.
- Detector exported degrees → either `circ_units = 'deg'` (the engine converts), or convert
  yourself with `deg2rad`.

## Pick the convention

Detectors disagree about which point of the slow oscillation their angle calls zero. Set
`circ_convention` to match your detector — in the GUI this is **Where zero sits**, next to
Angle units in the same "Phase / angle settings (required)" panel; picking *Custom offset*
enables the **Custom offset (deg)** field, which is `circ_zero_offset_deg`:

| `circ_convention` | Detector | Where zero sits | Rotation applied |
|---|---|---|---|
| `literature_uppeak0` | the published standard | up-state positive peak | 0° |
| `yasa_uppeak0` | YASA | up-state positive peak | 0° |
| `luna_zerocross0` | Luna `COUPL_ANGLE` | positive-to-negative zero crossing, so the up-state peak sits at **270°** | +90° |
| `turtlewave_pre_v4` | TurtleWave **before v4.0** | exactly **180° wrong** | +180° |
| `custom` | anything else | you declare it in `circ_zero_offset_deg` | as declared |

**TurtleWave exports predating v4.0 are exactly 180° wrong.** This was an upstream bug in
[TurtleWave-hdEEG](https://github.com/TancyKao/TurtleWave-hdEEG): the phase-bin centres spanned
`[0, 2π)` while the amplitude binning used `[−π, π)`. It was fixed in
[commit `d341690`](https://github.com/TancyKao/TurtleWave-hdEEG/commit/d3416909571b7f6fa39c624b72136acf8fc566ef)
on 16 July 2026. Recording **polarity inversion** adds a further 180° on top,
so an inverted recording exported from a fixed TurtleWave is back to correct, and a
non-inverted recording from an old TurtleWave is not.

This is not a hypothetical. It was verified on real data in this project: 514 rows with a grand
circular mean of **+158.7°**, which repairs to **−21.3°** once the 180° rotation is declared —
squarely inside the published range.

The engine will not let this through silently. It computes the pooled grand mean over frontal
channels and hard-errors with `core_snpm:circPhaseConventionInverted` if it lands in the inverted
window. When that fires you have two options, and only two:

1. Re-export from TurtleWave v4.0 or later, then run with `literature_uppeak0`.
2. Declare `circ_convention = 'turtlewave_pre_v4'`, which applies the +180° repair.

Do not "fix" it by rotating the CSV by hand — the applied rotation is recorded in the results
struct and printed in the report header, and a hand-rotated file loses that audit trail.

## Angular resolution

Some exports are whole-degree (Luna's `COUPL_ANGLE`, TurtleWave's degree column). What that costs
depends entirely on the statistic, so the gate is path-dependent:

- **`circ_phase_group_u2` — hard error** (`core_snpm:circResolutionTooCoarse`). Watson's U² is
  built on the empirical cumulative distribution, so ties are structural: whole-degree input
  moves the uncorrected p by up to **0.084** and drops the topographic rank correlation to 0.958.
- **`circ_phase_group` and the linearised measures — warning only**
  (`core_snpm:circResolutionCoarse`). Both are smooth functions of the angle with no ranks
  anywhere; a ≤0.5° perturbation against a 15–50° between-subject spread is a fraction of a
  percent.

Three remedies, in order of preference:

1. **Re-export at full precision.** Most detectors compute the angle in double precision and only
   round on write.
2. **Use `circ_phase_group` instead.** The Hotelling path accepts the data with a warning.
3. **Jitter, and declare it.** Add uniform noise on ±0.5° to break the ties. This is a last
   resort, it changes your data, and it must appear in your methods section.

## Run it

### From the GUI

Run `SnPMAnalysisGui` from the repo root, then:

1. **Analysis** — under the `--- Phase / angle data (advanced) ---` heading pick
   *Preferred phase: 2 groups (recommended)* (`circ_phase_group`),
   *Preferred phase: 2 groups (Watson U2, alternative)* (`circ_phase_group_u2`), or
   *Preferred phase vs a linear measure (paired by subject)* (`circ_corrAngLinear`).
2. **Files** — the pickers relabel themselves: **Group A / Group B Angles File** plus **Group A /
   Group B Event Count File** for the two group analyses, or **Angles File** + **Measure File
   (linear)** for circular–linear. Each picker has its own sheet dropdown for Excel input.
3. **Phase / angle settings (required)** — this panel appears with the analysis. Set **Angle
   units** and **Where zero sits**; both start on *(choose)* and the run will not start until you
   set them. *Custom offset (deg)* is enabled only for the Custom convention.
4. **Data Type and Tails are fixed** for these analyses (absolute, two-sided) and the controls are
   disabled — see the note at the bottom of the panel.
5. Set channels, permutations and the output folder as usual, then **Run**.

### Headless

From the repo root:

```matlab
params = struct( ...
    'comparison',       'circ_phase_group', ...     % or circ_phase_group_u2
    'data1_file',       'groupA_phase.csv', ...
    'data1_sheet',      'CSV File', ...
    'data2_file',       'groupB_phase.csv', ...
    'data2_sheet',      'CSV File', ...
    'count1_file',      'groupA_event_counts.csv', ... % REQUIRED ("Event Count File")
    'count2_file',      'groupB_event_counts.csv', ... % REQUIRED ("Event Count File")
    'channels',         'egi', ...
    'datatype',         'absolute', ...             % only 'absolute' is legal here
    'tail',             'both', ...                 % locked; see below
    'circ_units',       'rad', ...                  % REQUIRED, no default
    'circ_convention',  'literature_uppeak0', ...   % REQUIRED, no default
    'permutations',     10000, ...                  % 300 for a smoke test
    'output_path',      '/path/to/out');
results = core_snpm_analysis(params);
```

For circular–linear, swap the comparison and drop the count files:

```matlab
params.comparison  = 'circ_corrAngLinear';
params.data2_file  = 'behaviour.csv';
params.measure_col = 'measure';
params = rmfield(params, {'count1_file','count2_file'});
```

**`datatype` and `tail` are not free choices.** `logscale` and `normalize` are undefined on a
circular domain and are refused (`core_snpm:circDatatypeNotSupported`). All three statistics are
non-negative omnibus quantities with a single upper tail, so `left` / `right` are refused too
(`core_snpm:circTailNotSupported`) and the GUI Tails control is disabled.

Try it on the synthetic fixtures first — run `gen_synthetic_testdata` to write them, then use
`test_data/synthetic_gui/circ_phase_group{A,B}.csv` plus `circ_phase_counts{A,B}.csv`, radians,
`literature_uppeak0`. See that folder's `README.md`.

## Read the diagnostics before the p-values

The run prints three things that decide whether the p-values mean anything. Read them in this
order.

1. **The convention line** — `Circular import: units=… convention=… -> rotation ±X deg applied`,
   then `Convention check: pooled grand mean direction … = N deg`. The grand mean should land
   near the up-state peak. If it is near 180°, your convention is wrong.
2. **The event-count / group association** — `point-biserial r(log …, group): median …, max |r| …`.
   Above about **0.7** the covariate and the group contrast are nearly the same variable: the
   analysis will complete without error and will nonetheless be unable to separate phase from
   coupling strength. Report the two side by side instead.
3. **`Estimability: N/M channels analysed (K excluded)`** — and then the excluded-channel list.
   Look at **where** the exclusions are, not just how many. If the parieto-occipital channels are
   gone, you cannot claim a frontal effect is regionally specific.

## Read the outputs

Same output set as every other analysis, written to `output_path`:

- `<base>_<timestamp>.mat` — `results_struct`, including `results_struct.circ` (the applied
  rotation, the convention, the covariate report, the evaluable-channel mask) and
  `results_struct.circ_descriptive` (per-group mean direction, circular SD, Rayleigh — **labelled
  no-inference**).
- `<base>_<timestamp>.xlsx` — the significant-channel table, plus an **`excludedChannels`** sheet.
- `<base>_<timestamp>_report.html` — the combined report, with the applied rotation in the header
  and an excluded-channel banner.
- Topoplots.

**There is no global/omnibus test card**, deliberately: averaging angles arithmetically across
the head is meaningless, and a circularly correct whole-head mean would average over a real
25–50° anterior-posterior phase gradient.

**What a significant cluster means:** the two groups' phase distributions differ somewhere in
that region. That may be a mean-direction shift, a concentration difference, or both — the test
does not separate them, and it does **not** license a directional claim. For a directional
result, run the signed linearised measure through `unpairedT` / `ancova` instead. For
circular–linear, significance means the measure varies with phase, never that it is better at any
particular phase.

## Troubleshooting

**`core_snpm:circUnitsRequired`**
You did not set `circ_units`. There is no default on purpose.

**`core_snpm:circPhaseConventionInverted`**
Your angles are near-antiphase to the up-state peak. Re-export from TurtleWave v4.0+, or declare
`circ_convention = 'turtlewave_pre_v4'`. Do not rotate the CSV by hand.

**`core_snpm:circResolutionTooCoarse`**
Whole-degree angles on the U² path. Re-export at full precision, switch to `circ_phase_group`, or
jitter and declare it.

**`core_snpm:circCountsRequired`**
The two group analyses need `count1_file` and `count2_file` (the GUI's *Group A / Group B Event
Count File*). Without the event-count covariate there is no way to tell a genuine phase
difference from a difference in how well each group's angles were estimated.

**`core_snpm:circColumnMismatch` / subject-id mismatch**
The event-count file's channel columns or subject ids do not match its angle file exactly. The
error names the offending entries.

**`core_snpm:circNoEvaluableChannels`**
Every channel failed the ≥8-per-group gate, the completeness requirement, or the collinearity
check. Check the excluded-channel count and the point-biserial diagnostic.

**`core_snpm:circTailNotSupported` / `core_snpm:circDatatypeNotSupported`**
Set `tail = 'both'` and `datatype = 'absolute'`. Neither is a free choice here.

## Related

- [About circular statistics for phase](../explanation/circular-statistics-for-phase.md) — why
  the tests work this way, the power figures, and what you may and may not claim.
- [Analysis catalog](../reference/ANALYSIS_CATALOG.md) — the import contract per comparison key.
- [Interpreting cluster results](../explanation/interpreting-cluster-results.md).
