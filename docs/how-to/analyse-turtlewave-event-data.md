# Take TurtleWave detections into a group analysis

You have run [TurtleWave-hdEEG](https://github.com/TancyKao/TurtleWave-hdEEG) over a cohort and
want a corrected topographic map of where the groups differ. This page covers the hand-off:
what TurtleWave writes, how to turn it into the table this toolbox reads, which analysis each
kind of event measure belongs in, and the one trap that silently produces a wrong answer.

> **TurtleWave is a supported source, not a requirement.** hdEEG-SnPM consumes a subjects ×
> channels table of numbers. Where those numbers come from is your business — a spectral folder,
> your own detector, a spreadsheet you typed by hand. You do not need Python or TurtleWave
> installed to use this toolbox. If you *are* using TurtleWave, this page saves you the
> plumbing. For why the two are separate programs at all, see
> [Detection and statistics are two tools](../explanation/detection-and-statistics-two-tools.md).

## What TurtleWave produces

TurtleWave-hdEEG is a Python (≥ 3.10, MIT) high-density-EEG sleep event detector built on
Wonambi. It detects spindles, slow waves and K-complexes, computes phase–amplitude coupling,
runs in parallel for batch/HPC use, and ships two PyQt5 interfaces — `turtlewave_gui` for
detection and `eeg_review_gui` for review. Per subject it leaves, among other outputs:

| Output | What it holds |
|---|---|
| `wonambi/neural_events.db` | SQLite table of every detected event: channel, event type, frequency band, duration, peak-to-peak amplitude |
| `wonambi/*_clean_rebuilt.xml` | the staged hypnogram — epoch start/end and sleep stage |
| JSON detection results, XML annotations | per-channel event tables and annotations for review |

## What hdEEG-SnPM consumes

One number per subject per channel, wide: a `Subject` column plus one column per channel, with
channel columns named with real montage labels (`E1`…`E256`, `Cz`). Channel selection is
label-based and count-agnostic, so it does not matter that TurtleWave can export all 256
channels while your analysis montage is the 178-channel scalp set — the labels are matched
against the recording system and the data, chanlocs and neighbour matrix are subset into one
order. Hand it whatever survived your cleaning.

## Step 1 — bridge the database to group tables

`db_to_group_table.py` at the repository root does this. It is Python standard library only —
no MATLAB, no TurtleWave, no third-party packages — and it reads the SQLite database and the
staged XML directly rather than TurtleWave's parameter CSVs, so it survives changes to those.

```bash
python3 db_to_group_table.py <DATA_dir> <OUT_dir> --group CONTROL
```

`<DATA_dir>` is the cohort folder holding one subfolder per subject (each containing
`wonambi/neural_events.db` and `wonambi/*_clean_rebuilt.xml`). Always pass `--group`
explicitly — it carries a study-specific default.

It writes into `<OUT_dir>`:

| File | Layout |
|---|---|
| `eventStat_<param>_<event_band>.csv` | rows = subjects, columns = channels. `<param>` is `density`, `amplitude` or `duration` |
| `subjects.csv` | `Subject`, `group`, `session` — the metadata table every GLM preset joins on |
| `channels.csv` | the union of channel labels across subjects |

The event bands are fixed in the script: spindles 9–12 Hz, spindles 12–15 Hz, and slow waves
0.5–1.25 Hz. Density is **count divided by NREM2+NREM3 minutes**, taken from the hypnogram, so
it is already normalised for how much scorable sleep each subject had; amplitude and duration
are per-channel means.

## Step 2 — pick the analysis for the measure you have

Event data arrives at two granularities and they go to different engines.

**Aggregated per channel** — one number per subject per channel, which is what
`eventStat_*.csv` gives you (spindle density, amplitude, duration, slow-wave count). This is an
ordinary linear measure. It feeds:

| You want | Analysis |
|---|---|
| two independent groups | `unpairedT` |
| same people, two conditions or nights | `pairedT` |
| three or more groups | `anova1` (omnibus F + automatic pairwise post-hoc) |
| groups after removing age / sex / apnoea–hypopnoea index | `ancova` |
| relate density to a per-subject measure | `regression`, `correlationP`, `correlationS` |
| three or more conditions in the same people | `rmanova` |
| group × condition | `mixed2way` |

Run these from the GUI, or headless with `run_glm_analysis.m`, which joins one
`eventStat_*.csv` to `subjects.csv` on `Subject`.

**Event level** — one row per detected event, with subject, group and covariates repeated down
the rows. That is repeated-measures data and belongs in the per-channel linear mixed model,
`mixedmodel` (`core_snpm_lmm`), which is script-driven via `run_lmm_analysis.m`. Do not average
it down and then treat the averages as independent if the question is about within-subject
structure.

**Preferred phase** — the angle at which events couple to the slow oscillation is *not* a
linear number and must not go through the t-tests. It goes to the circular tier
(`circ_phase_group`, `circ_phase_group_u2`, `circ_corrAngLinear`); see
[Run a circular (phase) analysis](run-a-circular-phase-analysis.md). Coupling *strength* and
event *prevalence*, on the other hand, are ordinary numbers — use the linear table above and
save yourself the circular machinery.

## The phase-convention trap — read this before any coupling analysis

**TurtleWave exports predating v4.0 report preferred phase exactly 180° wrong.** The phase-bin
centres spanned `[0, 2π)` while the amplitude binning used `[−π, π)`. It was fixed upstream in
[commit `d341690`](https://github.com/TancyKao/TurtleWave-hdEEG/commit/d3416909571b7f6fa39c624b72136acf8fc566ef)
on 16 July 2026.

An inverted angle is still a perfectly well-formed angle. Nothing about it looks broken: you get
a tight distribution, a significant cluster, and a phase you can write a sentence about — and
the sentence says events couple to the down-state when they couple to the up-state. This is the
single most expensive mistake available on this path.

Two ways out, and only two:

1. Re-export from **TurtleWave v4.0 or later**, then run with `circ_convention =
   'literature_uppeak0'`.
2. If you are certain the export predates v4.0, declare it: `circ_convention =
   'turtlewave_pre_v4'`, which applies the +180° repair and records the rotation in the results
   struct and the report header. In the GUI this is **Where zero sits → "TurtleWave before v4.0
   (180 degrees off)"**.

Do not rotate the CSV by hand; that loses the audit trail.

The engine does try to catch you. `snpm_circ_check_inversion` pools both groups (so it cannot
contaminate the inference), takes the grand mean preferred phase over frontal estimable
channels, and hard-errors with `core_snpm:circPhaseConventionInverted` if it lands inside
[110°, 250°]. Real calibration from this project: an export sitting at +158.7° errors there and
repairs to −21.3° once declared. But it is a heuristic on one pooled number — a deviation
between 70° and 110° is only a warning, and less than 70° passes silently. **Declare the
convention; do not rely on the check to notice.**

Recording **polarity inversion** adds a further 180° on top, so an inverted recording exported
from a fixed TurtleWave is back to correct, and a non-inverted recording from an old TurtleWave
is not.

## Also export these columns

If you are running the coupling module across a cohort, export per subject per channel:
**preferred phase** (the analysed quantity), **`n_events`** (the required event-count file — it
is the precision covariate that stops a group difference in event count masquerading as a
difference in phase), and `rayleigh_z` / `rayleigh_p`. The full rationale, the required file
layout and the angular-resolution rule are in
[Run a circular (phase) analysis](run-a-circular-phase-analysis.md).

## What is not automated

There is no MATLAB importer for TurtleWave's CSV or SQLite event tables. `db_to_group_table.py`
covers the Wonambi database layout described above; any other export shape you assemble
yourself into the wide subjects × channels form. The columns are simple and the format is a
plain CSV.
