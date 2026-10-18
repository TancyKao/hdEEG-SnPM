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
python3 db_to_group_table.py <DATA_dir> <OUT_dir> --group CONTROL --stages NREM3
python3 db_to_group_table.py <DATA_dir> <OUT_dir> --group CONTROL --bands sw_0.5-1.25
```

`<DATA_dir>` is the cohort folder holding one subfolder per subject (each containing
`wonambi/neural_events.db` and `wonambi/*_clean_rebuilt.xml`). Always pass `--group`
explicitly — it carries a study-specific default.

| Option | Default | What it does |
|---|---|---|
| `--group LABEL` | `group1` | cohort label written into `subjects.csv` |
| `--stages LIST` | `NREM2,NREM3` | the sleep stages analysed, comma-separated (`NREM3`, `NREM2,NREM3`; the joint token `NREM2NREM3` is accepted too) |
| `--bands LIST` | all three | write only some event bands, by label (`spindle_9-12`, `spindle_12-15`, `sw_0.5-1.25`) |
| `--allow-hypnogram-denominator` | off | accept databases older than TurtleWave 4.3 and divide their density by raw hypnogram minutes (see below) |

It writes into `<OUT_dir>`:

| File | Layout |
|---|---|
| `eventStat_<param>_<event_band>.csv` | rows = subjects, columns = channels. `<param>` is `density`, `amplitude` or `duration` |
| `subjects.csv` | `Subject`, `group`, `session` — the metadata table every GLM preset joins on |
| `channels.csv` | the union of channel labels across subjects |
| `recordings.csv` | one row per subject and event band: `Subject`, `band`, `density_denominator` (`analysed_time` or `hypnogram`), `denominator_minutes`, `stages` (the stages density was divided over, joined with `+`), and `interpolated_channels` (the channels set to missing because TurtleWave interpolated them, repeated on each of the subject's rows). The three denominator columns are blank for a band never detected in that subject |

The event bands are fixed in the script: spindles 9–12 Hz, spindles 12–15 Hz, and slow waves
0.5–1.25 Hz. Amplitude and duration are per-channel means. The run-selection and density
rules below follow the hdEEG-SnPM Python reader and TurtleWave's own density function. Where
they overlap the three give the same numbers: on a TurtleWave 4.3+ database, with one
detection run per band inside `--stages`, density, amplitude and duration are equal, including
which channels are 0 and which are missing. They differ in these places:

- The bridge has a hypnogram fallback for old databases, behind
  `--allow-hypnogram-denominator`. The Python reader and TurtleWave have no density at all for
  such a database.
- The bridge cannot be told which detection method to read. Where the Python reader would
  accept a named method, the bridge stops with an error instead; leave that band out with
  `--bands` or use the Python reader.
- The bridge processes all three bands for every subject, so one band it refuses stops the
  whole run unless `--bands` leaves it out. A band with no event and no processed channel at
  all in a subject is written as a blank row, with a warning naming the subject and the band,
  and the run still finishes; the Python reader stops there.

**One detection run per band.** A TurtleWave database can hold more than one run of the same
band — for example a per-stage run stored under `NREM2` and a later joint run stored under
`NREM2NREM3` on the same night. Those runs re-detect the same spindles, so adding them together
counts most events twice. The bridge reads only the runs whose stored stages lie inside
`--stages`: with the default `NREM2,NREM3` a joint `NREM2NREM3` run, or a pair of per-stage
`NREM2` and `NREM3` runs, are both read. If two runs cover the same stage, or two detection
methods share the band, it stops with an error naming the database and the runs instead of
pooling them. Remove the unwanted run, choose `--stages` so only one run is covered, or leave
that band out with `--bands`. A joint run cannot be split: asking for `--stages NREM2` on a
database that only holds `NREM2NREM3` events stops with an error naming the database, the
stored stage token and the requested stages, as it does whenever a band's events all lie
outside `--stages`. A blank row there would look like missing data rather than a mismatch.

If a band has no event in scope at all, the channels TurtleWave processed for it get density
0. When TurtleWave's processing record lists more than one detection method for that band,
there is no way to tell whose channels those are, and the bridge stops rather than guess.

**Density is events per minute of analysed time.** TurtleWave 4.3 and later store the
artefact-free time of each stage (its `analysed_time` table). That table is shared by every
detector run on the recording, so it can hold stages a given band's run never searched. The
bridge therefore divides by the time of the stages the selected run actually searched, out
of those in `--stages`. It finds the searched stages the way TurtleWave's own density code
does: from TurtleWave's processing record for that method and band, or failing that from its
record of the detection run. A slow-wave run detected in NREM3 only and read with the default
`NREM2,NREM3` is therefore divided by NREM3 time, and the bridge prints
`density for sw_0.5-1.25 over NREM3 (run searched only these)` for that recording. Each of
those stages must have a stored time; if one is missing, the bridge stops and asks for the
right `--stages`.

Only when a database records neither (an old one) does the bridge fall back, with a warning,
to the stages its events are stored under, as TurtleWave and the Python reader do: a slow-wave
run whose events are all stored under `NREM3` is divided by NREM3 time even when read with the
default `NREM2,NREM3`. A stage the run searched but found nothing in leaves no stored event,
so its time is then missing from the denominator and the density can be too high. Re-running
detection with a current TurtleWave records the searched stages and removes the fallback.

Within one event band every subject must be divided over the same stages. If they are not, for
example one subject's run searched NREM2 and NREM3 and another's NREM1 to NREM3, the bridge
stops before writing anything and lists which subjects fall under which stage set, because a
rate over NREM2 and NREM3 and a rate over the whole of NREM are different quantities.

Older databases have no stored analysed time, and by default the bridge stops on them with an
error naming the database. The better fix is to re-run detection with a current TurtleWave.
If that is not possible, `--allow-hypnogram-denominator` divides their density by the raw
hypnogram minutes of the same stages, chosen by the same rule, read from the first rater in
the staging XML
(a warning is printed when the file holds several raters). Hypnogram time still contains the
artefact and arousal epochs the detector skipped, so this density is lower than the true rate
by an amount that depends on each recording's artefact load. That bias looks like a real
between-subject difference, so the bridge prints a warning for every such recording.

The two kinds of denominator are never mixed in one output folder: even with the flag, a run
whose recordings would use different denominators stops and lists which recordings use which.
`recordings.csv` records the denominator of every subject and band (`density_denominator` is
`analysed_time` or `hypnogram`, `stages` the stages it covers). Tables from separate runs, for
example one per `--group`, must not be merged into one comparison unless every subject in both
has the same `density_denominator` and the same `stages` for that band; check
`recordings.csv` before joining them.

**Zeros and missing values.** A channel TurtleWave processed but on which it found no event
gets density 0, not a blank. A channel TurtleWave recorded as interpolated (TurtleWave 4.5 and
later) is set to missing for every measure in that subject, and is listed in
`recordings.csv`, so the analysis tests it on the subjects where it was recorded or reports it
as excluded.

The bridge opens every database read-only, so it never writes into a TurtleWave database or
its write-ahead log.

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
`eventStat_*.csv` to `subjects.csv` on `Subject`. For `correlationP`/`correlationS`, Data 2
can be just a `Subject` column plus **one** behavioural measure column (a reaction time, a
memory score) — the value is broadcast across the channels; no need to copy it into a fake
channels table.

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
