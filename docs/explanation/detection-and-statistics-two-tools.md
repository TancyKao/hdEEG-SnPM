# Detection and statistics are two tools

> *Understanding-oriented.* Why sleep-event detection and group statistics live in two separate
> repositories, and what that means for you. To actually move data between them, see
> [Take TurtleWave detections into a group analysis](../how-to/analyse-turtlewave-event-data.md).

hdEEG-SnPM is the downstream half of a two-program pipeline. The upstream half is
[**TurtleWave-hdEEG**](https://github.com/TancyKao/TurtleWave-hdEEG) (Python ≥ 3.10, MIT), a
high-density-EEG sleep event detector built on Wonambi: it detects spindles, slow waves and
K-complexes, computes phase–amplitude coupling, batches in parallel for HPC, and ships two PyQt5
interfaces — `turtlewave_gui` for detection and `eeg_review_gui` for review. It answers *where
are the events in this recording*. hdEEG-SnPM answers *do these subject groups differ across the
scalp, corrected for testing every channel*.

## Why not one program

The two halves have almost nothing in common except the data that passes between them.

Detection is per-subject signal processing. It is I/O-bound and embarrassingly parallel, it wants
a review interface a human can scrub through recording by recording, and Python's scientific and
Wonambi ecosystems are where that work lives.

Group inference is per-channel permutation statistics over a *cohort*. It needs the montage,
the channel-adjacency structure, EEGLAB's `topoplot` for the maps, and MATLAB's Statistics
toolbox for `fitlm` / `fitlme` / `anovan` — which are also the ground truth every statistic in
this toolbox was checked against during development. It has no reason to see a raw recording at
all.

Splitting them at that seam means each half is usable on its own, and neither drags the other's
runtime into your environment. The interface between them is a file — a plain CSV of one number
per subject per channel — not a shared library, an API or a database schema. A file boundary is
easy to inspect, easy to version, and easy to produce from something that is not TurtleWave.

## TurtleWave is one source, not a dependency

This is worth stating plainly because the phrase "pipeline" invites the opposite reading.
**hdEEG-SnPM does not require TurtleWave, and does not require Python.** It takes a subjects ×
channels table. Common sources are:

- a spectral-power folder from a preprocessing tool (the built-in spectral-folder data source),
- TurtleWave event statistics via the `db_to_group_table.py` bridge,
- reconstructed cortical band power for the source-level path,
- any table you assemble yourself.

Nothing in the statistics knows or cares which of these produced the numbers. The only place the
upstream tool is named in the engine is the circular tier, and only because of a specific
historical defect, described next.

## The one place the coupling is not clean

Detectors disagree about which point of the slow oscillation their phase angle calls zero, and
one of them was, for a while, simply wrong: **TurtleWave before v4.0 exported preferred phase
exactly 180° inverted**, fixed upstream in
[commit `d341690`](https://github.com/TancyKao/TurtleWave-hdEEG/commit/d3416909571b7f6fa39c624b72136acf8fc566ef)
(16 July 2026).

That is why `circ_convention` exists as a required, defaultless input with a
`'turtlewave_pre_v4'` value, why the GUI offers "TurtleWave before v4.0 (180 degrees off)" in
its *Where zero sits* dropdown, and why the engine runs a pooled inversion check before it will
compute anything. An inverted angle is not malformed — it is a clean, tight, significant result
that means the opposite of what it says. A file boundary between two tools cannot protect you
from a convention mismatch; only declaring the convention can. The reasoning behind the check,
and its limits, is in
[Take TurtleWave detections into a group analysis](../how-to/analyse-turtlewave-event-data.md#the-phase-convention-trap--read-this-before-any-coupling-analysis)
and [About circular statistics for phase](circular-statistics-for-phase.md).
