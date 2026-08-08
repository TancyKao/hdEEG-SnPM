# Getting started: your first SnPM analysis

In this tutorial we'll run a complete hd-EEG group analysis from start to finish — using the
toolbox's own synthetic data, which has a known effect planted in it. By the end you'll have
run a permutation test, produced a topographic significance map, and seen the tool correctly
recover a cluster it was designed to find. We'll use the graphical interface and pick the
simplest analysis so nothing distracts from the flow.

You don't need any statistics background or your own data to follow along.

> This tutorial is written to be reliable, but it hasn't yet been re-verified click-by-click on
> a fresh desktop MATLAB. If a label differs slightly from what you see, the surrounding step
> still applies. Please tell us so we can correct it.

## What you'll need

- MATLAB with the Statistics and Machine Learning Toolbox.
- This repository, opened with MATLAB's current folder set to the repo root. (The tool adds its
  own paths from the current folder, so this matters.)
- The synthetic data in `test_data/synthetic_gui/`. **It is not in the repository** — generate
  it first by running `gen_synthetic_testdata` from the repo root, which writes that folder
  (and a `README.md` inside it mapping every file to its analysis and GUI settings).

## Step 1 — Launch the tool

At the MATLAB command window, from the repo root, type:

```matlab
SnPMAnalysisGui
```

A window titled **hdEEG-SnPM Toolbox** opens. Notice the **Data source** toggle near the top —
leave it on **Files**, which is what we want for CSV data.

## Step 2 — Choose the analysis

Find the **Comparison** dropdown and select **Unpaired t-test** (its internal key is
`unpairedT`). We're choosing this because it compares two independent groups and needs no extra
column mapping — the simplest possible first run.

Notice that the file pickers now expect **two** files (Group A and Group B). Some analyses hide
the second file; this one keeps it.

## Step 3 — Load the two data files

Load these two files from `test_data/synthetic_gui/`:

- **Group A** → `unpaired_groupA.csv`
- **Group B** → `unpaired_groupB.csv`

As each file loads, the tool reads its headers. Watch for a caption like **Detected: 256
channels, 1 metadata column**. Each file's columns are `E1 … E256` — real EGI channel labels —
so the tool classifies them as channels automatically.

## Step 4 — Set the montage and parameters

- **Recording system** → the EGI 256 (HydroCel) option (the legacy label "164 channels" resolves
  to the same montage). The synthetic files use `E1 … E256`, so this must match.
- **Data Type** → `absolute`.
- **Permutations** → `500`. Small on purpose: enough to see the effect, fast enough for a first
  run. Real analyses use 10000 or more.
- **Output** → choose any empty folder you can find again.

## Step 5 — Run it

The **Run** button enables once everything required is set (two files, matching channel count,
output path). Click **Run**.

The tool builds the channel-adjacency structure, computes a two-sample t at every channel, then
repeatedly shuffles the group labels to build a permutation null and correct for testing 256
channels at once. With 500 permutations this takes seconds to a minute.

## Step 6 — See what it found

Open your output folder. You should see three files sharing a timestamped base name:

- `<base>_<timestamp>_report.html` — open this first. It has a **TFCE / Cluster** toggle and
  shows the group-mean topographies and the significance map.
- `<base>_<timestamp>.xlsx` — the significant channels in a table.
- `<base>_<timestamp>.mat` — the full `results_struct` if you want to inspect it in MATLAB.

Look at the topographic map. **A significant cluster of channels appears around E129**, at the
back of the head. That is exactly the effect the synthetic data was built with — a 17-channel
neighbour cluster planted near E129, with noise everywhere else. The tool found the real effect
and (correctly) flagged nothing in the noise.

You've now run a full permutation-based SnPM analysis: loaded data, chosen a design, run a
corrected topographic test, and confirmed it recovers a known effect.

## Where to go next

- Try another synthetic file. The map of every file → analysis → GUI settings is in
  `test_data/synthetic_gui/README.md`, written by `gen_synthetic_testdata`. A good next
  step is **one-way ANOVA** (`glm_anova1.csv`, set Group = `group`), which introduces column
  roles.
- To run any analysis for real, follow the how-to:
  [Run an analysis in the GUI](../how-to/run-an-analysis-in-the-gui.md).
- To pick the right analysis for your own study, read
  [Choosing an analysis](../explanation/choosing-an-analysis.md).
- For the exact inputs and outputs of every analysis, see the
  [Analysis catalog](../reference/ANALYSIS_CATALOG.md).
