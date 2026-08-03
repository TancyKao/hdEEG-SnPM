# How-to guides

*Goal-oriented.* Short, direct recipes for a specific task. They assume you already know what
you want and roughly how the tool works — if you're new, do the
[getting-started tutorial](../tutorials/getting-started.md) first. For *why* a design is
appropriate, see [explanation](../explanation/); for exact inputs/outputs, see the
[reference](../reference/ANALYSIS_CATALOG.md).

## Guides

- **[Run an analysis in the GUI](run-an-analysis-in-the-gui.md)** — drive `SnPMAnalysisGui`
  for CSV files or a spectral folder, map column roles, and read the outputs.
- **[Run a circular (phase) analysis](run-a-circular-phase-analysis.md)** — compare preferred
  phase between two independent groups, or correlate it with a measure: the file layout including
  the required event-count files, the radians requirement, the detector zero-convention table
  (including the 180°-wrong TurtleWave pre-v4 export), the angular-resolution rule, and how to
  read the outputs.
- **[Run a source-level (2447-voxel) analysis](run-a-source-level-analysis.md)** — headless
  permutation t on reconstructed cortical band power: the `src0001..src2447` export contract,
  the magnitude / node-order guards, the significant-voxel deliverable, and the placeholder-
  coordinate precondition.

## Task guides that live next to the code

- **Run analyses headless (NCI / HPC)** — [`scripts/README_scripts.md`](../../scripts/README_scripts.md).
  Edit the `CONFIG` block of a `run_*.m` template and submit with `matlab -batch`. Covers the
  GLM/legacy runner, the LMM runner, the spectral-folder-as-factor runner (with band×stage
  sweep), the event/spectral report runners, and the `db_to_group_table.py` event-prep step.

## Open how-to gaps (not yet written)

- Prepare a TurtleWave event table for the GLM presets (waiting on the importer).
- Add a new recording system to the montage registry (developer task).
