# Design-assistant prompt — redesign the SnPM GUI layout (spacing + Data Files)

Hand this to the design/coding assistant **together with the current GUI file**
`SnPMAnalysisGui.m` (a MATLAB **App Designer** app, `classdef … < matlab.apps.AppBase`,
built programmatically in `createComponents` with `uifigure` + nested `uigridlayout`,
R2025a). This is a **layout** change only — keep all behaviour and engine calls.

(The earlier "smart inputs" prompt is archived as `DESIGN_PROMPT_gui_v1.md`; that redesign
is already implemented.)

---

> I have a working MATLAB **App Designer** GUI, `SnPMAnalysisGui.m` (title **"hdEEG-SnPM Toolbox"**) — a `classdef … < matlab.apps.AppBase` built **programmatically** in `createComponents` using `uifigure` + nested `uigridlayout` (MATLAB R2025a). It's a front-end for non-parametric EEG statistics; it only assembles a `params` struct and calls `core_snpm_analysis`. **I want you to redesign the layout for clarity and spacing — especially the Data Files area, which is too cramped — without changing any behaviour or the engine calls.**
>
> **Current layout (stacked panels in one `uigridlayout` column):**
> 1. **Choose your analysis** — a plain-language grouped dropdown (`ComparisonDropDown`) + a one-line blue hint box.
> 2. **Data files** — this panel is the problem. It holds a **"Data source: Files | Spectral folder"** toggle, and its contents swap by mode:
>    - *Files mode:* Data 1 / Data 2 file pickers (relabelled per analysis, e.g. "Condition A/B File"), each with a sheet dropdown, plus a "Detected: N channels…" caption.
>    - *Spectral folder mode:* a nested grid with **Add folder… / Remove selected** buttons, a **folder `uitable`** (columns: Folder path, Label — one row per group/condition), two multi-select listboxes (**Bands**, **Stages**), and a **Power type** dropdown.
>    The folder table + listboxes are badly squeezed; everything feels tight.
> 3. **Analysis parameters** — Recording system, Data Type, Tails, Permutations.
> 4. **Column roles** — Group / Predictor / Condition / Subject dropdowns + a Covariate multi-select (shown only for the GLM presets that need them).
> 5. **Output folder**, a **validation checklist** ("Before you can run:"), **Run / Reset / Export** buttons, and a **results text area**.
>
> **What I want:**
> - A cleaner visual hierarchy and **generous, breathing spacing** — the window currently feels packed.
> - The **Data Files panel given much more room**, with the spectral-folder folder-table and the Bands/Stages listboxes sized so they're comfortably usable (the table ~3–5 visible rows; listboxes tall enough to see all bands/stages).
> - Sensible default window size and **resizable** behaviour: use `uigridlayout` row/column weights (`'1x'`, `'fit'`) so panels grow with the window instead of fixed pixel heights that clip. The big Data-files and Results areas should expand; control rows stay compact.
> - Keep it **non-statistician friendly**: clear section titles, the per-analysis hint, and the inline validation.
>
> **Hard constraints (do not break):**
> - Additive only — **do not touch** `core_snpm_analysis.m` or any engine, and keep every `params` key the GUI produces (`comparison`, `data_file` or `data1_file`/`data2_file`, `data_sheet`/`data1_sheet`/`data2_sheet`, `channels`, `datatype`, `tail`, `permutations`, `meta_cols`, `group_col`, `predictor_col`, `condition_col`, `subject_col`, `covariate_cols`, plus the spectral-folder `folders`/`labels` path).
> - Keep all existing behaviour: the grouped `ComparisonDropDown`, the Files↔Spectral-folder toggle and its mode-swap, analysis-driven file relabelling, the folder-table add/remove, multi-select band/stage sweep, the per-preset role pickers, and `checkReadyToRun` validation. Reuse the existing callbacks/handles; this is a **layout** change, not a logic rewrite.
> - Same programmatic App Designer style (properties block + `createComponents` + callbacks). It must still **construct headless** (`app = SnPMAnalysisGui; delete(app)` with no error) and pass `checkcode('SnPMAnalysisGui.m')` with **no errors**.
>
> **Deliverable:** the updated `SnPMAnalysisGui.m` (or, if you'd rather show the design first, an ASCII/HTML mock-up of the new layout and the revised `createComponents` + layout-helper code). Walk me through the grid structure (row/column weights) you chose for the main window and the Data Files panel.

---

## Notes for whoever runs the result

- Attach the **current** `SnPM_2025/SnPMAnalysisGui.m` (~1,400 lines — normal for programmatic
  App Designer). The layout lives in `createComponents`; visibility/sizing logic is in
  `applyDataSourceLayout`, `applyAnalysisLayout`, `updateRoleFields`, `checkReadyToRun`.
- Verify the result headless: `matlab -batch "app=SnPMAnalysisGui; delete(app); disp ok"` and
  `matlab -batch "checkcode('SnPMAnalysisGui.m')"`. Interactive look-and-feel needs a display.
- Test data for a live check is in `example_data/`; spectral-folder mode can point at any two
  local BIDS spectral-power folders (one per condition).
