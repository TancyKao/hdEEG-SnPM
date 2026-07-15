# Design-assistant prompt — make the SnPM GUI inputs "smart" (analysis-driven)

Hand this to the design/coding assistant **together with the current GUI file** `SnPMAnalysisGui.m`
(a MATLAB **App Designer** app, `classdef ... < matlab.apps.AppBase`, built programmatically in
`createComponents` with `uigridlayout`). Keep that app as the base and evolve it **additively** — do not
change `core_snpm_analysis.m` or any engine; the GUI only assembles a `params` struct and calls
`runSnPMAnalysis`/`core_snpm_analysis`.

---

> Here is my working MATLAB App Designer GUI `SnPMAnalysisGui.m` (title "hdEEG-SnPM Toolbox"). It has two
> fixed file pickers, **Data 1 File** and **Data 2 File**, plus an analysis dropdown (`ComparisonDropDown`,
> already grouped into plain-language labels via `Items`/`ItemsData`), channel/datatype/tail/permutation
> controls, an output-folder picker, a covariate/column-roles panel, and a results text area.
>
> **The problem:** the two file slots mean different things per analysis and that confuses non-statistician
> users. For the GLM presets only one file is used (Data 2 is dead); for the legacy tests both are used but
> "Data 1/Data 2" doesn't say what they are. Users also have to **type** the metadata-column names and the
> column roles by hand. Make the inputs adapt to the chosen analysis so the user never has to know the
> convention. Keep everything else (engine calls, params keys, grouped dropdown, reports) working.
>
> **Implement these four behaviours, all driven by the selected analysis (`ComparisonDropDown.Value`, which
> is the internal key, e.g. `anova1`, `pairedT`):**
>
> **1. Analysis-driven file pickers (relabel + show/hide).** Replace the static "Data 1/Data 2 File" labels
> with labels that change per analysis, and hide the second picker entirely when it isn't used:
>
> | Analysis key | File 1 label | File 2 label |
> |---|---|---|
> | `anova1`,`ancova`,`regression`,`rmanova`,`mixed2way` | **Data File** | *(hidden)* |
> | `pairedT` | **Condition A File** | **Condition B File** |
> | `onesampleT` | **Condition A File** | **Condition B File** |
> | `unpairedT` | **Group A File** | **Group B File** |
> | `correlationP`,`correlationS` | **Measure 1 File (e.g. EEG)** | **Measure 2 File (e.g. behaviour)** |
> | `circ_wheeler_watson_Test`,`circ_WatsonsU2Test` | **Angles A File (radians)** | **Angles B File (radians)** |
>
> The single-file (GLM) presets must map their one file to `params.data_file`; the two-file tests map to
> `params.data1_file`/`params.data2_file` (as today). The Run button must enable for single-file presets with
> only File 1 + output set (this gate already exists — keep it).
>
> **2. Auto-detect channels vs metadata on file load.** When a file is selected, read its column headers.
> Any column whose name matches `^E\d+$` (or is `Cz`) is a **channel**; everything else is **metadata**.
> Auto-fill the "Meta cols" value with the detected metadata columns (the user can still override). Never
> make the user type the channel list. Show a small caption like "Detected: 256 channels, 4 metadata columns".
>
> **3. Column-role pickers become dropdowns from the file headers.** After a single-file preset's file loads,
> populate the role controls from the detected **metadata** column names instead of free-text:
> - **Group**, **Predictor**, **Condition**, **Subject** → single-select dropdowns of metadata column names.
> - **Covariate cols** → a multi-select list of metadata column names.
> Show **only the roles the selected preset needs** (this mapping is already implemented as `updateRoleFields`
> — keep it; just swap the EditFields for dropdowns/listbox):
>
> | Preset | Roles to show |
> |---|---|
> | `anova1` | Group |
> | `ancova` | Group, Covariates |
> | `regression` | Predictor, Covariates |
> | `rmanova` | Subject, Condition |
> | `mixed2way` | Group, Subject, Condition |
> | legacy/circular tests | none (hide the roles panel) |
>
> **4. Inline guidance + pre-run validation.** Under the analysis dropdown show a one-line plain-language hint
> for the current analysis (e.g. ANOVA → "Compares 3+ groups; pick the column that labels the groups").
> Before enabling Run, validate that every shown role has a selection and that channel count matches the
> montage (256 for "164 channels", 178 for "178 channels"); show a clear inline message if not.
>
> **Keep unchanged:** the grouped plain-language `ComparisonDropDown`; Channels / Data Type / Tails /
> Permutations controls; output-folder picker; results text area; the `params` keys consumed by
> `core_snpm_analysis` (`comparison`, `data_file` or `data1_file`/`data2_file`, `data_sheet`/`data1_sheet`/
> `data2_sheet`, `channels`, `datatype`, `tail`, `permutations`, `meta_cols`, `group_col`, `predictor_col`,
> `condition_col`, `subject_col`, `covariate_cols`); Excel/HTML report generation.
>
> Return the updated `SnPMAnalysisGui.m` (App Designer programmatic style, same property/`createComponents`
> structure). It must still construct headless (`app = SnPMAnalysisGui; delete(app)`) and pass `checkcode`
> with no errors.

---

## Notes for whoever runs the result

- Excel input: when a file is `.xlsx`, the column-header read (steps 2–3) should use the chosen sheet
  (`Data Sheet` dropdown); for `.csv` the sheet is `'CSV File'`.
- Test data to validate the redesigned GUI lives in `test_data/synthetic_gui/` (one file or file-pair per
  analysis, each with a planted cluster near channel E129; see its `README.md`). After wiring, load each and
  confirm the role pickers populate from the file and Run enables correctly.
- The current `updateRoleFields`, `ComparisonDropDownValueChanged`, and `checkReadyToRun` already encode the
  preset→roles map and the single-file Run gate — reuse them; the work is mostly swapping EditFields for
  dropdowns and relabelling/hiding the file pickers on analysis change.
