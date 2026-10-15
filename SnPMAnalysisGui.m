classdef SnPMAnalysisGui < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = private)
        UIFigure                    matlab.ui.Figure
        GridLayout                  matlab.ui.container.GridLayout
        TitleLabel                  matlab.ui.control.Label
        SubtitleLabel               matlab.ui.control.Label
        AnalysisNoteLabel           matlab.ui.control.Label   % "Inputs adapt to the analysis you choose ..."

        % 1. Choose your analysis
        AnalysisPanel               matlab.ui.container.Panel
        ComparisonDropDown          matlab.ui.control.DropDown
        ComparisonLabel             matlab.ui.control.Label
        HintBoxPanel                matlab.ui.container.Panel     % light-blue info box
        HintLabel                   matlab.ui.control.Label       % per-analysis plain-language hint

        % 2. Data files
        DataFilesPanel              matlab.ui.container.Panel
        Data1FileTitleLabel         matlab.ui.control.Label   % "Data File:" / "Condition A File:" etc (relabelled per analysis)
        Data1FileButton             matlab.ui.control.Button
        Data1FileLabel              matlab.ui.control.Label
        Data1SheetDropDown          matlab.ui.control.DropDown
        Data1SheetLabel             matlab.ui.control.Label

        Data2FileTitleLabel         matlab.ui.control.Label
        Data2FileButton             matlab.ui.control.Button
        Data2FileLabel              matlab.ui.control.Label
        Data2SheetTitleLabel        matlab.ui.control.Label
        Data2SheetDropDown          matlab.ui.control.DropDown
        Data2SheetLabel             matlab.ui.control.Label

        % Event-count files (circular phase-group analyses only): one wide
        % subjects x channels file per group, matching the angle file's subject
        % column and channel names. log(event count) enters the model as a
        % covariate -- without it, groups that differ in detected-event count
        % show inflated false positives (fewer events => noisier per-subject
        % phase estimate => genuinely wider observed spread, which the test
        % reads as a real group difference).
        Count1FileTitleLabel        matlab.ui.control.Label
        Count1FileButton            matlab.ui.control.Button
        Count1FileLabel             matlab.ui.control.Label
        Count1SheetTitleLabel       matlab.ui.control.Label
        Count1SheetDropDown         matlab.ui.control.DropDown

        Count2FileTitleLabel        matlab.ui.control.Label
        Count2FileButton            matlab.ui.control.Button
        Count2FileLabel             matlab.ui.control.Label
        Count2SheetTitleLabel       matlab.ui.control.Label
        Count2SheetDropDown         matlab.ui.control.DropDown

        DetectCaptionLabel          matlab.ui.control.Label   % "Detected: N channels, M metadata columns"

        % Data-source toggle + spectral-folder controls (spectral mode)
        DataSourceLabel             matlab.ui.control.Label
        DataSourceDropDown          matlab.ui.control.DropDown    % 'Files' | 'Spectral folder'
        SpecGrid                    matlab.ui.container.GridLayout % nested grid holding the spectral controls
        SpecAddFolderButton         matlab.ui.control.Button
        SpecRemoveFolderButton      matlab.ui.control.Button
        SpecFolderTable             matlab.ui.control.Table        % folders (path + level label), one per factor level
        SpecBandTitleLabel          matlab.ui.control.Label
        SpecStageTitleLabel         matlab.ui.control.Label
        SpecTypeTitleLabel          matlab.ui.control.Label
        SpecBandListBox             matlab.ui.control.ListBox     % multi-select bands (sweep)
        SpecStageListBox            matlab.ui.control.ListBox     % multi-select stages (sweep)
        SpecTypeListBox             matlab.ui.control.ListBox      % multi-select power types (absolute / normalized)

        % 3. Analysis Parameters Panel
        ParametersPanel             matlab.ui.container.Panel
        ChannelsDropDown            matlab.ui.control.DropDown
        ChannelsLabel               matlab.ui.control.Label
        DataTypeDropDown            matlab.ui.control.DropDown
        DataTypeLabel               matlab.ui.control.Label
        TailsDropDown               matlab.ui.control.DropDown
        TailsLabel                  matlab.ui.control.Label
        PermutationField            matlab.ui.control.NumericEditField
        PermutationLabel            matlab.ui.control.Label

        % Column-role controls for GLM presets (anova1/ancova/regression/rmanova/mixed2way).
        % Role pickers are dropdowns/listbox populated from the loaded file's metadata columns.
        RolesPanel                  matlab.ui.container.Panel
        MetaColsField               matlab.ui.control.EditField   % auto-filled detected metadata cols (editable override)
        GroupColField               matlab.ui.control.DropDown
        PredictorColField           matlab.ui.control.DropDown
        ConditionColField           matlab.ui.control.DropDown
        SubjectColField             matlab.ui.control.DropDown
        CovColsField                matlab.ui.control.ListBox     % multi-select covariate columns

        % Phase / angle settings Panel (circular analyses only).
        % Units and zero-phase convention are REQUIRED and deliberately have no
        % preselected value: the engine rejects a missing value rather than
        % guessing, because guessing wrong silently rotates every angle.
        CircPanel                   matlab.ui.container.Panel
        CircUnitsLabel              matlab.ui.control.Label
        CircUnitsDropDown           matlab.ui.control.DropDown
        CircConventionLabel         matlab.ui.control.Label
        CircConventionDropDown      matlab.ui.control.DropDown
        CircOffsetLabel             matlab.ui.control.Label
        CircOffsetField             matlab.ui.control.NumericEditField
        CircNoteLabel               matlab.ui.control.Label

        % Covariate-file Panel (legacy correlation covariates)
        CovariatePanel              matlab.ui.container.Panel
        CovariateFileButton         matlab.ui.control.Button
        CovariateFileLabel          matlab.ui.control.Label
        UseCovariatesCheckBox       matlab.ui.control.CheckBox

        % 4. Output Panel
        OutputPanel                 matlab.ui.container.Panel
        OutputPathButton            matlab.ui.control.Button
        OutputPathLabel             matlab.ui.control.Label

        % Pre-run validation checklist ("Before you can run:")
        ValidationPanel             matlab.ui.container.Panel
        ValidationLabel             matlab.ui.control.Label

        % Control Buttons
        RunAnalysisButton           matlab.ui.control.Button
        ResetButton                 matlab.ui.control.Button
        ExportPDFButton             matlab.ui.control.Button

        % Status
        StatusLabel                 matlab.ui.control.Label

        % Results Panel
        ResultsPanel                matlab.ui.container.Panel
        ResultsTextArea             matlab.ui.control.TextArea
    end

    % Properties for storing file paths and data
    properties (Access = private)
        data1_file = ''
        data2_file = ''
        count1_file = ''  % group A event-count file (circular phase-group analyses)
        count2_file = ''  % group B event-count file (circular phase-group analyses)
        output_path = ''
        covariate_file = ''  % Store covariate file path
        snpm_toolbox_path = ''
        eeglab_path = ''
        last_results = ''  % Store last results for PDF export
        lastComparison = 'pairedT'  % last real (non-header) comparison key selected
        RoleRows = {}  % {key, labelHandle, fieldHandle} for show/hide of GLM column-role fields
        % Header-detection state, kept PER SLOT. The primary file (Data 1, or
        % the first spectral folder) owns the montage validation and the role
        % pickers; Data 2 carries its own counts. They used to share one pair of
        % properties, so loading Data 2 overwrote what Data 1 had detected and
        % the montage check then validated the wrong file -- or, for a Data 2
        % holding no channel columns at all, silently stopped validating.
        DetectedNChan = 0  % channel columns in the PRIMARY file (montage validation)
        DetectedChanNames = {}  % channel column NAMES in the PRIMARY file (label-presence validation)
        DetectCaption1 = ''  % "Detected: ..." caption text for the primary file
        DetectedNChan2 = 0  % channel columns in the Data 2 file (0 = none, e.g. a behavioural table)
        % snpm_detect_behavior_table results per slot: .kind montage|behavior|
        % unknown, .channel_cols, .data_cols (every non-ID column, ID = the
        % FIRST 'subject' column only), .col_name. [] until a file is read;
        % always go through normalizeDetect() rather than touching them raw.
        DetectB1 = []
        DetectB2 = []
        SpectralFolder = ''  % selected BIDS spectral folder (spectral mode)
        SubjectsCsv = ''  % optional subject-metadata CSV joined on Subject (spectral mode)
        FullComparisonItems = {}  % saved ComparisonDropDown Items/ItemsData to restore when leaving spectral mode
        FullComparisonItemsData = {}
        % Channel-column counts per picker, used to catch an angle/count file
        % pair that does not describe the same channels.
        NChanAngle1 = 0
        NChanAngle2 = 0
        NChanCount1 = 0
        NChanCount2 = 0
        % Data Type / Tails selections saved when a circular analysis locks
        % them, so switching back to any other analysis restores the user's
        % choice rather than silently dropping the options.
        PreCircDataType = 'absolute'
        PreCircTail = 'both'
        % mixed2way balance check, cached. The key is file + sheet + the
        % Subject/Condition/Group roles, so the table is read only when one of
        % those changes, not on every readiness re-check. A new file pick
        % (detectAndFill on the primary slot) clears the key so a re-loaded,
        % edited file is read again.
        BalanceKey = ''
        BalanceMsg = {}
    end

    % Row indices of the collapsible auxiliary panels in app.GridLayout.
    % They sit AFTER the numbered sections (1-4 are always shown, in order).
    properties (Access = private, Constant)
        RolesRow = 8
        CovRow   = 9
        CircRow  = 10   % phase/angle settings (circular analyses only)
        RolesRowH = 170
        CovRowH   = 64
        % Rows of dataFilesGrid holding the event-count file pickers.
        CountRow1 = 4
        CountRow2 = 5
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % Initialize the app
            app.StatusLabel.Text = 'Ready to start analysis';
            app.RunAnalysisButton.Enable = 'off';

            % Set default values
            app.ChannelsDropDown.Value = 'egi';
            app.DataTypeDropDown.Value = 'absolute';
            app.ComparisonDropDown.Value = 'pairedT';
            app.TailsDropDown.Value = 'both';
            app.PermutationField.Value = 5000;

            % Clear sheet dropdowns
            app.Data1SheetDropDown.Items = {};
            app.Data2SheetDropDown.Items = {};

            app.lastComparison = 'pairedT';
            initRolePickers(app);    % empty role dropdowns until a file is loaded
            updateRoleFields(app);   % show only the column-role controls the analysis needs
            applyAnalysisLayout(app);% relabel file pickers + hint for the default analysis
            applyDataSourceLayout(app); % start in Files mode; collapse the spectral controls
            checkReadyToRun(app);    % populate the "Before you can run:" checklist
        end

        % Comparison changed: ignore header rows; refresh layout, role fields, hint, gate
        function ComparisonDropDownValueChanged(app, ~)
            val = app.ComparisonDropDown.Value;
            if ischar(val) && startsWith(val, '__h')   % inert header row
                app.ComparisonDropDown.Value = app.lastComparison;
                return
            end
            app.lastComparison = val;
            updateRoleFields(app);
            applyAnalysisLayout(app);
            checkReadyToRun(app);
        end

        % Re-detect columns when the Excel sheet changes
        function Data1SheetChanged(app, ~)
            if ~isempty(app.data1_file)
                [~, n] = detectAndFill(app, app.data1_file, app.Data1SheetDropDown.Value, true);
                app.NChanAngle1 = n;
                checkReadyToRun(app);
            end
        end
        function Data2SheetChanged(app, ~)
            if ~isempty(app.data2_file)
                [~, n, ~, B] = detectAndFill(app, app.data2_file, app.Data2SheetDropDown.Value, false);
                setData2Detection(app, n, B);
                app.NChanAngle2 = n;
                checkReadyToRun(app);
            end
        end
        function ChannelsChanged(app, ~)
            checkReadyToRun(app);   % montage changed -> re-validate channel count
        end

        % ---- Spectral-folder data source ----
        function DataSourceDropDownValueChanged(app, ~)
            applyDataSourceLayout(app);
            checkReadyToRun(app);
        end

        % Show file-source rows or the spectral nested grid; resize the panel;
        % restrict the comparison list to the spectral-supported analyses.
        function applyDataSourceLayout(app)
            spec = strcmp(app.DataSourceDropDown.Value, 'Spectral folder');
            g  = app.SpecGrid.Parent;                  % dataFilesGrid
            rh = g.RowHeight;
            % dataFilesGrid rows: 1 source toggle, 2 file-1, 3 file-2,
            % 4 count-1, 5 count-2, 6 spectral grid, 7 caption.
            % The count rows are sized by applyAnalysisLayout (circular only);
            % here they are only forced shut for the spectral source.
            if spec
                rh{2} = 0;     rh{3} = 0;     rh{6} = '1x';   % hide file rows; let the spectral grid grow
                rh{app.CountRow1} = 0; rh{app.CountRow2} = 0;
            else
                rh{2} = 'fit'; rh{3} = 'fit'; rh{6} = 0;      % file rows fit content; collapse spectral grid
            end
            g.RowHeight = rh;
            app.SpecGrid.Visible = spec;
            for h = {app.Data1FileTitleLabel, app.Data1FileButton, app.Data1FileLabel, app.Data1SheetDropDown, ...
                     app.Data2FileTitleLabel, app.Data2FileButton, app.Data2FileLabel, app.Data2SheetTitleLabel, app.Data2SheetDropDown}
                h{1}.Visible = ~spec;
            end
            orh = app.GridLayout.RowHeight;
            if spec, orh{5} = '2x'; else, orh{5} = 'fit'; end   % Data-files area expands in spectral mode
            app.GridLayout.RowHeight = orh;
            restrictComparisons(app, spec);
            if spec   % roles are automatic in spectral mode (factor = folder label)
                app.RolesPanel.Visible = false; app.CovariatePanel.Visible = false;
                app.CircPanel.Visible = false;
                orh = app.GridLayout.RowHeight;
                orh{app.RolesRow} = 0; orh{app.CovRow} = 0; orh{app.CircRow} = 0;
                app.GridLayout.RowHeight = orh;
            end
        end

        % In spectral mode show only the folder-as-factor analyses; restore the
        % full grouped list (saved once) when switching back to Files.
        function restrictComparisons(app, spec)
            if spec
                if isempty(app.FullComparisonItems)
                    app.FullComparisonItems     = app.ComparisonDropDown.Items;
                    app.FullComparisonItemsData = app.ComparisonDropDown.ItemsData;
                end
                app.ComparisonDropDown.Items     = {'Two independent samples t-test (2 groups)', 'Paired samples t-test (2 conditions)', ...
                    'One-sample t-test (vs 0)', 'One-way ANOVA (3+ groups)', 'Repeated-measures ANOVA (3+ conditions, within)'};
                app.ComparisonDropDown.ItemsData = {'unpairedT', 'pairedT', 'onesampleT', 'anova1', 'rmanova'};
                if ~ismember(app.ComparisonDropDown.Value, app.ComparisonDropDown.ItemsData)
                    app.ComparisonDropDown.Value = 'pairedT';
                    app.lastComparison = 'pairedT';
                end
            elseif ~isempty(app.FullComparisonItems)
                app.ComparisonDropDown.Items     = app.FullComparisonItems;
                app.ComparisonDropDown.ItemsData = app.FullComparisonItemsData;
                if ismember(app.lastComparison, app.FullComparisonItemsData)
                    app.ComparisonDropDown.Value = app.lastComparison;
                end
            end
            updateRoleFields(app);
            applyAnalysisLayout(app);
        end

        % Add a spectral-power folder (one level of the design factor). The first
        % folder fills the band/stage/type pickers + montage channel detection.
        function SpecAddFolderButtonPushed(app, ~)
            d = uigetdir(pwd, 'Select a spectral-power folder (one group / condition)');
            if isequal(d, 0), return; end
            if isempty(app.SpecFolderTable.Data)
                try
                    info = scan_spectral_folder(d);
                catch ME
                    uialert(app.UIFigure, ME.message, 'Could not read spectral folder'); return;
                end
                app.SpecBandListBox.Items  = info.bands;   app.SpecBandListBox.Value  = info.bands(1);
                app.SpecStageListBox.Items = info.stages;  app.SpecStageListBox.Value = info.stages(1);
                app.DetectedNChan     = numel(info.chanlabels);
                app.DetectedChanNames = info.chanlabels;
                app.DetectCaption1    = sprintf('Detected: %d channels, %d bands, %d stages', ...
                    numel(info.chanlabels), numel(info.bands), numel(info.stages));
                app.DetectCaptionLabel.Text = app.DetectCaption1;
            end
            [~, base] = fileparts(d);
            app.SpecFolderTable.Data = [app.SpecFolderTable.Data; {d, base}];
            checkReadyToRun(app);
        end

        function SpecRemoveFolderButtonPushed(app, ~)
            sel = app.SpecFolderTable.Selection;
            D = app.SpecFolderTable.Data;
            if isempty(sel) && ~isempty(D), sel = size(D,1); end   % default: last row
            D(sel, :) = [];
            app.SpecFolderTable.Data = D;
            checkReadyToRun(app);
        end

        function SpecSelectionChanged(app, ~)
            checkReadyToRun(app);
        end

        % Loop the selected bands x stages, run each cell, collect a results grid.
        % Folders + labels come from the folder table; the factor is the label.
        function runSpectralSweep(app)
            comparison = app.ComparisonDropDown.Value;
            types   = cellstr(app.SpecTypeListBox.Value);
            bands   = cellstr(app.SpecBandListBox.Value);
            stages  = cellstr(app.SpecStageListBox.Value);
            outRoot = app.output_path;
            D = app.SpecFolderTable.Data;
            FOLDERS = D(:,1)'; LABELS = D(:,2)';

            base = struct('channels',app.ChannelsDropDown.Value, ...
                'datatype',app.DataTypeDropDown.Value, 'permutations',app.PermutationField.Value, ...
                'tail',app.TailsDropDown.Value);

            grid = {}; lines = {};
            for ti = 1:numel(types)
                type = types{ti};
                for bi = 1:numel(bands)
                    for si = 1:numel(stages)
                        stg = stages{si};
                        celldir = fullfile(outRoot, sprintf('%s_%s_%s_%s', comparison, type, bands{bi}, stg));
                        if ~exist(celldir,'dir'), mkdir(celldir); end
                        o = base; o.type = type; o.band = bands{bi}; o.stages = {stg}; o.output_path = celldir;
                        o.folders = FOLDERS; o.labels = LABELS;
                        app.StatusLabel.Text = sprintf('Running %s | %s | %s / %s ...', comparison, type, bands{bi}, stg); drawnow;
                        [params, ns] = spectral_to_snpm_params('', comparison, o);
                        rs = core_snpm_analysis(params);
                        u=numel(rs.uncorrsigch); t=numel(rs.correctTFCEsigch); c=numel(rs.SnPMsigch); p=specMinP(rs);
                        grid(end+1,:) = {comparison, bands{bi}, type, stg, ns, u, t, c, p}; %#ok<AGROW>
                        lines{end+1} = sprintf('%-9s %-10s %-9s %-7s  n=%d  TFCE=%d  cluster=%d  minP=%.3f', ...
                            comparison, type, bands{bi}, stg, ns, t, c, p); %#ok<AGROW>
                    end
                end
            end
            G = cell2table(grid, 'VariableNames', {'comparison','band','type','stage','n','nUncorr','nTFCE','nCluster','minClusterP'});
            gf = fullfile(outRoot, sprintf('SWEEP_grid_%s.csv', comparison));
            writetable(G, gf);
            app.last_results = [{sprintf('Spectral sweep complete: %d cell(s)', size(grid,1)); ['Grid CSV: ' gf]; ''}; lines(:)];
            app.ResultsTextArea.Value = app.last_results;

            % ---- faceted dashboard: one HTML (abs+relative x all stages x Uncorrected/TFCE/Cluster) ----
            % 2-level contrasts only; runs its own full sweep via export_report (all bands/stages/powers).
            dashMsg = '';
            if ismember(comparison, {'pairedT','onesampleT','unpairedT'})
                try
                    app.StatusLabel.Text = 'Building faceted dashboard (abs/relative x all stages)...'; drawnow;
                    root = fileparts(mfilename('fullpath'));
                    dash = fullfile(outRoot, 'dashboard'); if ~exist(dash,'dir'), mkdir(dash); end
                    copyfile(fullfile(root,'templates','sleep_eeg_report.html'), fullfile(dash,'sleep_eeg_report.html'));
                    dopts = struct('comparison',comparison, 'permutations',app.PermutationField.Value, ...
                        'powers',{{'absolute','normalised'}}, 'condA',LABELS{1}, 'condB',LABELS{2});
                    dopts.folders = FOLDERS; dopts.labels = LABELS;
                    export_report('', dash, dopts);
                    dashMsg = fullfile(dash, 'sleep_eeg_report_filled.html');
                    app.last_results = [app.last_results; {''; ['Dashboard: ' dashMsg]}];
                    app.ResultsTextArea.Value = app.last_results;
                catch ME
                    app.last_results = [app.last_results; {''; ['Dashboard failed: ' ME.message]}];
                    app.ResultsTextArea.Value = app.last_results;
                end
            end
            if isempty(dashMsg)
                app.StatusLabel.Text = sprintf('Spectral sweep complete: %d cell(s).', size(grid,1));
            else
                app.StatusLabel.Text = sprintf('Done: %d cell(s) + faceted dashboard.', size(grid,1));
            end
        end

        % Relabel the file pickers, hide the 2nd picker for single-file presets, set the hint,
        % and show only the panels the selected analysis needs.
        function applyAnalysisLayout(app)
            key = app.ComparisonDropDown.Value;
            [l1, l2, hint, c1, c2] = analysis_labels(key);
            app.Data1FileTitleLabel.Text = l1;
            app.HintLabel.Text = hint;
            twoFile = ~isempty(l2);
            if twoFile, app.Data2FileTitleLabel.Text = l2; end
            for h = {app.Data2FileTitleLabel, app.Data2FileButton, app.Data2FileLabel, ...
                     app.Data2SheetTitleLabel, app.Data2SheetDropDown}
                h{1}.Visible = twoFile;
            end

            % Roles panel only for the GLM presets; covariate-file panel only for correlation.
            isGLM  = ismember(key, {'anova1','ancova','regression','rmanova','mixed2way'});
            isCorr = ischar(key) && contains(key, 'correlation');
            isCirc = ischar(key) && startsWith(key, 'circ');
            app.RolesPanel.Visible     = isGLM;
            app.CovariatePanel.Visible = isCorr;
            app.CircPanel.Visible      = isCirc;

            % ---- Event-count file pickers (circular phase-group analyses) ----
            % Required whenever the design compares phase distributions between
            % groups; not used by the circular-linear correlation, whose pairing
            % is by subject. Collapsed and hidden otherwise.
            spec = strcmp(app.DataSourceDropDown.Value, 'Spectral folder');
            needCounts = isCirc && ~isempty(c1) && ~spec;
            if needCounts
                app.Count1FileTitleLabel.Text = c1;
                app.Count2FileTitleLabel.Text = c2;
            end
            for h = {app.Count1FileTitleLabel, app.Count1FileButton, app.Count1FileLabel, ...
                     app.Count1SheetTitleLabel, app.Count1SheetDropDown, ...
                     app.Count2FileTitleLabel, app.Count2FileButton, app.Count2FileLabel, ...
                     app.Count2SheetTitleLabel, app.Count2SheetDropDown}
                h{1}.Visible = needCounts;
            end
            drh = app.SpecGrid.Parent.RowHeight;
            if needCounts
                drh{app.CountRow1} = 'fit'; drh{app.CountRow2} = 'fit';
            else
                drh{app.CountRow1} = 0;     drh{app.CountRow2} = 0;
            end
            app.SpecGrid.Parent.RowHeight = drh;

            % ---- Data Type / Tails: locked for circular analyses ----
            % Angles are a circular domain: log10 of a radian goes complex and a
            % cross-channel z-score is meaningless, so only 'absolute' is valid.
            % Both circular statistics are non-negative and inherently
            % upper-tailed, so only 'both' is valid; the engine rejects the rest.
            % Switching AWAY must restore the full option sets, otherwise every
            % other analysis in the toolbox silently loses those options.
            if isCirc
                if strcmp(char(app.DataTypeDropDown.Enable), 'on')
                    app.PreCircDataType = app.DataTypeDropDown.Value;   % save once, on entry
                    app.PreCircTail     = app.TailsDropDown.Value;
                end
                app.DataTypeDropDown.Items  = {'absolute'};
                app.DataTypeDropDown.Value  = 'absolute';
                app.DataTypeDropDown.Enable = 'off';
                app.DataTypeDropDown.Tooltip = ...
                    'Angles are used as-is; log/z-score are undefined on a circular domain.';
                app.TailsDropDown.Items  = {'both'};
                app.TailsDropDown.Value  = 'both';
                app.TailsDropDown.Enable = 'off';
                app.TailsDropDown.Tooltip = ...
                    'The circular statistics are non-negative and inherently upper-tailed.';
            else
                app.DataTypeDropDown.Items = {'absolute', 'logscale', 'normalize'};
                if ismember(app.PreCircDataType, app.DataTypeDropDown.Items)
                    app.DataTypeDropDown.Value = app.PreCircDataType;
                else
                    app.DataTypeDropDown.Value = 'absolute';
                end
                app.DataTypeDropDown.Enable  = 'on';
                app.DataTypeDropDown.Tooltip = '';
                app.TailsDropDown.Items = {'both', 'left', 'right'};
                if ismember(app.PreCircTail, app.TailsDropDown.Items)
                    app.TailsDropDown.Value = app.PreCircTail;
                else
                    app.TailsDropDown.Value = 'both';
                end
                app.TailsDropDown.Enable  = 'on';
                app.TailsDropDown.Tooltip = '';
            end

            % Recording-system picker is analysis-aware: source-space montages
            % (is_source) have no GLM guards in core_snpm_glm, so exclude them
            % for the GLM presets. They are excluded for circular analyses too:
            % no exporter produces a per-voxel preferred phase, so the option
            % would only ever be a mistake. If a source system was selected when
            % the user switches to such an analysis, fall back to the scalp
            % default (egi) and note it in the hint so the reset is visible.
            reg = snpm_montage_registry();
            if isGLM || isCirc, keep = ~[reg.is_source]; else, keep = true(1, numel(reg)); end
            prevSys = app.ChannelsDropDown.Value;
            app.ChannelsDropDown.Items     = {reg(keep).display};
            app.ChannelsDropDown.ItemsData = {reg(keep).key};
            if ismember(prevSys, app.ChannelsDropDown.ItemsData)
                app.ChannelsDropDown.Value = prevSys;   % keep the user's choice
            else
                app.ChannelsDropDown.Value = 'egi';     % scalp default
                app.HintLabel.Text = [hint ...
                    '   (Source montage is scalp-only for this analysis ' char(8212) ' reset to EGI 256.)'];
            end

            % Collapse the rows of the hidden panels so there is no empty gap;
            % a visible auxiliary panel sizes to its content ('fit').
            rh = app.GridLayout.RowHeight;
            if isGLM,  rh{app.RolesRow} = 'fit'; else, rh{app.RolesRow} = 0; end
            if isCorr, rh{app.CovRow}   = 'fit'; else, rh{app.CovRow}   = 0; end
            if isCirc, rh{app.CircRow}  = 'fit'; else, rh{app.CircRow}  = 0; end
            app.GridLayout.RowHeight = rh;
        end

        % Enable the custom zero-offset field only for the 'custom' convention.
        function CircConventionChanged(app, ~)
            isCustom = strcmp(app.CircConventionDropDown.Value, 'custom');
            app.CircOffsetField.Enable = matlab.lang.OnOffSwitchState(isCustom);
            app.CircOffsetLabel.Enable = matlab.lang.OnOffSwitchState(isCustom);
            checkReadyToRun(app);
        end

        % Show only the GLM column-role controls relevant to the selected analysis
        function updateRoleFields(app)
            needed = struct( ...
                'anova1',     {{'meta','group'}}, ...
                'ancova',     {{'meta','group','cov'}}, ...
                'regression', {{'meta','predictor','cov'}}, ...
                'rmanova',    {{'meta','subject','condition'}}, ...
                'mixed2way',  {{'meta','group','subject','condition'}});
            key = app.ComparisonDropDown.Value;
            if isfield(needed, key), want = needed.(key); else, want = {}; end  % legacy tests: hide all
            for r = 1:size(app.RoleRows,1)
                vis = ismember(app.RoleRows{r,1}, want);
                app.RoleRows{r,2}.Visible = vis;   % label
                app.RoleRows{r,3}.Visible = vis;   % control
            end
        end

        % Initialise role pickers to an empty "(load file first)" state
        function initRolePickers(app)
            for dd = {app.GroupColField, app.PredictorColField, app.ConditionColField, app.SubjectColField}
                dd{1}.Items = {'(load file first)'};
                dd{1}.ItemsData = {''};
                dd{1}.Value = '';
            end
            app.CovColsField.Items = {};
            app.CovColsField.Value = {};
            app.MetaColsField.Value = '';
            app.DetectCaptionLabel.Text = '';
            app.DetectedNChan = 0;
            app.DetectedChanNames = {};
            app.DetectB1 = [];
            app.DetectCaption1 = '';
            app.BalanceKey = '';
            app.BalanceMsg = {};
            clearData2Detection(app);
        end

        % Forget everything detected in the Data 2 slot.
        function clearData2Detection(app)
            app.DetectedNChan2 = 0;
            app.DetectB2 = [];
        end

        % Record what detectAndFill found in the Data 2 slot (see the property
        % block: Data 2 must not overwrite the primary file's detection).
        function setData2Detection(app, nchan, B)
            app.DetectedNChan2 = nchan;
            app.DetectB2 = B;
        end

        % True when Data 2's columns line up with the ones core_snpm_analysis
        % will actually read out of it. snpm_setup_channels matches column
        % positions on DATA 1's header and applies those same positions to
        % Data 2 (snpm_setup_channels.m:49-65) -- it never reads Data 2's own
        % names. This mirrors the engine's channelOrderMismatch guard exactly
        % (the Data 2 alignment block, core_snpm_analysis.m:245-456; the order
        % guard itself at :326-425): the same data-column lists (ID column
        % = the FIRST 'subject' match only), the same expected labels (the
        % montage's LOADED chanlocs, registry labels only as a fallback), the
        % same short-file padding, and strcmpi -- a pure case difference
        % ('e1' vs 'E1') is not a mis-pairing and must not block a run the
        % engine accepts. The engine runs that block for EVERY two-file legacy
        % analysis, not just the correlations, so this helper is used for the
        % t-tests too.
        function tf = data2ColumnsAlign(app)
            tf = true;
            B1 = normalizeDetect(app.DetectB1);
            B2 = normalizeDetect(app.DetectB2);
            c1 = B1.data_cols; c2 = B2.data_cols;
            if isempty(c1) || isempty(c2), return; end
            pos = find(ismember(c1, montageOrderLabels(app.ChannelsDropDown.Value)));
            if isempty(pos), return; end          % nothing selected: not our error
            c2p = repmat({''}, 1, max(pos));      % pad: a short file mismatches, never errors
            n_avail = min(numel(c2), max(pos));
            c2p(1:n_avail) = c2(1:n_avail);
            tf = all(strcmpi(c1(pos), c2p(pos)));
        end

        % Read a file's headers and classify channel vs metadata columns.
        % Returns the split for the caller's slot; the shared detection state
        % (montage validation) and the role pickers are written ONLY for the
        % primary file, so a Data 2 / event-count pick cannot clobber them.
        %   META  non-channel column names   NCHAN number of channel columns
        %   CHANS channel column names       B     the classifier result
        function [meta, nchan, chans, B] = detectAndFill(app, file, sheet, isPrimary)
            % Read the header EXACTLY as the engine will. core_snpm_analysis
            % loads both tables with a plain readtable (core_snpm_analysis.m:116
            % and :136), whose default 'modify' naming rule rewrites any header
            % that is not a valid MATLAB identifier -- 'MRRT (ms)' becomes
            % 'MRRT_ms_'. detectImportOptions defaults to that same rule, so
            % NOT asking for 'preserve' is what keeps these names identical to
            % readtable's (verified equal for CSV and .xlsx). Asking for
            % 'preserve' broke three things at once: the classifier judged a
            % rewritten header against its raw spelling, data2ColumnsAlign
            % compared Data 1's raw names with labels the engine sees modified,
            % and every role name this function writes into the role pickers
            % (and MetaColsField) went to the engine in a spelling readtable
            % had already changed, so the column could not be found.
            try
                if strcmp(sheet, 'CSV File')
                    o = detectImportOptions(file);
                else
                    o = detectImportOptions(file, 'Sheet', sheet);
                end
                vn = o.VariableNames;
            catch
                vn = {};
            end
            % Classification is delegated to the shared helper
            % dependencies/snpm_detect_behavior_table.m -- the same function the
            % engine uses to decide what a correlation Data 2 is -- so given the
            % same column names (see above) the GUI's idea of "this is a
            % channel", of which column is the ID, and of
            % montage-vs-behaviour-vs-unknown cannot drift from the core's.
            try
                B = snpm_detect_behavior_table(vn);
            catch ME
                % snpm:behavior:noDataColumns is the only error the helper can
                % raise on a well-formed cellstr (multipleColumns is RETIRED --
                % the ambiguous case now comes back as kind 'unknown'). It means
                % an unreadable header or a file holding nothing but an ID
                % column. A file pick must not throw: keep the message and let
                % the readiness gate report it where it matters.
                if ~strcmp(ME.identifier, 'snpm:behavior:noDataColumns')
                    rethrow(ME);
                end
                B = normalizeDetect([]);
                B.error_message = ME.message;
            end
            chans = B.channel_cols;
            meta  = vn(~ismember(vn, chans));
            nchan = numel(chans);
            if isPrimary
                app.BalanceKey = '';    % file (re)read: force the mixed2way balance re-check
                app.DetectedNChan = nchan;
                app.DetectedChanNames = chans;
                app.DetectB1 = B;
                app.DetectCaption1 = sprintf('Detected: %d channels, %d metadata columns', nchan, numel(meta));
                app.MetaColsField.Value = strjoin(meta, ',');
                fillRoleDropdown(app.GroupColField,     meta, firstMatch(meta,'group'));
                fillRoleDropdown(app.PredictorColField, meta, '');
                fillRoleDropdown(app.ConditionColField, meta, firstMatch(meta,'cond'));
                fillRoleDropdown(app.SubjectColField,   meta, firstMatch(meta,'subject'));
                app.CovColsField.Items = meta;
                app.CovColsField.Value = {};
            end
        end

        % Button pushed function: Data1FileButton
        function Data1FileButtonPushed(app, event)
            [file, path] = uigetfile({'*.xlsx;*.xls;*.csv', 'Excel and CSV Files (*.xlsx,*.xls,*.csv)'}, ...
                                   'Select Data 1 File');
            drawnow;
            figure(app.UIFigure);
            if file ~= 0
                app.data1_file = fullfile(path, file);
                app.Data1FileLabel.Text = file;
                app.StatusLabel.Text = 'Loading Data 1 file...';

                % Load sheet names for Excel files
                if endsWith(file, {'.xlsx', '.xls'})
                    try
                        sheets = sheetnames(app.data1_file);
                        app.Data1SheetDropDown.Items = sheets;
                        app.Data1SheetDropDown.Value = sheets{1};
                        app.StatusLabel.Text = 'Data 1 file loaded successfully';
                    catch ME
                        app.StatusLabel.Text = 'Error reading Excel file';
                        uialert(app.UIFigure, ME.message, 'File Error');
                    end
                else
                    % For CSV files, no sheet selection needed
                    app.Data1SheetDropDown.Items = {'CSV File'};
                    app.Data1SheetDropDown.Value = 'CSV File';
                    app.StatusLabel.Text = 'Data 1 CSV file loaded successfully';
                end

                [~, n] = detectAndFill(app, app.data1_file, app.Data1SheetDropDown.Value, true);
                app.NChanAngle1 = n;
                checkReadyToRun(app);
            end
        end

        % Button pushed function: Data2FileButton
        function Data2FileButtonPushed(app, event)
            [file, path] = uigetfile({'*.xlsx;*.xls;*.csv', 'Excel and CSV Files (*.xlsx,*.xls,*.csv)'}, ...
                                   'Select Data 2 File');
            drawnow;
            figure(app.UIFigure);
            if file ~= 0
                app.data2_file = fullfile(path, file);
                app.Data2FileLabel.Text = file;
                app.StatusLabel.Text = 'Loading Data 2 file...';

                % Load sheet names for Excel files
                if endsWith(file, {'.xlsx', '.xls'})
                    try
                        sheets = sheetnames(app.data2_file);
                        app.Data2SheetDropDown.Items = sheets;
                        app.Data2SheetDropDown.Value = sheets{1};
                        app.StatusLabel.Text = 'Data 2 file loaded successfully';
                    catch ME
                        app.StatusLabel.Text = 'Error reading Excel file';
                        uialert(app.UIFigure, ME.message, 'File Error');
                    end
                else
                    % For CSV files, no sheet selection needed
                    app.Data2SheetDropDown.Items = {'CSV File'};
                    app.Data2SheetDropDown.Value = 'CSV File';
                    app.StatusLabel.Text = 'Data 2 CSV file loaded successfully';
                end

                [~, n, ~, B] = detectAndFill(app, app.data2_file, app.Data2SheetDropDown.Value, false);
                setData2Detection(app, n, B);
                app.NChanAngle2 = n;
                checkReadyToRun(app);
            end
        end

        % Button pushed function: Count1FileButton (group A event counts)
        function Count1FileButtonPushed(app, event)
            pickCountFile(app, 1);
        end

        % Button pushed function: Count2FileButton (group B event counts)
        function Count2FileButtonPushed(app, event)
            pickCountFile(app, 2);
        end

        function Count1SheetChanged(app, ~)
            if ~isempty(app.count1_file)
                [~, n] = detectAndFill(app, app.count1_file, app.Count1SheetDropDown.Value, false);
                app.NChanCount1 = n;
                checkReadyToRun(app);
            end
        end
        function Count2SheetChanged(app, ~)
            if ~isempty(app.count2_file)
                [~, n] = detectAndFill(app, app.count2_file, app.Count2SheetDropDown.Value, false);
                app.NChanCount2 = n;
                checkReadyToRun(app);
            end
        end

        % Shared picker for the two event-count files (mirrors the Data 1 / 2
        % pickers: browse, list Excel sheets, detect channel columns, revalidate).
        function pickCountFile(app, which)
            [file, path] = uigetfile({'*.xlsx;*.xls;*.csv', 'Excel and CSV Files (*.xlsx,*.xls,*.csv)'}, ...
                                   sprintf('Select the Group %c Event Count File', 'A' + which - 1));
            drawnow;
            figure(app.UIFigure);
            if isequal(file, 0), return; end
            full = fullfile(path, file);
            if which == 1
                dd = app.Count1SheetDropDown; lb = app.Count1FileLabel;
            else
                dd = app.Count2SheetDropDown; lb = app.Count2FileLabel;
            end
            if endsWith(file, {'.xlsx', '.xls'})
                try
                    sheets = sheetnames(full);
                    dd.Items = sheets; dd.Value = sheets{1};
                catch ME
                    % leave the slot empty rather than half-set, so the gate
                    % keeps reporting the file as missing
                    app.StatusLabel.Text = 'Error reading Excel file';
                    uialert(app.UIFigure, ME.message, 'File Error');
                    return;
                end
            else
                dd.Items = {'CSV File'}; dd.Value = 'CSV File';
            end
            if which == 1, app.count1_file = full; else, app.count2_file = full; end
            lb.Text = file;
            [~, n] = detectAndFill(app, full, dd.Value, false);
            if which == 1, app.NChanCount1 = n; else, app.NChanCount2 = n; end
            app.StatusLabel.Text = sprintf('Group %c event-count file loaded', 'A' + which - 1);
            checkReadyToRun(app);
        end

        % Button pushed function: OutputPathButton
        function OutputPathButtonPushed(app, event)
            folder = uigetdir('', 'Select Output Directory');
            drawnow;
            figure(app.UIFigure);
            if folder ~= 0
                app.output_path = folder;
                app.OutputPathLabel.Text = folder;
                app.StatusLabel.Text = 'Output path selected';
                checkReadyToRun(app);
            end
        end

        % Button pushed function: CovariateFileButton
        function CovariateFileButtonPushed(app, event)
            [file, path] = uigetfile('*.csv', 'Select Covariate CSV File');
            drawnow;
            figure(app.UIFigure);
            if file ~= 0
                app.covariate_file = fullfile(path, file);
                app.CovariateFileLabel.Text = file;
                app.StatusLabel.Text = 'Loading covariate file...';

                % Validate CSV file format
                try
                    covariate_data = readtable(app.covariate_file);
                    if size(covariate_data, 2) < 2
                        uialert(app.UIFigure, 'CSV file must have at least 2 columns (Subject ID and covariates)', 'Invalid File Format');
                        app.covariate_file = '';
                        app.CovariateFileLabel.Text = 'No covariate file selected';
                        app.StatusLabel.Text = 'Ready to start analysis';
                        return;
                    end

                    % Display preview of covariate columns
                    col_names = covariate_data.Properties.VariableNames;
                    msg = sprintf('Detected %d columns:\n%s\n\nFirst column should be Subject ID.\nRemaining columns will be used as covariates.', ...
                        length(col_names), strjoin(col_names, ', '));
                    uialert(app.UIFigure, msg, 'Covariate File Loaded', 'Icon', 'info');

                    app.StatusLabel.Text = sprintf('Covariate file loaded: %d subjects, %d covariates', ...
                        size(covariate_data, 1), size(covariate_data, 2) - 1);

                    % Auto-enable the checkbox
                    app.UseCovariatesCheckBox.Value = true;

                catch ME
                    uialert(app.UIFigure, ['Error reading CSV file: ' ME.message], 'File Error');
                    app.covariate_file = '';
                    app.CovariateFileLabel.Text = 'No covariate file selected';
                    app.StatusLabel.Text = 'Ready to start analysis';
                end
            end
        end

        % Button pushed function: RunAnalysisButton
        function RunAnalysisButtonPushed(app, event)
            % Disable the run button to prevent multiple runs
            app.RunAnalysisButton.Enable = 'off';
            app.StatusLabel.Text = 'Running analysis...';

            % Clear previous results
            app.ResultsTextArea.Value = '';

            try
                % Spectral-folder data source: loop the selected bands x stages
                if strcmp(app.DataSourceDropDown.Value, 'Spectral folder')
                    runSpectralSweep(app);
                    app.RunAnalysisButton.Enable = 'on';
                    return;
                end

                % Prepare parameters
                params = struct();
                params.data1_file = app.data1_file;
                params.data2_file = app.data2_file;
                params.data1_sheet = app.Data1SheetDropDown.Value;
                params.data2_sheet = app.Data2SheetDropDown.Value;
                params.output_path = app.output_path;
                params.channels= app.ChannelsDropDown.Value;
                params.datatype = app.DataTypeDropDown.Value;
                params.comparison = app.ComparisonDropDown.Value;
                params.tail = app.TailsDropDown.Value;
                params.permutations = app.PermutationField.Value;
                params.covariate_file = app.covariate_file;
                params.use_covariates = app.UseCovariatesCheckBox.Value;

                % ---- GLM presets: single Data-1 table + column roles ----
                glmPresets = {'anova1','ancova','regression','rmanova','mixed2way'};
                if ismember(params.comparison, glmPresets)
                    params.data_file  = app.data1_file;
                    params.data_sheet = app.Data1SheetDropDown.Value;
                    params.meta_cols  = splitCols(app.MetaColsField.Value);
                    if isempty(params.meta_cols)
                        uialert(app.UIFigure,'Enter the non-channel "Meta cols" (e.g. Subject,group,...) for this analysis.','Meta columns required');
                        app.RunAnalysisButton.Enable='on'; app.StatusLabel.Text='Ready to start analysis'; return;
                    end
                    if ~isempty(app.GroupColField.Value),     params.group_col     = app.GroupColField.Value; end
                    if ~isempty(app.PredictorColField.Value), params.predictor_col = app.PredictorColField.Value; end
                    if ~isempty(app.ConditionColField.Value), params.condition_col = app.ConditionColField.Value; end
                    if ~isempty(app.SubjectColField.Value),   params.subject_col   = app.SubjectColField.Value; end
                    cc = app.CovColsField.Value; if ischar(cc), cc = {cc}; end
                    cc = cc(~cellfun(@isempty, cc)); if ~isempty(cc), params.covariate_cols = cc; end
                    app.StatusLabel.Text = 'Running GLM analysis...'; drawnow;
                    results = runSnPMAnalysis(app, params);
                    app.StatusLabel.Text = 'Analysis completed successfully!';
                    app.ResultsTextArea.Value = results; app.last_results = results;
                    app.RunAnalysisButton.Enable = 'on';
                    return;
                end

                % ---- Circular (phase / angle) analyses ----
                % Units and convention are always passed; the engine rejects a
                % missing value rather than assuming one. Event-count files are
                % passed for the phase-group designs only.
                if ischar(params.comparison) && startsWith(params.comparison, 'circ')
                    params.circ_units      = app.CircUnitsDropDown.Value;
                    params.circ_convention = app.CircConventionDropDown.Value;
                    if strcmp(params.circ_convention, 'custom')
                        params.circ_zero_offset_deg = app.CircOffsetField.Value;
                    end
                    [~, ~, ~, cl1] = analysis_labels(params.comparison);
                    if ~isempty(cl1)
                        params.count1_file  = app.count1_file;
                        params.count2_file  = app.count2_file;
                        params.count1_sheet = app.Count1SheetDropDown.Value;
                        params.count2_sheet = app.Count2SheetDropDown.Value;
                    end
                end

                if contains(params.comparison, 'correlation')
                    app.StatusLabel.Text = 'Validating subject matching for correlation...';
                    drawnow;

                    % Quick check if both files have subject columns
                    try
                        if strcmp(params.data1_sheet, 'CSV File')
                            temp1 = readtable(params.data1_file, 'ReadVariableNames', true);
                        else
                            temp1 = readtable(params.data1_file, 'Sheet', params.data1_sheet, 'ReadVariableNames', true);
                        end

                        if strcmp(params.data2_sheet, 'CSV File')
                            temp2 = readtable(params.data2_file, 'ReadVariableNames', true);
                        else
                            temp2 = readtable(params.data2_file, 'Sheet', params.data2_sheet, 'ReadVariableNames', true);
                        end

                        has_subj1 = any(contains(lower(temp1.Properties.VariableNames), 'subject'));
                        has_subj2 = any(contains(lower(temp2.Properties.VariableNames), 'subject'));

                        if ~has_subj1 || ~has_subj2
                            uialert(app.UIFigure, 'Correlation analysis requires Subject columns in both files!', 'Missing Subject Data');
                            app.RunAnalysisButton.Enable = 'on';
                            app.StatusLabel.Text = 'Ready to start analysis';
                            return;
                        end

                        % Quick preview of subject matching
                        subj1 = temp1{:, contains(lower(temp1.Properties.VariableNames), 'subject')};
                        subj2 = temp2{:, contains(lower(temp2.Properties.VariableNames), 'subject')};

                        %if isnumeric(subj1), subj1 = string(subj1); else, subj1 = string(subj1); end
                        %if isnumeric(subj2), subj2 = string(subj2); else, subj2 = string(subj2); end

                        common_subj = intersect(subj1, subj2);

                        if length(common_subj) < 3
                            uialert(app.UIFigure, sprintf('Correlation analysis requires at least 3 matching subjects. Found: %d', length(common_subj)), 'Insufficient Matching Subjects');
                            app.RunAnalysisButton.Enable = 'on';
                            app.StatusLabel.Text = 'Ready to start analysis';
                            return;
                        end

                        if length(common_subj) < length(subj1) || length(common_subj) < length(subj2)
                            % Show warning but continue
                            choice = uiconfirm(app.UIFigure, ...
                                sprintf('Subject mismatch detected!\nData1: %d subjects, Data2: %d subjects\nMatching: %d subjects\n\nContinue with matching subjects only?', ...
                                        length(subj1), length(subj2), length(common_subj)), ...
                                'Subject Mismatch Warning', ...
                                'Options', {'Continue', 'Cancel'}, ...
                                'DefaultOption', 1, ...
                                'CancelOption', 2);

                            if strcmp(choice, 'Cancel')
                                app.RunAnalysisButton.Enable = 'on';
                                app.StatusLabel.Text = 'Ready to start analysis';
                                return;
                            end
                        end

                    catch ME
                        uialert(app.UIFigure, ['Error checking subject data: ' ME.message], 'Validation Error');
                        app.RunAnalysisButton.Enable = 'on';
                        app.StatusLabel.Text = 'Ready to start analysis';
                        return;
                    end
                end

                app.StatusLabel.Text = 'Processing data...';
                drawnow;

                % Call the main analysis function
                results = runSnPMAnalysis(app, params);

                app.StatusLabel.Text = 'Analysis completed successfully!';

                % Display results
                app.ResultsTextArea.Value = results;
                app.last_results = results;

            catch ME
                app.StatusLabel.Text = 'Analysis failed - see error details';
                if strcmp(ME.identifier, 'core_snpm_glm:mixed2wayUnbalanced')
                    % A data problem, not a code fault: the engine's message
                    % names the subjects and the fix, so no stack trace.
                    app.ResultsTextArea.Value = {['Error: ' ME.message]};
                    app.BalanceKey = '';   % file may have changed since it was checked
                else
                    app.ResultsTextArea.Value = {['Error: ' ME.message], '', 'Stack trace:', ME.getReport()};
                end
                uialert(app.UIFigure, ME.message, 'Analysis Error');
            end

            % Re-enable the run button. After an unbalanced-table refusal the
            % readiness gate re-reads the table (cache cleared above) and puts
            % the incomplete subjects in the checklist, keeping Run disabled.
            app.RunAnalysisButton.Enable = 'on';
            if isempty(app.BalanceKey) && strcmp(app.ComparisonDropDown.Value, 'mixed2way')
                checkReadyToRun(app);
            end
        end

        % Button pushed function: ResetButton
        function ResetButtonPushed(app, event)
            % Reset all fields to default
            app.data1_file = '';
            app.data2_file = '';
            app.count1_file = '';
            app.count2_file = '';
            app.output_path = '';
            app.covariate_file = '';

            app.Data1FileLabel.Text = 'No file selected';
            app.Data2FileLabel.Text = 'No file selected';
            app.Count1FileLabel.Text = 'No file selected';
            app.Count2FileLabel.Text = 'No file selected';
            app.OutputPathLabel.Text = 'No path selected';
            app.CovariateFileLabel.Text = 'No covariate file selected';
            app.UseCovariatesCheckBox.Value = false;

            % clear GLM column-role controls + detection state
            initRolePickers(app);

            app.Data1SheetDropDown.Items = {};
            app.Data2SheetDropDown.Items = {};
            app.Count1SheetDropDown.Items = {};
            app.Count2SheetDropDown.Items = {};

            % circular settings back to "not chosen" (they are required inputs)
            app.CircUnitsDropDown.Value = '';
            app.CircConventionDropDown.Value = '';
            app.CircOffsetField.Value = 0;
            app.CircOffsetField.Enable = 'off';
            app.CircOffsetLabel.Enable = 'off';
            app.NChanAngle1 = 0; app.NChanAngle2 = 0;
            app.NChanCount1 = 0; app.NChanCount2 = 0;

            app.ChannelsDropDown.Value = 'egi';
            app.PreCircDataType = 'absolute';
            app.PreCircTail = 'both';
            app.ComparisonDropDown.Value = 'pairedT';
            app.lastComparison = 'pairedT';
            updateRoleFields(app);
            applyAnalysisLayout(app);   % restores + resets Data Type / Tails
            app.PermutationField.Value = 5000;

            app.StatusLabel.Text = 'Ready to start analysis';
            app.ResultsTextArea.Value = '';
            app.last_results = '';

            checkReadyToRun(app);
        end

        % Button pushed function: ExportPDFButton
        function ExportPDFButtonPushed(app, event)
            if isempty(app.last_results)
                uialert(app.UIFigure, 'No results to export. Please run analysis first.', 'Export Error');
                return;
            end

            % The reports are already generated during analysis
            % Just inform the user
            uialert(app.UIFigure, 'HTML reports have been generated in the output directory!', 'Reports Ready', 'Icon', 'info');
        end

        % Helper function to export results to PDF
        function exportResultsToPDF(app, filepath)
            % Create a temporary figure for PDF export
            temp_fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [0 0 8.5 11]);

            % Create text annotation with results
            results_str = strjoin(app.last_results, '\n');
            % Fix underscore issue by replacing underscores with escaped underscores
            results_str = strrep(results_str, '_', '\_');

            % Add title
            annotation(temp_fig, 'textbox', [0.1 0.9 0.8 0.08], ...
                'String', 'SnPM Analysis Results', ...
                'FontSize', 16, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'EdgeColor', 'none', ...
                'Interpreter', 'tex');  % Added interpreter specification

            % Add results text
            annotation(temp_fig, 'textbox', [0.1 0.1 0.8 0.75], ...
                'String', results_str, ...
                'FontSize', 10, 'FontName', 'Courier', ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'left', ...
                'EdgeColor', 'black', 'LineWidth', 1, ...
                'Interpreter', 'tex');  % Added interpreter specification

            % Export to PDF
            exportgraphics(temp_fig, filepath, 'ContentType', 'vector');

            % Close temporary figure
            close(temp_fig);
        end

        % Helper function to check if ready to run analysis (analysis-aware).
        % Drives the "Before you can run:" checklist and the Run-button gate.
        function checkReadyToRun(app)
            key = app.ComparisonDropDown.Value;
            isGLM = ismember(key, {'anova1','ancova','regression','rmanova','mixed2way'});
            isCirc = ischar(key) && startsWith(key, 'circ');
            [~, l2, ~, c1] = analysis_labels(key);
            twoFile = ~isempty(l2);
            needCounts = isCirc && ~isempty(c1);

            missing = {};
            spec = strcmp(app.DataSourceDropDown.Value, 'Spectral folder');

            if spec
                D = app.SpecFolderTable.Data;
                nF = size(D,1);
                need2 = any(strcmp(key, {'pairedT','onesampleT','unpairedT'}));
                if nF == 0
                    missing{end+1} = 'Add the spectral folders (one per group/condition)';
                elseif need2 && nF ~= 2
                    missing{end+1} = sprintf('%s needs exactly 2 folders (have %d)', key, nF);
                elseif any(strcmp(key,{'anova1','rmanova'})) && nF < 2
                    missing{end+1} = sprintf('%s needs at least 2 folders (have %d)', key, nF);
                end
                if nF > 0
                    labs = D(:,2);
                    if any(cellfun(@(x) isempty(x) || (ischar(x)&&isempty(strtrim(x))), labs))
                        missing{end+1} = 'Give every folder a label';
                    elseif numel(unique(labs)) < numel(labs)
                        missing{end+1} = 'Folder labels must be unique';
                    end
                end
                if isempty(app.SpecBandListBox.Value),   missing{end+1} = 'Select at least one band'; end
                if isempty(app.SpecStageListBox.Value),  missing{end+1} = 'Select at least one stage'; end
                if isempty(app.SpecTypeListBox.Value),   missing{end+1} = 'Select at least one power type'; end
            else
                % File 1 (always needed)
                lab1 = strtrim(erase(app.Data1FileTitleLabel.Text, ':'));
                if isempty(app.data1_file), missing{end+1} = ['Select the ' lab1]; end

                % File 2 (two-file legacy/circular/correlation tests)
                if twoFile
                    lab2 = strtrim(erase(app.Data2FileTitleLabel.Text, ':'));
                    if isempty(app.data2_file), missing{end+1} = ['Select the ' lab2]; end
                end

                % Event-count files: hard requirement for the circular
                % phase-group analyses. The engine errors without them, so
                % gating here is the difference between a clear checklist item
                % and a run that dies after the user clicks Run.
                if needCounts
                    labc1 = strtrim(erase(app.Count1FileTitleLabel.Text, ':'));
                    labc2 = strtrim(erase(app.Count2FileTitleLabel.Text, ':'));
                    if isempty(app.count1_file), missing{end+1} = ['Select the ' labc1]; end
                    if isempty(app.count2_file), missing{end+1} = ['Select the ' labc2]; end
                    % An angle/count pair describing different channel sets is a
                    % wrong-file mistake, not a formatting one -- catch it here.
                    if app.NChanAngle1 > 0 && app.NChanCount1 > 0 && app.NChanAngle1 ~= app.NChanCount1
                        missing{end+1} = sprintf(['Group A count file has %d channel columns but the ' ...
                            'angle file has %d ' char(8212) ' they must describe the same channels'], ...
                            app.NChanCount1, app.NChanAngle1);
                    end
                    if app.NChanAngle2 > 0 && app.NChanCount2 > 0 && app.NChanAngle2 ~= app.NChanCount2
                        missing{end+1} = sprintf(['Group B count file has %d channel columns but the ' ...
                            'angle file has %d ' char(8212) ' they must describe the same channels'], ...
                            app.NChanCount2, app.NChanAngle2);
                    end
                end
            end

            % Angle units + zero-phase convention: required, no default. Both
            % silently rotate or rescale every angle if wrong, and neither is
            % recoverable from the output, so neither may be left unset.
            if isCirc
                if isempty(app.CircUnitsDropDown.Value)
                    missing{end+1} = 'Choose the angle units (radians or degrees)';
                end
                if isempty(app.CircConventionDropDown.Value)
                    missing{end+1} = 'Choose where zero sits in your exported phase (the convention)';
                end
            end

            % Output folder
            if isempty(app.output_path), missing{end+1} = 'Choose an output folder'; end

            % The PRIMARY file must contain channels of the selected recording
            % system. Keyed on the primary slot only: Data 2 may legitimately
            % hold no channels at all (a correlation against one behavioural
            % score), and an event-count file is checked against its angle file
            % above, not against the montage.
            capLines = {};
            if ~isempty(app.DetectCaption1), capLines{end+1} = app.DetectCaption1; end
            if app.DetectedNChan > 0
                m = snpm_montage_registry(app.ChannelsDropDown.Value);
                present = nnz(ismember(m.labels, app.DetectedChanNames));
                if present == 0
                    missing{end+1} = sprintf('No %s channels found in the file (expected labels like %s)', ...
                        m.display, m.labels{min(2,numel(m.labels))});
                else
                    capLines{1} = sprintf('Detected: %d channels  (matches %d of %d %s channels)', ...
                        app.DetectedNChan, present, numel(m.labels), m.display);
                end
            end

            % Data 2 has three recognised shapes, so it gets its own check.
            % This block mirrors the engine's Data 2 alignment block branch for
            % branch (core_snpm_analysis.m:245-456) -- what the engine accepts,
            % warns about and refuses -- so the checklist never blocks a run the
            % engine would happily do, nor waves through one it would kill after
            % the user clicked Run. It switches on the SHARED classifier's
            % .kind, not on a guess of its own:
            %   'behavior'  Subject + ONE measure column: one value per subject,
            %               broadcast across Data 1's channels. CORRELATION ONLY
            %               -- a t-test compares two channel maps, so there is
            %               nothing for one score per subject to pair with and
            %               the engine refuses it (snpm:data2:behaviorNotSupported,
            %               core_snpm_analysis.m:272-283).
            %   'montage'   names recognised as channel labels: matched BY NAME.
            %               Same channel count, then same channels in the same
            %               order (see data2ColumnsAlign).
            %   'unknown'   too few names recognised to match by name -- what a
            %               lowercased ('e1') or suffixed ('E1_density') montage
            %               export looks like. The engine matches it BY POSITION
            %               and warns, so an equal data-column count is ALLOWED
            %               here with an unverified-pairing caption; only a
            %               different count is an error.
            %   'none'      headers unreadable, or nothing but an ID column.
            %
            % The engine guards EVERY two-file legacy analysis this way, not
            % just the correlations: the block is gated on ~onesample_vs_zero
            % (core_snpm_analysis.m:245), which is true for pairedT, unpairedT,
            % both correlations, and a onesampleT that was given a Data 2. So
            % the gate here is the same set, keyed on a Data 2 actually being
            % loaded -- app.data2_file is passed to the engine unconditionally,
            % so a file left over from a previous analysis is guarded too,
            % exactly as the engine would guard it. Before this, a mis-ordered
            % pairedT pair showed ready and only failed after Run.
            if ~spec && ismember(key, {'correlationP','correlationS','pairedT','unpairedT','onesampleT'}) ...
                    && ~isempty(app.data2_file)
                isCorr = ismember(key, {'correlationP','correlationS'});
                if twoFile
                    lab2 = strtrim(erase(app.Data2FileTitleLabel.Text, ':'));
                else
                    % onesampleT hides the Data 2 picker but nothing clears
                    % app.data2_file, and it is still handed to the engine --
                    % which then runs the legacy paired-difference-vs-0 path on
                    % it instead of the vs-zero path. Naming the hidden picker
                    % would send the user looking for a control that is not on
                    % screen, so name the situation and the way out.
                    lab2 = ['second data file left over from a previous analysis ' ...
                            '(press Reset to clear it)'];
                end
                % The behavioural format is a CORRELATION format. Offering it
                % as a fix for a t-test width mismatch would send the user
                % straight into snpm:data2:behaviorNotSupported.
                if isCorr
                    bothFormats = ['export both files on the same montage, or give Data 2 just a ' ...
                                   'Subject column and ONE behavioural measure column'];
                else
                    bothFormats = 'export both files on the same montage, one column per channel.';
                end
                B1 = normalizeDetect(app.DetectB1);
                B2 = normalizeDetect(app.DetectB2);
                n1_chan = numel(B1.channel_cols); n1_data = numel(B1.data_cols);
                n2_chan = numel(B2.channel_cols); n2_data = numel(B2.data_cols);
                switch B2.kind
                    case 'behavior'
                        if isCorr
                            capLines{end+1} = sprintf(['Behavioral measure: %s (1 value per subject, ' ...
                                'compared against every channel of Data 1)'], B2.col_name);
                        else
                            % One score per subject can only be broadcast across
                            % channels for a correlation. This analysis compares
                            % two channel maps, so there is nothing to pair it
                            % with; the engine refuses the run outright
                            % (snpm:data2:behaviorNotSupported).
                            missing{end+1} = sprintf(['%s is a Subject column plus ONE measure ' ...
                                'column (''%s''). A single score per subject is only valid for a ' ...
                                'correlation ' char(8212) ' this analysis compares two channel ' ...
                                'maps, so Data 2 needs one column per channel. Export Data 2 on ' ...
                                'the same montage as Data 1, or change the analysis to ' ...
                                'Correlation (Pearson) or Correlation (Spearman).'], ...
                                lab2, B2.col_name);
                        end

                    case 'montage'
                        if n1_chan > 0 && n1_chan ~= n2_chan
                            missing{end+1} = sprintf(['%s has %d channel columns but Data 1 has %d ' ...
                                char(8212) ' the two files are paired channel by channel (channel i ' ...
                                'of Data 1 with channel i of Data 2), so %s'], ...
                                lab2, n2_chan, n1_chan, bothFormats);
                        elseif n1_chan > 0 && n2_chan > 0 && ~data2ColumnsAlign(app)
                            % Same COUNT is not same CHANNELS. Channel selection
                            % is matched on Data 1's header and the same column
                            % positions are applied to Data 2, so a different
                            % order means channel i is compared against some
                            % other channel and the whole map is wrong with
                            % nothing in the output to show it. The engine raises
                            % snpm:correlation:channelOrderMismatch (kept under
                            % that identifier, but raised for every two-file
                            % legacy analysis) for headless runs; this is the
                            % interactive half of that guard.
                            missing{end+1} = sprintf(['%s has the same number of channels as Data 1 ' ...
                                'but not in the same order ' char(8212) ' the two files are read ' ...
                                'position by position, so channel i of Data 1 would be paired with ' ...
                                'a different channel of Data 2 and the whole map would be wrong ' ...
                                'with nothing in the output to show it. Re-export Data 2 with its ' ...
                                'columns in Data 1''s order.'], lab2);
                        else
                            capLines{end+1} = sprintf('Data 2: %d channel columns, matched by name', n2_chan);
                        end

                    case 'unknown'
                        if n1_data > 0 && n1_data ~= n2_data
                            missing{end+1} = sprintf(['%s has %d data columns but Data 1 has %d, and ' ...
                                'Data 2''s columns are not recognised channel labels so they cannot ' ...
                                'be matched by name either ' char(8212) ' %s'], ...
                                lab2, n2_data, n1_data, bothFormats);
                        else
                            capLines{end+1} = ['Data 2 matched by position (column names not ' ...
                                'recognised ' char(8212) ' unverified pairing). Check both files ' ...
                                'list their channels in the same order.'];
                        end

                    otherwise   % 'none': unreadable, or an ID column and nothing else
                        if ~isempty(B2.error_message)
                            missing{end+1} = B2.error_message;
                        else
                            missing{end+1} = sprintf('Could not read the columns of the %s', lab2);
                        end
                end
            end
            if isempty(capLines)
                app.DetectCaptionLabel.Text = '';
            else
                app.DetectCaptionLabel.Text = capLines;   % one line per loaded slot
            end

            % Source-space montage is scalp-only for the GLM presets (no source
            % guards in core_snpm_glm) and for the circular analyses (no
            % exporter produces a per-voxel preferred phase) -- backstop the
            % dropdown filter in case a stale value survives.
            if isGLM || isCirc
                msys = snpm_montage_registry(app.ChannelsDropDown.Value);
                if isfield(msys,'is_source') && msys.is_source
                    missing{end+1} = sprintf(['%s is scalp-only ' char(8212) ...
                        ' the source montage is not supported for this analysis'], key);
                end
            end

            % Every visible role picker (except covariates, which are optional) must be set
            if isGLM
                for r = 1:size(app.RoleRows,1)
                    rkey = app.RoleRows{r,1}; h = app.RoleRows{r,3};
                    if any(strcmp(rkey, {'meta','cov'})), continue; end
                    if strcmp(char(h.Visible),'on') && isempty(h.Value)
                        lab = strtrim(erase(app.RoleRows{r,2}.Text, ':'));
                        missing{end+1} = ['Choose the ' lab];
                    end
                end
            end

            % mixed2way needs a balanced table (one row per subject per
            % condition, one group per subject); the engine refuses anything
            % else (core_snpm_glm:mixed2wayUnbalanced). Checked only once the
            % file and the Subject + Condition roles are set; cached.
            if strcmp(key, 'mixed2way') && ~spec && ~isempty(app.data1_file) ...
                    && ~isempty(app.SubjectColField.Value) && ~isempty(app.ConditionColField.Value)
                missing = [missing, mixed2wayBalanceMessages(app)];
            end

            if isempty(missing)
                app.RunAnalysisButton.Enable = 'on';
                app.ValidationPanel.Title = 'Ready to run';
                app.ValidationPanel.ForegroundColor = [0.15 0.55 0.15];
                app.ValidationLabel.FontColor = [0.15 0.55 0.15];
                app.ValidationLabel.Text = 'All required inputs are set. Click Run Analysis.';
            else
                app.RunAnalysisButton.Enable = 'off';
                app.ValidationPanel.Title = 'Before you can run:';
                app.ValidationPanel.ForegroundColor = [0.80 0.20 0.20];
                app.ValidationLabel.FontColor = [0.80 0.20 0.20];
                bullet = char(8226);   % real Unicode bullet (single-quoted escapes are literal in MATLAB)
                app.ValidationLabel.Text = cellfun(@(s) [bullet '  ' s], missing, 'UniformOutput', false);
            end
        end

        % Cached mixed2way balance check (see BalanceKey). Returns a cell of
        % checklist lines, empty when balanced or when the table could not be
        % read (a read failure must not block Run; the engine then reports).
        function msgs = mixed2wayBalanceMessages(app)
            key = strjoin({char(app.data1_file), char(app.Data1SheetDropDown.Value), ...
                char(app.SubjectColField.Value), char(app.ConditionColField.Value), ...
                char(app.GroupColField.Value)}, char(0));
            if ~strcmp(key, app.BalanceKey)
                app.BalanceMsg = mixed2way_balance_messages(app.data1_file, ...
                    app.Data1SheetDropDown.Value, app.SubjectColField.Value, ...
                    app.ConditionColField.Value, app.GroupColField.Value);
                app.BalanceKey = key;
            end
            msgs = app.BalanceMsg;
        end

        % Main analysis function that calls your existing code
        function results = runSnPMAnalysis(app, params)
            try
                % Call the extracted core analysis function (without progress callback)
                [results_struct, results_text] = core_snpm_analysis(params);

                % Return formatted results
                results = results_text;

            catch ME
                rethrow(ME);
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % dependencies/ holds the montage registry + label lists used to
            % build the recording-system dropdown and validate files.
            addpath(fullfile(fileparts(mfilename('fullpath')), 'dependencies'));

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.94 0.94 0.94];
            app.UIFigure.Position = [100 80 880 780];   % laptop-friendly; scroll covers the rest
            app.UIFigure.Name = 'hdEEG-SnPM Toolbox';
            app.UIFigure.Icon = 'none';
            app.UIFigure.WindowStyle = 'normal';

            % Create GridLayout (single stacked column of numbered sections)
            app.GridLayout = uigridlayout(app.UIFigure);
            % Scrollable so the full form (taller than a laptop screen) is always
            % reachable: a vertical scrollbar appears whenever the summed row
            % heights exceed the window. Requires determinate row heights (no
            % top-level '1x'), so the Results row below is a fixed height.
            app.GridLayout.Scrollable = 'on';
            app.GridLayout.ColumnWidth = {'1x'};
            % rows: 1 title, 2 subtitle, 3 note, 4 analysis(1), 5 DATA(2, expands),
            %       6 params(3), 7 output(4), 8 roles(GLM), 9 covariate(corr),
            %       10 phase/angle(circular), 11 validation, 12 buttons,
            %       13 status, 14 results.
            % Sections 1-4 always render in order; roles/covariate/phase are
            % auxiliary boxes below them, start collapsed (0), and are sized by
            % the layout helpers (see RolesRow / CovRow / CircRow).
            % Weighted rows so panels grow with the window instead of clipping.
            % 'fit' = size to content (compact control rows); '1x' / '2x' = expand.
            % Last row (Results) is a fixed height, not '1x': a '1x' row would
            % expand to fill and suppress overflow, so the grid would never
            % scroll. Fixed height keeps the total content height determinate.
            app.GridLayout.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 0, 0, 0, 'fit', 'fit', 'fit', 180};
            app.GridLayout.RowSpacing = 10;
            app.GridLayout.Padding = [16 16 16 16];

            % Create Title
            app.TitleLabel = uilabel(app.GridLayout);
            app.TitleLabel.Layout.Row = 1;
            app.TitleLabel.Layout.Column = 1;
            app.TitleLabel.FontSize = 20;
            app.TitleLabel.FontWeight = 'bold';
            app.TitleLabel.FontColor = [0.2 0.4 0.8];
            app.TitleLabel.HorizontalAlignment = 'center';
            app.TitleLabel.Text = 'hdEEG-SnPM Toolbox';

            % Create Subtitle
            app.SubtitleLabel = uilabel(app.GridLayout);
            app.SubtitleLabel.Layout.Row = 2;
            app.SubtitleLabel.Layout.Column = 1;
            app.SubtitleLabel.FontSize = 12;
            app.SubtitleLabel.FontAngle = 'italic';
            app.SubtitleLabel.HorizontalAlignment = 'center';
            app.SubtitleLabel.Text = 'Non-parametric (SnPM) topographic statistics for high-density EEG: spectral, event & GLM/LMM analyses';

            % Adaptive-inputs note
            app.AnalysisNoteLabel = uilabel(app.GridLayout);
            app.AnalysisNoteLabel.Layout.Row = 3;
            app.AnalysisNoteLabel.Layout.Column = 1;
            app.AnalysisNoteLabel.FontSize = 11;
            app.AnalysisNoteLabel.FontColor = [0.2 0.4 0.8];
            app.AnalysisNoteLabel.HorizontalAlignment = 'center';
            app.AnalysisNoteLabel.Text = 'Inputs adapt to the analysis you choose - no need to remember the file convention.';

            % ============ 1. Choose your analysis ============
            app.AnalysisPanel = uipanel(app.GridLayout);
            app.AnalysisPanel.Layout.Row = 4;
            app.AnalysisPanel.Layout.Column = 1;
            app.AnalysisPanel.Title = '1.  Choose your analysis';
            app.AnalysisPanel.FontWeight = 'bold';
            app.AnalysisPanel.FontSize = 12;

            analysisGrid = uigridlayout(app.AnalysisPanel);
            analysisGrid.ColumnWidth = {120, '1x'};
            analysisGrid.RowHeight = {30, 46};
            analysisGrid.RowSpacing = 8;
            analysisGrid.Padding = [10 10 10 10];

            app.ComparisonLabel = uilabel(analysisGrid);
            app.ComparisonLabel.Layout.Row = 1;
            app.ComparisonLabel.Layout.Column = 1;
            app.ComparisonLabel.Text = 'Comparison:';

            app.ComparisonDropDown = uidropdown(analysisGrid);
            app.ComparisonDropDown.Layout.Row = 1;
            app.ComparisonDropDown.Layout.Column = 2;
            % Plain-language labels (Items) mapped to internal keys (ItemsData).
            % Header rows have '__h*__' data and are inert (reverted in the callback).
            app.ComparisonDropDown.Items = { ...
                '--- Compare 2 groups ---', ...
                '   Two independent samples t-test (2 groups)', ...
                '   Paired samples t-test (2 conditions)', ...
                '   One-sample t-test (vs 0)', ...
                '--- Compare 3+ groups / conditions ---', ...
                '   One-way ANOVA (3+ groups)', ...
                '   One-way ANCOVA (3+ groups + covariates)', ...
                '   Repeated-measures ANOVA (3+ conditions, within)', ...
                '   Two-way mixed ANOVA (group × condition)', ...
                '--- Relate to a measure ---', ...
                '   Pearson correlation (linear)', ...
                '   Spearman correlation (rank)', ...
                '   Linear regression (continuous predictor)', ...
                '--- Phase / angle data (advanced) ---', ...
                '   Preferred phase: 2 groups (recommended)', ...
                '   Preferred phase: 2 groups (Watson U2, alternative)', ...
                '   Preferred phase vs a linear measure (paired by subject)'};
            app.ComparisonDropDown.ItemsData = { ...
                '__h1__', 'unpairedT', 'pairedT', 'onesampleT', ...
                '__h2__', 'anova1', 'ancova', 'rmanova', 'mixed2way', ...
                '__h3__', 'correlationP', 'correlationS', 'regression', ...
                '__h4__', 'circ_phase_group', 'circ_phase_group_u2', 'circ_corrAngLinear'};
            app.ComparisonDropDown.ValueChangedFcn = createCallbackFcn(app, @ComparisonDropDownValueChanged, true);

            % Light-blue hint box with a per-analysis one-line explanation
            app.HintBoxPanel = uipanel(analysisGrid);
            app.HintBoxPanel.Layout.Row = 2;
            app.HintBoxPanel.Layout.Column = [1 2];
            app.HintBoxPanel.BackgroundColor = [0.90 0.94 1.00];
            app.HintBoxPanel.BorderType = 'line';

            hintGrid = uigridlayout(app.HintBoxPanel);
            hintGrid.ColumnWidth = {'1x'};
            hintGrid.RowHeight = {'1x'};
            hintGrid.Padding = [8 4 8 4];

            app.HintLabel = uilabel(hintGrid);
            app.HintLabel.Layout.Row = 1;
            app.HintLabel.Layout.Column = 1;
            app.HintLabel.WordWrap = 'on';
            app.HintLabel.VerticalAlignment = 'center';
            app.HintLabel.FontColor = [0.10 0.30 0.55];
            app.HintLabel.Text = '';

            % ============ 2. Data files ============
            app.DataFilesPanel = uipanel(app.GridLayout);
            app.DataFilesPanel.Layout.Row = 5;
            app.DataFilesPanel.Layout.Column = 1;
            app.DataFilesPanel.Title = '2.  Data files';
            app.DataFilesPanel.FontWeight = 'bold';
            app.DataFilesPanel.FontSize = 12;

            dataFilesGrid = uigridlayout(app.DataFilesPanel);
            dataFilesGrid.ColumnWidth = {175, 80, '1x', 50, 115};
            % rows: 1 source toggle, 2 file-1, 3 file-2, 4 count-1, 5 count-2,
            %       6 spectral grid (expands), 7 caption.
            % Rows 4-5 are shown only for the circular phase-group analyses
            % (applyAnalysisLayout); row 6 grows in spectral mode.
            dataFilesGrid.RowHeight = {'fit', 'fit', 'fit', 0, 0, 0, 'fit'};
            dataFilesGrid.RowSpacing = 10;
            dataFilesGrid.ColumnSpacing = 8;
            dataFilesGrid.Padding = [12 12 12 12];

            % Row 1: data-source toggle (Files | Spectral folder)
            app.DataSourceLabel = uilabel(dataFilesGrid);
            app.DataSourceLabel.Layout.Row = 1; app.DataSourceLabel.Layout.Column = 1;
            app.DataSourceLabel.Text = 'Data source:'; app.DataSourceLabel.FontWeight = 'bold';

            app.DataSourceDropDown = uidropdown(dataFilesGrid);
            app.DataSourceDropDown.Layout.Row = 1; app.DataSourceDropDown.Layout.Column = [2 3];
            app.DataSourceDropDown.Items = {'Files', 'Spectral folder'};
            app.DataSourceDropDown.Value = 'Files';
            app.DataSourceDropDown.ValueChangedFcn = createCallbackFcn(app, @DataSourceDropDownValueChanged, true);

            % --- File-source rows (rows 2-3): shown in 'Files' mode ---
            % File 1 row (title relabelled per analysis in applyAnalysisLayout)
            app.Data1FileTitleLabel = uilabel(dataFilesGrid);
            app.Data1FileTitleLabel.Layout.Row = 2; app.Data1FileTitleLabel.Layout.Column = 1;
            app.Data1FileTitleLabel.Text = 'Data File:';
            app.Data1FileTitleLabel.WordWrap = 'on';
            app.Data1FileTitleLabel.VerticalAlignment = 'center';

            app.Data1FileButton = uibutton(dataFilesGrid);
            app.Data1FileButton.Layout.Row = 2; app.Data1FileButton.Layout.Column = 2;
            app.Data1FileButton.Text = 'Browse...';
            app.Data1FileButton.ButtonPushedFcn = createCallbackFcn(app, @Data1FileButtonPushed, true);

            app.Data1FileLabel = uilabel(dataFilesGrid);
            app.Data1FileLabel.Layout.Row = 2; app.Data1FileLabel.Layout.Column = 3;
            app.Data1FileLabel.Text = 'No file selected'; app.Data1FileLabel.FontAngle = 'italic';

            data1SheetLabel = uilabel(dataFilesGrid);
            data1SheetLabel.Layout.Row = 2; data1SheetLabel.Layout.Column = 4;
            data1SheetLabel.Text = 'Sheet:';

            app.Data1SheetDropDown = uidropdown(dataFilesGrid);
            app.Data1SheetDropDown.Layout.Row = 2; app.Data1SheetDropDown.Layout.Column = 5;
            app.Data1SheetDropDown.ValueChangedFcn = createCallbackFcn(app, @Data1SheetChanged, true);

            % File 2 row (hidden for single-file GLM presets)
            app.Data2FileTitleLabel = uilabel(dataFilesGrid);
            app.Data2FileTitleLabel.Layout.Row = 3; app.Data2FileTitleLabel.Layout.Column = 1;
            app.Data2FileTitleLabel.Text = 'Data 2 File:';
            app.Data2FileTitleLabel.WordWrap = 'on';
            app.Data2FileTitleLabel.VerticalAlignment = 'center';

            app.Data2FileButton = uibutton(dataFilesGrid);
            app.Data2FileButton.Layout.Row = 3; app.Data2FileButton.Layout.Column = 2;
            app.Data2FileButton.Text = 'Browse...';
            app.Data2FileButton.ButtonPushedFcn = createCallbackFcn(app, @Data2FileButtonPushed, true);

            app.Data2FileLabel = uilabel(dataFilesGrid);
            app.Data2FileLabel.Layout.Row = 3; app.Data2FileLabel.Layout.Column = 3;
            app.Data2FileLabel.Text = 'No file selected'; app.Data2FileLabel.FontAngle = 'italic';

            app.Data2SheetTitleLabel = uilabel(dataFilesGrid);
            app.Data2SheetTitleLabel.Layout.Row = 3; app.Data2SheetTitleLabel.Layout.Column = 4;
            app.Data2SheetTitleLabel.Text = 'Sheet:';

            app.Data2SheetDropDown = uidropdown(dataFilesGrid);
            app.Data2SheetDropDown.Layout.Row = 3; app.Data2SheetDropDown.Layout.Column = 5;
            app.Data2SheetDropDown.ValueChangedFcn = createCallbackFcn(app, @Data2SheetChanged, true);

            % --- Event-count rows (4-5): circular phase-group analyses only ---
            % One wide subjects x channels file per group, same subject column
            % and same channel names as the matching angle file.
            countTip = ['Number of detected events (e.g. spindles) behind each subject''s ' ...
                        'preferred phase, one column per channel. Required: log(count) is ' ...
                        'entered as a covariate, because a subject with few events has a ' ...
                        'noisier and genuinely wider-looking phase estimate, which the test ' ...
                        'would otherwise read as a real group difference.'];

            app.Count1FileTitleLabel = uilabel(dataFilesGrid);
            app.Count1FileTitleLabel.Layout.Row = 4; app.Count1FileTitleLabel.Layout.Column = 1;
            app.Count1FileTitleLabel.Text = 'Group A Event Count File:';
            app.Count1FileTitleLabel.WordWrap = 'on';
            app.Count1FileTitleLabel.VerticalAlignment = 'center';
            app.Count1FileTitleLabel.Tooltip = countTip;
            app.Count1FileTitleLabel.Visible = 'off';

            app.Count1FileButton = uibutton(dataFilesGrid);
            app.Count1FileButton.Layout.Row = 4; app.Count1FileButton.Layout.Column = 2;
            app.Count1FileButton.Text = 'Browse...';
            app.Count1FileButton.Tooltip = countTip;
            app.Count1FileButton.Visible = 'off';
            app.Count1FileButton.ButtonPushedFcn = createCallbackFcn(app, @Count1FileButtonPushed, true);

            app.Count1FileLabel = uilabel(dataFilesGrid);
            app.Count1FileLabel.Layout.Row = 4; app.Count1FileLabel.Layout.Column = 3;
            app.Count1FileLabel.Text = 'No file selected'; app.Count1FileLabel.FontAngle = 'italic';
            app.Count1FileLabel.Visible = 'off';

            app.Count1SheetTitleLabel = uilabel(dataFilesGrid);
            app.Count1SheetTitleLabel.Layout.Row = 4; app.Count1SheetTitleLabel.Layout.Column = 4;
            app.Count1SheetTitleLabel.Text = 'Sheet:';
            app.Count1SheetTitleLabel.Visible = 'off';

            app.Count1SheetDropDown = uidropdown(dataFilesGrid);
            app.Count1SheetDropDown.Layout.Row = 4; app.Count1SheetDropDown.Layout.Column = 5;
            app.Count1SheetDropDown.Visible = 'off';
            app.Count1SheetDropDown.ValueChangedFcn = createCallbackFcn(app, @Count1SheetChanged, true);

            app.Count2FileTitleLabel = uilabel(dataFilesGrid);
            app.Count2FileTitleLabel.Layout.Row = 5; app.Count2FileTitleLabel.Layout.Column = 1;
            app.Count2FileTitleLabel.Text = 'Group B Event Count File:';
            app.Count2FileTitleLabel.WordWrap = 'on';
            app.Count2FileTitleLabel.VerticalAlignment = 'center';
            app.Count2FileTitleLabel.Tooltip = countTip;
            app.Count2FileTitleLabel.Visible = 'off';

            app.Count2FileButton = uibutton(dataFilesGrid);
            app.Count2FileButton.Layout.Row = 5; app.Count2FileButton.Layout.Column = 2;
            app.Count2FileButton.Text = 'Browse...';
            app.Count2FileButton.Tooltip = countTip;
            app.Count2FileButton.Visible = 'off';
            app.Count2FileButton.ButtonPushedFcn = createCallbackFcn(app, @Count2FileButtonPushed, true);

            app.Count2FileLabel = uilabel(dataFilesGrid);
            app.Count2FileLabel.Layout.Row = 5; app.Count2FileLabel.Layout.Column = 3;
            app.Count2FileLabel.Text = 'No file selected'; app.Count2FileLabel.FontAngle = 'italic';
            app.Count2FileLabel.Visible = 'off';

            app.Count2SheetTitleLabel = uilabel(dataFilesGrid);
            app.Count2SheetTitleLabel.Layout.Row = 5; app.Count2SheetTitleLabel.Layout.Column = 4;
            app.Count2SheetTitleLabel.Text = 'Sheet:';
            app.Count2SheetTitleLabel.Visible = 'off';

            app.Count2SheetDropDown = uidropdown(dataFilesGrid);
            app.Count2SheetDropDown.Layout.Row = 5; app.Count2SheetDropDown.Layout.Column = 5;
            app.Count2SheetDropDown.Visible = 'off';
            app.Count2SheetDropDown.ValueChangedFcn = createCallbackFcn(app, @Count2SheetChanged, true);

            % --- Spectral-folder controls (row 6 nested grid): shown in 'Spectral folder' mode ---
            app.SpecGrid = uigridlayout(dataFilesGrid);
            app.SpecGrid.Layout.Row = 6; app.SpecGrid.Layout.Column = [1 5];
            app.SpecGrid.ColumnWidth = {'1x', '1x', '1x'};
            % rows: 1 buttons, 2 hint, 3 folder table (expands),
            %       4 band/stage/type headers, 5 band/stage/type listboxes (expand).
            app.SpecGrid.RowHeight = {'fit', 'fit', '1x', 'fit', '1x'};
            app.SpecGrid.Padding = [0 0 0 0]; app.SpecGrid.RowSpacing = 8; app.SpecGrid.ColumnSpacing = 12;
            app.SpecGrid.Visible = 'off';

            % Row 1: add / remove folder buttons
            app.SpecAddFolderButton = uibutton(app.SpecGrid);
            app.SpecAddFolderButton.Layout.Row = 1; app.SpecAddFolderButton.Layout.Column = 1;
            app.SpecAddFolderButton.Text = 'Add folder...';
            app.SpecAddFolderButton.ButtonPushedFcn = createCallbackFcn(app, @SpecAddFolderButtonPushed, true);

            app.SpecRemoveFolderButton = uibutton(app.SpecGrid);
            app.SpecRemoveFolderButton.Layout.Row = 1; app.SpecRemoveFolderButton.Layout.Column = 2;
            app.SpecRemoveFolderButton.Text = 'Remove selected';
            app.SpecRemoveFolderButton.ButtonPushedFcn = createCallbackFcn(app, @SpecRemoveFolderButtonPushed, true);

            % Row 2: hint
            specHint = uilabel(app.SpecGrid);
            specHint.Layout.Row = 2; specHint.Layout.Column = [1 3];
            specHint.Text = 'Each folder = one group / condition (a level of the design factor).';
            specHint.FontAngle = 'italic'; specHint.FontColor = [0.4 0.4 0.4];

            % Row 3: folders table (path + editable label) - grows with the window (3-5 rows)
            app.SpecFolderTable = uitable(app.SpecGrid);
            app.SpecFolderTable.Layout.Row = 3; app.SpecFolderTable.Layout.Column = [1 3];
            app.SpecFolderTable.ColumnName = {'Folder', 'Label'};
            app.SpecFolderTable.ColumnEditable = [false true];
            app.SpecFolderTable.ColumnWidth = {'1x', 120};
            app.SpecFolderTable.SelectionType = 'row';
            app.SpecFolderTable.Data = cell(0,2);
            app.SpecFolderTable.CellEditCallback = createCallbackFcn(app, @SpecSelectionChanged, true);

            % Row 4: band / stage / power-type headers
            app.SpecBandTitleLabel = uilabel(app.SpecGrid);
            app.SpecBandTitleLabel.Layout.Row = 4; app.SpecBandTitleLabel.Layout.Column = 1;
            app.SpecBandTitleLabel.Text = 'Bands (multi-select):';

            app.SpecStageTitleLabel = uilabel(app.SpecGrid);
            app.SpecStageTitleLabel.Layout.Row = 4; app.SpecStageTitleLabel.Layout.Column = 2;
            app.SpecStageTitleLabel.Text = 'Stages (multi-select):';

            app.SpecTypeTitleLabel = uilabel(app.SpecGrid);
            app.SpecTypeTitleLabel.Layout.Row = 4; app.SpecTypeTitleLabel.Layout.Column = 3;
            app.SpecTypeTitleLabel.Text = 'Power type (multi-select):';

            % Row 5: band / stage / power-type multi-select listboxes - tall enough to show all items
            app.SpecBandListBox = uilistbox(app.SpecGrid);
            app.SpecBandListBox.Layout.Row = 5; app.SpecBandListBox.Layout.Column = 1;
            app.SpecBandListBox.Multiselect = 'on';
            app.SpecBandListBox.Items = {}; app.SpecBandListBox.Value = {};
            app.SpecBandListBox.ValueChangedFcn = createCallbackFcn(app, @SpecSelectionChanged, true);

            app.SpecStageListBox = uilistbox(app.SpecGrid);
            app.SpecStageListBox.Layout.Row = 5; app.SpecStageListBox.Layout.Column = 2;
            app.SpecStageListBox.Multiselect = 'on';
            app.SpecStageListBox.Items = {}; app.SpecStageListBox.Value = {};
            app.SpecStageListBox.ValueChangedFcn = createCallbackFcn(app, @SpecSelectionChanged, true);

            app.SpecTypeListBox = uilistbox(app.SpecGrid);
            app.SpecTypeListBox.Layout.Row = 5; app.SpecTypeListBox.Layout.Column = 3;
            app.SpecTypeListBox.Multiselect = 'on';
            app.SpecTypeListBox.Items = {'absolute', 'normalized'};
            app.SpecTypeListBox.Value = {'absolute'};
            app.SpecTypeListBox.ValueChangedFcn = createCallbackFcn(app, @SpecSelectionChanged, true);

            % Detected channels / metadata caption (auto-filled on file/folder load)
            app.DetectCaptionLabel = uilabel(dataFilesGrid);
            app.DetectCaptionLabel.Layout.Row = 7;
            app.DetectCaptionLabel.Layout.Column = [1 5];
            app.DetectCaptionLabel.Text = '';
            app.DetectCaptionLabel.FontAngle = 'italic';
            app.DetectCaptionLabel.FontColor = [0.2 0.4 0.8];
            app.DetectCaptionLabel.FontSize = 11;

            % ============ 3. Analysis parameters ============
            app.ParametersPanel = uipanel(app.GridLayout);
            app.ParametersPanel.Layout.Row = 6;
            app.ParametersPanel.Layout.Column = 1;
            app.ParametersPanel.Title = '3.  Analysis parameters';
            app.ParametersPanel.FontWeight = 'bold';
            app.ParametersPanel.FontSize = 12;

            paramGrid = uigridlayout(app.ParametersPanel);
            paramGrid.ColumnWidth = {140, '1x', 110, '1x'};   % wider 1st col so "Recording system:" isn't truncated
            paramGrid.RowHeight = {30, 30};
            paramGrid.RowSpacing = 8;
            paramGrid.ColumnSpacing = 12;
            paramGrid.Padding = [10 10 10 10];

            channelsLabel = uilabel(paramGrid);
            channelsLabel.Layout.Row = 1;
            channelsLabel.Layout.Column = 1;
            channelsLabel.Text = 'Recording system:';

            app.ChannelsDropDown = uidropdown(paramGrid);
            app.ChannelsDropDown.Layout.Row = 1;
            app.ChannelsDropDown.Layout.Column = 2;
            % Plain-language system names (Items) -> internal keys (ItemsData),
            % from the single source of truth snpm_montage_registry.
            reg = snpm_montage_registry();
            app.ChannelsDropDown.Items = {reg.display};
            app.ChannelsDropDown.ItemsData = {reg.key};
            app.ChannelsDropDown.ValueChangedFcn = createCallbackFcn(app, @ChannelsChanged, true);

            dataTypeLabel = uilabel(paramGrid);
            dataTypeLabel.Layout.Row = 1;
            dataTypeLabel.Layout.Column = 3;
            dataTypeLabel.Text = 'Data Type:';

            app.DataTypeDropDown = uidropdown(paramGrid);
            app.DataTypeDropDown.Layout.Row = 1;
            app.DataTypeDropDown.Layout.Column = 4;
            app.DataTypeDropDown.Items = {'absolute', 'logscale', 'normalize'};

            tailsLabel = uilabel(paramGrid);
            tailsLabel.Layout.Row = 2;
            tailsLabel.Layout.Column = 1;
            tailsLabel.Text = 'Tails:';

            app.TailsDropDown = uidropdown(paramGrid);
            app.TailsDropDown.Layout.Row = 2;
            app.TailsDropDown.Layout.Column = 2;
            app.TailsDropDown.Items = {'both', 'left', 'right'};

            permutationLabel = uilabel(paramGrid);
            permutationLabel.Layout.Row = 2;
            permutationLabel.Layout.Column = 3;
            permutationLabel.Text = 'Permutations:';

            app.PermutationField = uieditfield(paramGrid, 'numeric');
            app.PermutationField.Layout.Row = 2;
            app.PermutationField.Layout.Column = 4;
            app.PermutationField.ValueDisplayFormat = '%.0f';

            % ============ Column roles (GLM presets only) ============
            app.RolesPanel = uipanel(app.GridLayout);
            app.RolesPanel.Layout.Row = app.RolesRow;
            app.RolesPanel.Layout.Column = 1;
            app.RolesPanel.Title = 'Column roles (pick which columns define the design)';
            app.RolesPanel.FontWeight = 'bold';
            app.RolesPanel.FontSize = 12;

            rolesGrid = uigridlayout(app.RolesPanel);
            rolesGrid.ColumnWidth = {120, '1x', 120, '1x'};
            rolesGrid.RowHeight = {30, 30, 64};   % taller last row so the Covariate-cols listbox shows several items
            rolesGrid.RowSpacing = 8;
            rolesGrid.ColumnSpacing = 12;
            rolesGrid.Padding = [10 10 10 10];

            lbMeta = uilabel(rolesGrid); lbMeta.Layout.Row=1; lbMeta.Layout.Column=1; lbMeta.Text='Meta cols:';
            app.MetaColsField = uieditfield(rolesGrid,'text'); app.MetaColsField.Layout.Row=1; app.MetaColsField.Layout.Column=2;
            app.MetaColsField.Placeholder = 'auto-detected on file load (editable)';
            lbGroup = uilabel(rolesGrid); lbGroup.Layout.Row=1; lbGroup.Layout.Column=3; lbGroup.Text='Group col:';
            app.GroupColField = uidropdown(rolesGrid); app.GroupColField.Layout.Row=1; app.GroupColField.Layout.Column=4;
            lbPred = uilabel(rolesGrid); lbPred.Layout.Row=2; lbPred.Layout.Column=1; lbPred.Text='Predictor:';
            app.PredictorColField = uidropdown(rolesGrid); app.PredictorColField.Layout.Row=2; app.PredictorColField.Layout.Column=2;
            lbCond = uilabel(rolesGrid); lbCond.Layout.Row=2; lbCond.Layout.Column=3; lbCond.Text='Condition col:';
            app.ConditionColField = uidropdown(rolesGrid); app.ConditionColField.Layout.Row=2; app.ConditionColField.Layout.Column=4;
            lbSubj = uilabel(rolesGrid); lbSubj.Layout.Row=3; lbSubj.Layout.Column=1; lbSubj.Text='Subject col:';
            app.SubjectColField = uidropdown(rolesGrid); app.SubjectColField.Layout.Row=3; app.SubjectColField.Layout.Column=2;
            lbCov = uilabel(rolesGrid); lbCov.Layout.Row=3; lbCov.Layout.Column=3; lbCov.Text='Covariate cols:';
            app.CovColsField = uilistbox(rolesGrid); app.CovColsField.Layout.Row=3; app.CovColsField.Layout.Column=4;
            app.CovColsField.Multiselect = 'on'; app.CovColsField.Items = {}; app.CovColsField.Value = {};
            % role pickers fire the readiness re-check when changed
            for dd = {app.GroupColField, app.PredictorColField, app.ConditionColField, app.SubjectColField}
                dd{1}.ValueChangedFcn = createCallbackFcn(app, @ChannelsChanged, true);
            end

            % role key -> {label, field}; updateRoleFields toggles Visible per analysis
            app.RoleRows = { ...
                'meta',      lbMeta,  app.MetaColsField; ...
                'group',     lbGroup, app.GroupColField; ...
                'predictor', lbPred,  app.PredictorColField; ...
                'condition', lbCond,  app.ConditionColField; ...
                'subject',   lbSubj,  app.SubjectColField; ...
                'cov',       lbCov,   app.CovColsField};

            % ============ Phase / angle settings (circular analyses only) ============
            % Units and convention are REQUIRED with no preselected value. A
            % preselected default would be a guess about how the user's detector
            % exported its angles, and a wrong guess rotates every angle in the
            % dataset without anything in the output looking wrong.
            app.CircPanel = uipanel(app.GridLayout);
            app.CircPanel.Layout.Row = app.CircRow;
            app.CircPanel.Layout.Column = 1;
            app.CircPanel.Title = 'Phase / angle settings (required)';
            app.CircPanel.FontWeight = 'bold';
            app.CircPanel.FontSize = 12;
            app.CircPanel.Visible = 'off';

            circGrid = uigridlayout(app.CircPanel);
            circGrid.ColumnWidth = {150, '1x', 150, '1x'};
            circGrid.RowHeight = {30, 30, 'fit'};
            circGrid.RowSpacing = 8;
            circGrid.ColumnSpacing = 12;
            circGrid.Padding = [10 10 10 10];

            app.CircUnitsLabel = uilabel(circGrid);
            app.CircUnitsLabel.Layout.Row = 1; app.CircUnitsLabel.Layout.Column = 1;
            app.CircUnitsLabel.Text = 'Angle units:';

            app.CircUnitsDropDown = uidropdown(circGrid);
            app.CircUnitsDropDown.Layout.Row = 1; app.CircUnitsDropDown.Layout.Column = 2;
            app.CircUnitsDropDown.Items     = {'(choose)', 'radians (-pi to pi or 0 to 2pi)', 'degrees (0 to 360)'};
            app.CircUnitsDropDown.ItemsData = {'', 'radians', 'degrees'};
            app.CircUnitsDropDown.Value = '';
            app.CircUnitsDropDown.Tooltip = ['The units your detector exported. There is no default: ' ...
                'reading degrees as radians silently scrambles every angle.'];
            app.CircUnitsDropDown.ValueChangedFcn = createCallbackFcn(app, @ChannelsChanged, true);

            app.CircConventionLabel = uilabel(circGrid);
            app.CircConventionLabel.Layout.Row = 1; app.CircConventionLabel.Layout.Column = 3;
            app.CircConventionLabel.Text = 'Where zero sits:';

            app.CircConventionDropDown = uidropdown(circGrid);
            app.CircConventionDropDown.Layout.Row = 1; app.CircConventionDropDown.Layout.Column = 4;
            app.CircConventionDropDown.Items = { ...
                '(choose)', ...
                'Zero = slow-oscillation up-state peak (standard)', ...
                'Luna COUPL_ANGLE (zero = SO zero crossing)', ...
                'YASA (zero = SO up-state peak)', ...
                'TurtleWave before v4.0 (180 degrees off)', ...
                'Custom offset (enter below)'};
            app.CircConventionDropDown.ItemsData = { ...
                '', 'literature_uppeak0', 'luna_zerocross0', 'yasa_uppeak0', ...
                'turtlewave_pre_v4', 'custom'};
            app.CircConventionDropDown.Value = '';
            app.CircConventionDropDown.Tooltip = ['Which point of the slow oscillation your ' ...
                'exported angle calls zero. This is a live trap, not a hypothetical: Luna''s ' ...
                'COUPL_ANGLE puts the up-state peak at 270 degrees, and any TurtleWave export ' ...
                'from before v4.0 is exactly 180 degrees wrong. Choosing the wrong one rotates ' ...
                'every angle and moves the reported preferred phase.'];
            app.CircConventionDropDown.ValueChangedFcn = createCallbackFcn(app, @CircConventionChanged, true);

            app.CircOffsetLabel = uilabel(circGrid);
            app.CircOffsetLabel.Layout.Row = 2; app.CircOffsetLabel.Layout.Column = 1;
            app.CircOffsetLabel.Text = 'Custom offset (deg):';
            app.CircOffsetLabel.Enable = 'off';

            app.CircOffsetField = uieditfield(circGrid, 'numeric');
            app.CircOffsetField.Layout.Row = 2; app.CircOffsetField.Layout.Column = 2;
            app.CircOffsetField.Value = 0;
            app.CircOffsetField.Enable = 'off';
            app.CircOffsetField.Tooltip = ['Degrees to add so that zero lands on the ' ...
                'slow-oscillation up-state peak. Only used for the Custom convention.'];

            app.CircNoteLabel = uilabel(circGrid);
            app.CircNoteLabel.Layout.Row = 3; app.CircNoteLabel.Layout.Column = [1 4];
            app.CircNoteLabel.WordWrap = 'on';
            app.CircNoteLabel.FontSize = 11;
            app.CircNoteLabel.FontColor = [0.55 0.30 0.05];
            app.CircNoteLabel.Text = ['Angles are converted to the standard convention (zero = ' ...
                'up-state peak) before testing, so this setting changes the result. Data Type and ' ...
                'Tails are fixed for circular analyses: angles cannot be log-scaled or z-scored, ' ...
                'and the statistic is inherently two-sided.'];

            % ============ Covariate file (correlation analyses only) ============
            app.CovariatePanel = uipanel(app.GridLayout);
            app.CovariatePanel.Layout.Row = app.CovRow;
            app.CovariatePanel.Layout.Column = 1;
            app.CovariatePanel.Title = 'Covariates (optional, for partial correlation)';
            app.CovariatePanel.FontWeight = 'bold';
            app.CovariatePanel.FontSize = 12;

            covariateGrid = uigridlayout(app.CovariatePanel);
            covariateGrid.ColumnWidth = {'1x', 110, '1x'};
            covariateGrid.RowHeight = {28};

            app.UseCovariatesCheckBox = uicheckbox(covariateGrid);
            app.UseCovariatesCheckBox.Layout.Row = 1;
            app.UseCovariatesCheckBox.Layout.Column = 1;
            app.UseCovariatesCheckBox.Text = 'Include covariate file';

            app.CovariateFileButton = uibutton(covariateGrid);
            app.CovariateFileButton.Layout.Row = 1;
            app.CovariateFileButton.Layout.Column = 2;
            app.CovariateFileButton.Text = 'Select CSV...';
            app.CovariateFileButton.ButtonPushedFcn = createCallbackFcn(app, @CovariateFileButtonPushed, true);

            app.CovariateFileLabel = uilabel(covariateGrid);
            app.CovariateFileLabel.Layout.Row = 1;
            app.CovariateFileLabel.Layout.Column = 3;
            app.CovariateFileLabel.Text = 'No covariate file selected';
            app.CovariateFileLabel.FontAngle = 'italic';
            app.CovariateFileLabel.FontSize = 11;
            app.CovariateFileLabel.WordWrap = 'on';
            app.CovariateFileLabel.VerticalAlignment = 'center';

            % ============ 4. Output folder ============
            app.OutputPanel = uipanel(app.GridLayout);
            app.OutputPanel.Layout.Row = 7;
            app.OutputPanel.Layout.Column = 1;
            app.OutputPanel.Title = '4.  Output folder';
            app.OutputPanel.FontWeight = 'bold';
            app.OutputPanel.FontSize = 12;

            outputGrid = uigridlayout(app.OutputPanel);
            outputGrid.ColumnWidth = {'1x'};
            outputGrid.RowHeight = {32, 'fit'};
            outputGrid.RowSpacing = 6;
            outputGrid.Padding = [10 10 10 10];

            app.OutputPathButton = uibutton(outputGrid);
            app.OutputPathButton.Layout.Row = 1;
            app.OutputPathButton.Layout.Column = 1;
            app.OutputPathButton.Text = 'Select Directory...';
            app.OutputPathButton.ButtonPushedFcn = createCallbackFcn(app, @OutputPathButtonPushed, true);

            app.OutputPathLabel = uilabel(outputGrid);
            app.OutputPathLabel.Layout.Row = 2;
            app.OutputPathLabel.Layout.Column = 1;
            app.OutputPathLabel.Text = 'No path selected';
            app.OutputPathLabel.FontAngle = 'italic';
            app.OutputPathLabel.FontSize = 11;
            app.OutputPathLabel.WordWrap = 'on';
            app.OutputPathLabel.VerticalAlignment = 'top';

            % ============ Before you can run: (validation checklist) ============
            app.ValidationPanel = uipanel(app.GridLayout);
            app.ValidationPanel.Layout.Row = 11;
            app.ValidationPanel.Layout.Column = 1;
            app.ValidationPanel.Title = 'Before you can run:';
            app.ValidationPanel.FontWeight = 'bold';
            app.ValidationPanel.FontSize = 12;
            app.ValidationPanel.ForegroundColor = [0.80 0.20 0.20];

            validationGrid = uigridlayout(app.ValidationPanel);
            validationGrid.ColumnWidth = {'1x'};
            validationGrid.RowHeight = {'fit'};
            validationGrid.Padding = [10 8 10 8];

            app.ValidationLabel = uilabel(validationGrid);
            app.ValidationLabel.Layout.Row = 1;
            app.ValidationLabel.Layout.Column = 1;
            app.ValidationLabel.WordWrap = 'on';
            app.ValidationLabel.VerticalAlignment = 'top';
            app.ValidationLabel.FontColor = [0.80 0.20 0.20];
            app.ValidationLabel.Text = '';

            % ============ Run / Reset buttons ============
            buttonGrid = uigridlayout(app.GridLayout);
            buttonGrid.Layout.Row = 12;
            buttonGrid.Layout.Column = 1;
            buttonGrid.ColumnWidth = {'1x', '1x'};
            buttonGrid.RowHeight = {42};
            buttonGrid.ColumnSpacing = 12;

            app.RunAnalysisButton = uibutton(buttonGrid);
            app.RunAnalysisButton.Layout.Row = 1;
            app.RunAnalysisButton.Layout.Column = 1;
            app.RunAnalysisButton.Text = 'Run Analysis';
            app.RunAnalysisButton.FontSize = 14;
            app.RunAnalysisButton.FontWeight = 'bold';
            app.RunAnalysisButton.BackgroundColor = [0.2 0.6 0.2];
            app.RunAnalysisButton.FontColor = [1 1 1];
            app.RunAnalysisButton.ButtonPushedFcn = createCallbackFcn(app, @RunAnalysisButtonPushed, true);

            app.ResetButton = uibutton(buttonGrid);
            app.ResetButton.Layout.Row = 1;
            app.ResetButton.Layout.Column = 2;
            app.ResetButton.Text = 'Reset';
            app.ResetButton.FontSize = 14;
            app.ResetButton.FontWeight = 'bold';
            app.ResetButton.BackgroundColor = [0.8 0.4 0.4];
            app.ResetButton.FontColor = [1 1 1];
            app.ResetButton.ButtonPushedFcn = createCallbackFcn(app, @ResetButtonPushed, true);

            % Status Label
            app.StatusLabel = uilabel(app.GridLayout);
            app.StatusLabel.Layout.Row = 13;
            app.StatusLabel.Layout.Column = 1;
            app.StatusLabel.HorizontalAlignment = 'center';
            app.StatusLabel.FontAngle = 'italic';

            % Results Panel
            app.ResultsPanel = uipanel(app.GridLayout);
            app.ResultsPanel.Layout.Row = 14;
            app.ResultsPanel.Layout.Column = 1;
            app.ResultsPanel.Title = 'Results';
            app.ResultsPanel.FontWeight = 'bold';

            resultsGrid = uigridlayout(app.ResultsPanel);
            resultsGrid.ColumnWidth = {'1x'};
            resultsGrid.RowHeight = {'1x'};

            app.ResultsTextArea = uitextarea(resultsGrid);
            app.ResultsTextArea.Layout.Row = 1;
            app.ResultsTextArea.Layout.Column = 1;
            app.ResultsTextArea.Editable = 'off';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % Scripted input setting (used by the headless test harness).
    % These do exactly what the corresponding Browse... buttons do once the
    % file dialog has returned: record the path, list the sheets, re-detect the
    % channel columns and re-run the readiness gate. They exist because the
    % readiness gate is a correctness guard and there is no display in CI, so
    % the only way to prove it blocks the invalid configurations and allows the
    % valid one is to drive it without a file dialog. No UI code calls them.
    methods (Access = public, Hidden)

        % SLOT is 'data1' | 'data2' | 'count1' | 'count2' | 'output'.
        % For 'output', FILE is the folder path and SHEET is ignored.
        function setInputPath(app, slot, file, sheet)
            if nargin < 4, sheet = ''; end
            if strcmp(slot, 'output')
                app.output_path = file;
                app.OutputPathLabel.Text = file;
                checkReadyToRun(app);
                return;
            end
            switch slot
                case 'data1',  dd = app.Data1SheetDropDown;  lb = app.Data1FileLabel;
                case 'data2',  dd = app.Data2SheetDropDown;  lb = app.Data2FileLabel;
                case 'count1', dd = app.Count1SheetDropDown; lb = app.Count1FileLabel;
                case 'count2', dd = app.Count2SheetDropDown; lb = app.Count2FileLabel;
                otherwise, error('SnPMAnalysisGui:badSlot', 'Unknown input slot "%s".', slot);
            end
            if isempty(file)   % clearing an input
                dd.Items = {}; lb.Text = 'No file selected'; n = 0;
                switch slot
                    case 'data1', initRolePickers(app);   % also drops the primary detection
                    case 'data2', clearData2Detection(app);
                end
            else
                if isempty(sheet)
                    if endsWith(file, {'.xlsx','.xls'})
                        sn = sheetnames(file); sheet = char(sn(1));
                    else
                        sheet = 'CSV File';
                    end
                end
                dd.Items = {sheet}; dd.Value = sheet;
                [~, b, e] = fileparts(file); lb.Text = [b e];
                [~, n, ~, B] = detectAndFill(app, file, sheet, strcmp(slot,'data1'));
                if strcmp(slot, 'data2'), setData2Detection(app, n, B); end
            end
            switch slot
                case 'data1',  app.data1_file  = file; app.NChanAngle1 = n;
                case 'data2',  app.data2_file  = file; app.NChanAngle2 = n;
                case 'count1', app.count1_file = file; app.NChanCount1 = n;
                case 'count2', app.count2_file = file; app.NChanCount2 = n;
            end
            checkReadyToRun(app);
        end

        % Select a comparison exactly as the dropdown callback would.
        function setComparison(app, key)
            app.ComparisonDropDown.Value = key;
            ComparisonDropDownValueChanged(app, []);
        end

        % Read-only snapshot of the state the readiness gate acts on.
        function s = uiState(app)
            s = struct( ...
                'comparison',      app.ComparisonDropDown.Value, ...
                'datatypeItems',   {app.DataTypeDropDown.Items}, ...
                'datatypeValue',   app.DataTypeDropDown.Value, ...
                'datatypeEnable',  char(app.DataTypeDropDown.Enable), ...
                'tailItems',       {app.TailsDropDown.Items}, ...
                'tailValue',       app.TailsDropDown.Value, ...
                'tailEnable',      char(app.TailsDropDown.Enable), ...
                'systemKeys',      {app.ChannelsDropDown.ItemsData}, ...
                'systemValue',     app.ChannelsDropDown.Value, ...
                'comparisonKeys',  {app.ComparisonDropDown.ItemsData}, ...
                'label1',          app.Data1FileTitleLabel.Text, ...
                'label2',          app.Data2FileTitleLabel.Text, ...
                'countLabel1',     app.Count1FileTitleLabel.Text, ...
                'countLabel2',     app.Count2FileTitleLabel.Text, ...
                'countVisible',    char(app.Count1FileButton.Visible), ...
                'circPanelVisible',char(app.CircPanel.Visible), ...
                'circOffsetEnable',char(app.CircOffsetField.Enable), ...
                'hint',            app.HintLabel.Text, ...
                'detectCaption',   {cellstr(app.DetectCaptionLabel.Text)}, ...
                'metaCols',        app.MetaColsField.Value, ...
                'chanNames',       {app.DetectedChanNames}, ...
                'nChan1',          app.DetectedNChan, ...
                'nChan2',          app.DetectedNChan2, ...
                'kind1',           normalizeDetect(app.DetectB1).kind, ...
                'kind2',           normalizeDetect(app.DetectB2).kind, ...
                'runEnable',       char(app.RunAnalysisButton.Enable), ...
                'validationTitle', app.ValidationPanel.Title, ...
                'validation',      {cellstr(app.ValidationLabel.Text)});
        end

        % Set a circular setting ('units' | 'convention' | 'offset' | 'system')
        % through the same callback path the widget uses.
        function setCircOption(app, what, value)
            switch what
                case 'units',      app.CircUnitsDropDown.Value = value; ChannelsChanged(app, []);
                case 'convention', app.CircConventionDropDown.Value = value; CircConventionChanged(app, []);
                case 'offset',     app.CircOffsetField.Value = value; checkReadyToRun(app);
                case 'system',     app.ChannelsDropDown.Value = value; ChannelsChanged(app, []);
                otherwise, error('SnPMAnalysisGui:badOption', 'Unknown option "%s".', what);
            end
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SnPMAnalysisGui

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end

% ---- local helpers (class file local functions) ----
function c = splitCols(s)
    % split a comma/space-separated column list into a cellstr, dropping blanks
    if isempty(s), c = {}; return; end
    parts = strtrim(strsplit(s, {',',';'}));
    c = parts(~cellfun(@isempty, parts));
end

function info = scan_spectral_folder(folder)
    % Discover bands (from EEG.features/EEG.bands of the first file), stages
    % (from filenames) and channel labels in one flat spectral-power folder
    % (a cohort/condition folder of sub-*_..._powerspect.mat files).
    files = dir(fullfile(folder, '*powerspect*.mat'));   % PSG _powerspect.mat + KDT _powerspect_final.mat
    assert(~isempty(files), 'No *powerspect*.mat files found in the selected folder.');
    stg = regexp({files.name}, 'desc-([^_.]+)', 'tokens', 'once');
    stg = unique([stg{:}], 'stable');
    order  = {'n1','n2','n3','rem','eyesopen','eyesclosed'};
    stages = [order(ismember(order, stg)), setdiff(stg, order, 'stable')];
    S = load(fullfile(files(1).folder, files(1).name));
    B = [];
    if isfield(S.EEG,'bands') && ~isempty(S.EEG.bands),         B = S.EEG.bands;
    elseif isfield(S.EEG,'features') && ~isempty(S.EEG.features), B = S.EEG.features; end
    assert(~isempty(B), 'First file has no EEG.bands / EEG.features band-power struct.');
    info = struct('bands', {unique({B.label}, 'stable')}, 'stages', {stages}, ...
                  'chanlabels', {{S.EEG.chanlocs.labels}});
end

function p = specMinP(rs)
    p = NaN;
    if isfield(rs,'Clusters') && ~isempty(rs.Clusters) && isfield(rs.Clusters,'p')
        ps = [rs.Clusters.p]; if ~isempty(ps), p = min(ps); end
    end
end

function [l1, l2, hint, c1, c2] = analysis_labels(key)
    % Per-analysis file-picker labels (l2 empty => hide the 2nd picker) and a hint.
    % c1/c2 are the event-count file-picker labels; both empty => no count files.
    c1 = ''; c2 = '';
    switch key
        case 'anova1'
            l1='Data File:'; l2=''; hint='Compares 3+ groups. Pick the column that labels the groups.';
        case 'ancova'
            l1='Data File:'; l2=''; hint='Compares groups while controlling covariates. Set Group and Covariate cols.';
        case 'regression'
            l1='Data File:'; l2=''; hint='Relates a continuous measure to EEG. Set the Predictor (and any covariates).';
        case 'rmanova'
            l1='Data File:'; l2=''; hint='Compares 3+ within-subject conditions (same subjects across conditions). Set Subject and Condition.';
        case 'mixed2way'
            l1='Data File:'; l2=''; hint=['Two-way mixed ANOVA: does the condition effect differ between groups? Every subject needs ' ...
                  'exactly one row for each condition. Set Group, Subject and Condition.'];
        case 'pairedT'
            l1='Condition A File:'; l2='Condition B File:'; hint='Two conditions, same subjects. Load the A and B files (warm = A > B).';
        case 'onesampleT'
            l1='Condition File:'; l2=''; hint='Tests one condition against 0 (e.g. an overnight change score). Load a single file.';
        case 'unpairedT'
            l1='Group A File:'; l2='Group B File:'; hint='Two independent groups. Load the group A and B files.';
        case 'correlationP'
            l1='Measure 1 File (e.g. EEG):'; l2='Measure 2 File (e.g. behaviour):';
            hint=['Linear correlation across subjects. Both files need a Subject column. ' ...
                  'Measure 2 is either a full channels file (correlated channel by channel) ' ...
                  'or a single behavioural column ' char(8212) ' e.g. Subject,MRRT, one value ' ...
                  'per subject, compared against every channel. With a single score, leave ' ...
                  'Data Type on absolute: normalize re-scales each subject''s channels to ' ...
                  'average zero, which forces as many positive as negative channels.'];
        case 'correlationS'
            l1='Measure 1 File (e.g. EEG):'; l2='Measure 2 File (e.g. behaviour):';
            hint=['Rank correlation (robust to outliers). Both files need a Subject column. ' ...
                  'Measure 2 is either a full channels file (correlated channel by channel) ' ...
                  'or a single behavioural column ' char(8212) ' e.g. Subject,MRRT, one value ' ...
                  'per subject, compared against every channel. With a single score, leave ' ...
                  'Data Type on absolute: normalize re-scales each subject''s channels to ' ...
                  'average zero, which forces as many positive as negative channels.'];
        case {'circ_phase_group','circ_phase_group_u2'}
            % Wording is deliberately group-explicit. The retired labels said
            % "Load both condition files", which invited the single worst
            % mistake available here: running within-subject paired data
            % through an independent-samples test.
            l1='Group A Angles File:'; l2='Group B Angles File:';
            c1='Group A Event Count File:'; c2='Group B Event Count File:';
            hint=['Two INDEPENDENT groups of subjects (not two conditions of the same subjects). ' ...
                  'Each angle file holds one preferred coupling phase per subject per channel. ' ...
                  'The matching event-count file holds the number of detected events behind each ' ...
                  'of those angles ' char(8212) ' it is required, and controls for the fact that ' ...
                  'fewer events give a noisier, wider-looking phase estimate.'];
            if strcmp(key, 'circ_phase_group_u2')
                hint=[hint '   Watson U2 is the alternative test and needs full-precision ' ...
                      'angles: a whole-degree export (Luna COUPL_ANGLE, TurtleWave''s degree ' ...
                      'column) is rejected. Prefer the recommended test unless you have a ' ...
                      'specific reason.'];
            end
        case 'circ_corrAngLinear'
            l1='Angles File:'; l2='Measure File (linear):';
            hint=['Circular-linear correlation: preferred phase versus a linear behavioural ' ...
                  'measure, PAIRED WITHIN SUBJECT. Both files need a matching Subject column. ' ...
                  'No event-count file is used.'];
        otherwise
            l1='Data 1 File:'; l2='Data 2 File:'; hint='';
    end
end

function B = normalizeDetect(stored)
% A snpm_detect_behavior_table result that is always safe to read. Slots that
% have never been filled, and files whose headers could not be read at all,
% come back as kind 'none' with empty lists rather than [] -- so no caller has
% to test isempty before looking at .kind or .data_cols.
    if isstruct(stored) && isfield(stored, 'kind')
        B = stored;
        if ~isfield(B, 'error_message'), B.error_message = ''; end
        return;
    end
    B = struct('kind', 'none', 'is_behavior', false, 'col_name', '', ...
               'channel_cols', {{}}, 'nonchannel_cols', {{}}, 'data_cols', {{}}, ...
               'subject_cols', {{}}, 'channel_fraction', 0, 'error_message', '');
end

function labels = montageOrderLabels(key)
% The channel labels the engine will actually match Data 1's columns against:
% the montage's LOADED chanlocs, exactly as core_snpm_analysis does for its
% order guard, with the registry's own list as the documented fallback. Cached
% per recording system because the readiness gate runs on every widget change
% and this would otherwise hit the disk each time.
    persistent CACHE
    if isempty(CACHE), CACHE = struct(); end
    fld = matlab.lang.makeValidName(char(key));
    if isfield(CACHE, fld), labels = CACHE.(fld); return; end
    m = snpm_montage_registry(key);
    try
        cl = load(m.chanloc_file);
        labels = {cl.(m.chanloc_var).labels};
    catch
        labels = m.labels;
    end
    CACHE.(fld) = labels;
end

function fillRoleDropdown(dd, meta, def)
    % Populate a role dropdown from metadata column names; '' = unselected.
    dd.Items     = [{'(select)'}, meta(:)'];
    dd.ItemsData = [{''},        meta(:)'];
    if ~isempty(def) && any(strcmp(meta, def)), dd.Value = def; else, dd.Value = ''; end
end

function msgs = mixed2way_balance_messages(file, sheet, subjCol, condCol, grpCol)
% Readiness-gate copy of core_snpm_glm.m's require_mixed2way_balance. The
% definition is the engine's, line for line: the columns come from the same
% readtable type detection (detectImportOptions defaults, only the role columns
% selected), levels are stable-unique values of the raw column -- so a missing
% condition (NaN / '' / <missing>) is its own level and a subject is flagged
% whenever its row count in any level is not exactly 1 (missing OR duplicated).
% The engine raises the same error id when a subject sits in more than one
% group; that is checked too when a Group role is set. Any read or type problem
% returns {} so Run stays available and the engine reports.
    msgs = {};
    try
        if strcmp(sheet, 'CSV File')
            o = detectImportOptions(file);
        else
            o = detectImportOptions(file, 'Sheet', sheet);
        end
        cols = {subjCol, condCol};
        if ~isempty(grpCol), cols{end+1} = grpCol; end
        cols = unique(cols, 'stable');
        if ~all(ismember(cols, o.VariableNames)), return; end
        o.SelectedVariableNames = cols;
        T = readtable(file, o);

        [subj, ~, si] = unique(T.(subjCol), 'stable');
        [cond, ~, ci] = unique(T.(condCol), 'stable');
        counts = accumarray([si(:), ci(:)], 1, [numel(subj), numel(cond)]);
        bad = find(any(counts ~= 1, 2));
        if ~isempty(bad)
            % An empty Condition cell is its own "condition" to the engine, so
            % one blank cell flags nearly every subject. Name the real cause.
            nEmpty = nnz(ismissing(T.(condCol)));
            if nEmpty > 0
                msgs{end+1} = sprintf(['The %s column has %d empty cell(s); every row ' ...
                    'needs a condition.'], condCol, nEmpty);
            end
            msgs{end+1} = sprintf(['%s not have exactly one row for every condition ' ...
                '(missing or repeated): %s. Complete or remove them; unbalanced data need ' ...
                'the mixed model (scripts/run_lmm_analysis.m).'], ...
                count_phrase(numel(bad)), id_list(subj, bad));
        end
        if ~isempty(grpCol)
            [~, ~, gi] = unique(T.(grpCol), 'stable');
            ngrp = accumarray(si(:), gi(:), [numel(subj), 1], @(v) numel(unique(v)));
            multi = find(ngrp > 1);
            if ~isempty(multi)
                if isscalar(multi), who = '1 subject appears'; else, who = sprintf('%d subjects appear', numel(multi)); end
                msgs{end+1} = sprintf(['%s in more than one group: %s. Each subject ' ...
                    'must belong to exactly one group.'], who, id_list(subj, multi));
            end
        end
    catch
        msgs = {};
    end
end

function s = count_phrase(n)
    if n == 1, s = '1 subject does'; else, s = sprintf('%d subjects do', n); end
end

function s = id_list(levels, idx)
% Up to three subject IDs, formatted like the engine's level_to_str, then
% "and N more" if there are more.
    k = min(numel(idx), 3);
    names = cell(1, k);
    for j = 1:k
        if iscell(levels), v = levels{idx(j)}; else, v = levels(idx(j)); end
        if isnumeric(v) || islogical(v), names{j} = num2str(v); else, names{j} = char(string(v)); end
    end
    s = strjoin(names, ', ');
    if numel(idx) > k, s = sprintf('%s and %d more', s, numel(idx) - k); end
end

function n = firstMatch(meta, pat)
    % First metadata column whose name contains pat (case-insensitive), else ''.
    if isempty(meta), n = ''; return; end
    idx = find(contains(lower(meta), lower(pat)), 1);
    if isempty(idx), n = ''; else, n = meta{idx}; end
end
