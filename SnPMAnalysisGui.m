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
        output_path = ''
        covariate_file = ''  % Store covariate file path
        snpm_toolbox_path = ''
        eeglab_path = ''
        last_results = ''  % Store last results for PDF export
        lastComparison = 'pairedT'  % last real (non-header) comparison key selected
        RoleRows = {}  % {key, labelHandle, fieldHandle} for show/hide of GLM column-role fields
        DetectedNChan = 0  % channel columns detected in the last loaded file (for montage validation)
        DetectedChanNames = {}  % channel column NAMES detected (for label-presence validation)
        SpectralFolder = ''  % selected BIDS spectral folder (spectral mode)
        SubjectsCsv = ''  % optional subject-metadata CSV joined on Subject (spectral mode)
        FullComparisonItems = {}  % saved ComparisonDropDown Items/ItemsData to restore when leaving spectral mode
        FullComparisonItemsData = {}
    end

    % Row indices of the collapsible auxiliary panels in app.GridLayout.
    % They sit AFTER the numbered sections (1-4 are always shown, in order).
    properties (Access = private, Constant)
        RolesRow = 8
        CovRow   = 9
        RolesRowH = 170
        CovRowH   = 64
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
                detectAndFill(app, app.data1_file, app.Data1SheetDropDown.Value, true);
                checkReadyToRun(app);
            end
        end
        function Data2SheetChanged(app, ~)
            if ~isempty(app.data2_file)
                detectAndFill(app, app.data2_file, app.Data2SheetDropDown.Value, false);
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
            if spec
                rh{2} = 0;     rh{3} = 0;     rh{4} = '1x';   % hide file rows; let the spectral grid grow
            else
                rh{2} = 'fit'; rh{3} = 'fit'; rh{4} = 0;      % file rows fit content; collapse spectral grid
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
                orh = app.GridLayout.RowHeight;
                orh{app.RolesRow} = 0; orh{app.CovRow} = 0;
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
                app.DetectCaptionLabel.Text = sprintf('Detected: %d channels, %d bands, %d stages', ...
                    numel(info.chanlabels), numel(info.bands), numel(info.stages));
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
            [l1, l2, hint] = analysis_labels(key);
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
            app.RolesPanel.Visible     = isGLM;
            app.CovariatePanel.Visible = isCorr;

            % Collapse the rows of the hidden panels so there is no empty gap;
            % a visible auxiliary panel sizes to its content ('fit').
            rh = app.GridLayout.RowHeight;
            if isGLM,  rh{app.RolesRow} = 'fit'; else, rh{app.RolesRow} = 0; end
            if isCorr, rh{app.CovRow}   = 'fit'; else, rh{app.CovRow}   = 0; end
            app.GridLayout.RowHeight = rh;
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
        end

        % Read a file's headers, classify channel vs metadata columns, fill the UI
        function meta = detectAndFill(app, file, sheet, isPrimary)
            meta = {};
            try
                if strcmp(sheet, 'CSV File')
                    o = detectImportOptions(file, 'VariableNamingRule','preserve');
                else
                    o = detectImportOptions(file, 'Sheet', sheet, 'VariableNamingRule','preserve');
                end
                vn = o.VariableNames;
            catch
                vn = {};
            end
            % A column is a channel if its name is an EGI E-number, or is in any
            % recording system's label set (so Compumedics names like Fp1 count too).
            reg = snpm_montage_registry();
            allLabels = [reg.labels];
            ischan = cellfun(@(s) ~isempty(regexp(s,'^E\d+$','once')) || strcmpi(s,'Cz'), vn) ...
                     | ismember(vn, allLabels);
            chans = vn(ischan);
            meta  = vn(~ischan);
            app.DetectedNChan = numel(chans);
            app.DetectedChanNames = chans;
            app.DetectCaptionLabel.Text = sprintf('Detected: %d channels, %d metadata columns', numel(chans), numel(meta));
            if isPrimary
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

                detectAndFill(app, app.data1_file, app.Data1SheetDropDown.Value, true);
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

                detectAndFill(app, app.data2_file, app.Data2SheetDropDown.Value, false);
                checkReadyToRun(app);
            end
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
                app.ResultsTextArea.Value = {['Error: ' ME.message], '', 'Stack trace:', ME.getReport()};
                uialert(app.UIFigure, ME.message, 'Analysis Error');
            end

            % Re-enable the run button
            app.RunAnalysisButton.Enable = 'on';
        end

        % Button pushed function: ResetButton
        function ResetButtonPushed(app, event)
            % Reset all fields to default
            app.data1_file = '';
            app.data2_file = '';
            app.output_path = '';
            app.covariate_file = '';

            app.Data1FileLabel.Text = 'No file selected';
            app.Data2FileLabel.Text = 'No file selected';
            app.OutputPathLabel.Text = 'No path selected';
            app.CovariateFileLabel.Text = 'No covariate file selected';
            app.UseCovariatesCheckBox.Value = false;

            % clear GLM column-role controls + detection state
            initRolePickers(app);

            app.Data1SheetDropDown.Items = {};
            app.Data2SheetDropDown.Items = {};

            app.ChannelsDropDown.Value = 'egi';
            app.DataTypeDropDown.Value = 'absolute';
            app.ComparisonDropDown.Value = 'pairedT';
            app.lastComparison = 'pairedT';
            updateRoleFields(app);
            applyAnalysisLayout(app);
            app.TailsDropDown.Value = 'both';
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
            [~, l2] = analysis_labels(key);
            twoFile = ~isempty(l2);

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
            end

            % Output folder
            if isempty(app.output_path), missing{end+1} = 'Choose an output folder'; end

            % File must contain channels of the selected recording system
            if app.DetectedNChan > 0
                m = snpm_montage_registry(app.ChannelsDropDown.Value);
                present = nnz(ismember(m.labels, app.DetectedChanNames));
                if present == 0
                    missing{end+1} = sprintf('No %s channels found in the file (expected labels like %s)', ...
                        m.display, m.labels{min(2,numel(m.labels))});
                else
                    app.DetectCaptionLabel.Text = sprintf('Detected: %d channels  (matches %d of %d %s channels)', ...
                        app.DetectedNChan, present, numel(m.labels), m.display);
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
            % rows: 1 title, 2 subtitle, 3 note, 4 analysis(1), 5 data(2), 6 params(3),
            %       7 output(4), 8 roles(GLM), 9 covariate(corr), 10 validation,
            %       11 buttons, 12 status, 13 results
            % Sections 1-4 always render in order; roles/covariate are auxiliary boxes below them.
            % Weighted rows so panels grow with the window instead of clipping.
            % 'fit' = size to content (compact control rows); '1x' / '2x' = expand.
            % rows: 1 title, 2 subtitle, 3 note, 4 analysis, 5 DATA(expands),
            %       6 params, 7 output, 8 roles, 9 covariate, 10 validation,
            %       11 buttons, 12 status, 13 RESULTS(expands).
            % Roles(8)/covariate(9) start collapsed (0) and are sized by the layout helpers.
            % Last row (Results) is a fixed height, not '1x': a '1x' row would
            % expand to fill and suppress overflow, so the grid would never
            % scroll. Fixed height keeps the total content height determinate.
            app.GridLayout.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 0, 0, 'fit', 'fit', 'fit', 180};
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
                '   Wheeler-Watson test', ...
                '   Watson U2 test'};
            app.ComparisonDropDown.ItemsData = { ...
                '__h1__', 'unpairedT', 'pairedT', 'onesampleT', ...
                '__h2__', 'anova1', 'ancova', 'rmanova', 'mixed2way', ...
                '__h3__', 'correlationP', 'correlationS', 'regression', ...
                '__h4__', 'circ_wheeler_watson_Test', 'circ_WatsonsU2Test'};
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
            % rows: 1 source toggle, 2 file-1, 3 file-2, 4 spectral grid (expands), 5 caption.
            % Row 4 starts collapsed; the layout helper grows it in spectral mode.
            dataFilesGrid.RowHeight = {'fit', 'fit', 'fit', 0, 'fit'};
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

            % --- Spectral-folder controls (row 4 nested grid): shown in 'Spectral folder' mode ---
            app.SpecGrid = uigridlayout(dataFilesGrid);
            app.SpecGrid.Layout.Row = 4; app.SpecGrid.Layout.Column = [1 5];
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
            app.DetectCaptionLabel.Layout.Row = 5;
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
            app.ValidationPanel.Layout.Row = 10;
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
            buttonGrid.Layout.Row = 11;
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
            app.StatusLabel.Layout.Row = 12;
            app.StatusLabel.Layout.Column = 1;
            app.StatusLabel.HorizontalAlignment = 'center';
            app.StatusLabel.FontAngle = 'italic';

            % Results Panel
            app.ResultsPanel = uipanel(app.GridLayout);
            app.ResultsPanel.Layout.Row = 13;
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

function [l1, l2, hint] = analysis_labels(key)
    % Per-analysis file-picker labels (l2 empty => hide the 2nd picker) and a hint.
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
            l1='Data File:'; l2=''; hint='Two-way mixed ANOVA: does the condition effect differ between groups? Set Group, Subject and Condition.';
        case 'pairedT'
            l1='Condition A File:'; l2='Condition B File:'; hint='Two conditions, same subjects. Load the A and B files (warm = A > B).';
        case 'onesampleT'
            l1='Condition File:'; l2=''; hint='Tests one condition against 0 (e.g. an overnight change score). Load a single file.';
        case 'unpairedT'
            l1='Group A File:'; l2='Group B File:'; hint='Two independent groups. Load the group A and B files.';
        case 'correlationP'
            l1='Measure 1 File (e.g. EEG):'; l2='Measure 2 File (e.g. behaviour):'; hint='Linear correlation across subjects, channel by channel. Both files need a Subject column.';
        case 'correlationS'
            l1='Measure 1 File (e.g. EEG):'; l2='Measure 2 File (e.g. behaviour):'; hint='Rank correlation (robust to outliers). Both files need a Subject column.';
        case {'circ_wheeler_watson_Test','circ_WatsonsU2Test'}
            l1='Angles A File (radians):'; l2='Angles B File (radians):'; hint='Circular statistics for phase/angle data in radians. Load both condition files.';
        otherwise
            l1='Data 1 File:'; l2='Data 2 File:'; hint='';
    end
end

function fillRoleDropdown(dd, meta, def)
    % Populate a role dropdown from metadata column names; '' = unselected.
    dd.Items     = [{'(select)'}, meta(:)'];
    dd.ItemsData = [{''},        meta(:)'];
    if ~isempty(def) && any(strcmp(meta, def)), dd.Value = def; else, dd.Value = ''; end
end

function n = firstMatch(meta, pat)
    % First metadata column whose name contains pat (case-insensitive), else ''.
    if isempty(meta), n = ''; return; end
    idx = find(contains(lower(meta), lower(pat)), 1);
    if isempty(idx), n = ''; else, n = meta{idx}; end
end
