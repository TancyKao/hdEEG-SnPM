function [results_struct, results_text] = core_snpm_analysis(params)
% Core SnPM analysis extracted from SnPM_main4csv.m with Covariate Support
%
% INPUT:
%   params - structure with analysis parameters:
%       .data1_file, .data2_file - Data file paths
%       .data1_sheet, .data2_sheet - Sheet names
%       .output_path - Output directory
%       .channels - Channel configuration
%       .datatype - Data type (absolute/logscale/normalize)
%       .comparison - Analysis type (pairedT/unpairedT/correlationS/correlationP/etc.)
%       .tail - Test tail (both/left/right)
%       .permutations - Number of permutations
%       .covariate_file - Path to covariate CSV file (optional)
%       .use_covariates - Boolean flag to use covariates (optional)
%
% OUTPUT:
%   results_struct - structured results
%   results_text - formatted text results

    % Add paths (make these configurable)
    % Use current working directory as base path for multi-user/multi-PC compatibility
    current_path = pwd;
    
    if isfield(params, 'snpm_path') && ~isempty(params.snpm_path)
        addpath(genpath(params.snpm_path));
        fprintf('Added SnPM path: %s\n', params.snpm_path);
    else
        addpath(genpath(current_path), '-end');
        fprintf('Added SnPM path (current directory): %s\n', current_path);
    end
    
    if isfield(params, 'eeglab_path') && ~isempty(params.eeglab_path)
        addpath(params.eeglab_path);
        fprintf('Added EEGLAB path: %s\n', params.eeglab_path);
    else
        % EEGLAB is located in the SnPM_2025 folder
        eeglab_path = fullfile(current_path, 'eeglab2022.1');
        if exist(eeglab_path, 'dir')
            addpath(genpath(eeglab_path), '-end');
            fprintf('Added EEGLAB path: %s\n', eeglab_path);
        else
            warning('EEGLAB folder not found at: %s', eeglab_path);
        end
    end
    
    % Mixed linear model path uses trial-level (long-format) data and a
    % different statistical engine; route it to its own pipeline. Paths have
    % already been set above, so core_snpm_lmm can rely on them.
    if strcmp(params.comparison, 'mixedmodel')
        [results_struct, results_text] = core_snpm_lmm(params);
        return
    end

    % GLM group-analysis presets (>2 groups, ANCOVA, regression, RM-ANOVA,
    % group x condition) route to their own pipeline; legacy tests below
    % are untouched.
    if ismember(params.comparison, {'anova1', 'ancova', 'regression', 'rmanova', 'mixed2way'})
        [results_struct, results_text] = core_snpm_glm(params);
        return
    end

    % Extract covariate parameters (with defaults)
    if isfield(params, 'covariate_file')
        covariate_file = params.covariate_file;
    else
        covariate_file = '';
    end
    
    if isfield(params, 'use_covariates')
        use_covariates = params.use_covariates;
    else
        use_covariates = false;
    end
    
    % Load data tables
    if strcmp(params.data1_sheet, 'CSV File')
        data1Table = readtable(params.data1_file);
    else
        data1Table = readtable(params.data1_file, 'Sheet', params.data1_sheet);
    end
    
    % One-sample t-test needs only one condition (tested against 0). When no
    % Data 2 is supplied, data_y becomes a zero matrix (set after channel
    % selection, below): the onesampleTboth permutation is a sign-flip test, so
    % data_x vs 0 is exactly a one-sample test. A copy of data1Table carries the
    % right shape/columns through preprocessing. A supplied Data 2 keeps the
    % legacy paired-difference-vs-0 behaviour.
    onesample_vs_zero = false;
    if strcmp(params.comparison, 'onesampleT') && ...
            (~isfield(params, 'data2_file') || isempty(params.data2_file))
        data2Table = data1Table;
        onesample_vs_zero = true;
        if ~isfield(params, 'data2_sheet') || isempty(params.data2_sheet)
            params.data2_sheet = 'zero';   % keeps downstream sprintf/labels tidy
        end
    elseif strcmp(params.data2_sheet, 'CSV File')
        data2Table = readtable(params.data2_file);
    else
        data2Table = readtable(params.data2_file, 'Sheet', params.data2_sheet);
    end

    %% Load and Process Covariates
    covariates = [];
    covariate_names = {};
    covariate_subjects = {};
    
    if use_covariates && ~isempty(covariate_file) && exist(covariate_file, 'file')
        try
            fprintf('Loading covariates from: %s\n', covariate_file);
            covariate_table = readtable(covariate_file);
            
            % Assume first column is Subject ID, rest are covariates
            covariate_subjects_raw = covariate_table{:, 1};
            covariate_data = covariate_table{:, 2:end};
            covariate_names = covariate_table.Properties.VariableNames(2:end);
            
            % Convert subject IDs to strings for matching
            if isnumeric(covariate_subjects_raw)
                covariate_subjects = arrayfun(@(x) sprintf('sub%03d', x), covariate_subjects_raw, 'UniformOutput', false);
            else
                covariate_subjects = cellstr(string(covariate_subjects_raw));
            end
            
            fprintf('Loaded %d covariates for %d subjects: %s\n', ...
                size(covariate_data, 2), size(covariate_data, 1), strjoin(covariate_names, ', '));
            
        catch ME
            warning('SnPM:CovariateLoadError', 'Error loading covariates: %s. Proceeding without covariates.', ME.message);
            use_covariates = false;
            covariates = [];
            covariate_names = {};
            covariate_subjects = {};
        end
    end
    
    
    if contains(params.comparison, 'correlation')
        % Find subject columns
        subj_col1 = find(contains(lower(data1Table.Properties.VariableNames), 'subject'), 1);
        subj_col2 = find(contains(lower(data2Table.Properties.VariableNames), 'subject'), 1);
        
        if isempty(subj_col1) || isempty(subj_col2)
            error('Correlation analysis requires Subject columns in both files');
        end
        
        % Get subject IDs
        subjects1_raw = data1Table{:, subj_col1};
        subjects2_raw = data2Table{:, subj_col2};
        
        % Convert to strings for consistent matching
        % Handle both numeric and string subject IDs
        if isnumeric(subjects1_raw)
            subjects1 = arrayfun(@(x) sprintf('sub%03d', x), subjects1_raw, 'UniformOutput', false);
            subjects1 = string(subjects1);
        elseif iscell(subjects1_raw)
            subjects1 = string(subjects1_raw);
        else
            subjects1 = string(subjects1_raw);
        end
        
        if isnumeric(subjects2_raw)
            subjects2 = arrayfun(@(x) sprintf('sub%03d', x), subjects2_raw, 'UniformOutput', false);
            subjects2 = string(subjects2);
        elseif iscell(subjects2_raw)
            subjects2 = string(subjects2_raw);
        else
            subjects2 = string(subjects2_raw);
        end
        
        % Debug: Print subject IDs
        fprintf('Data1 subjects (%d): %s\n', length(subjects1), strjoin(subjects1(1:min(5,end)), ', '));
        fprintf('Data2 subjects (%d): %s\n', length(subjects2), strjoin(subjects2(1:min(5,end)), ', '));
        
        % Find common subjects and their indices
        [common_subjects, idx1, idx2] = intersect(subjects1, subjects2, 'stable');
        
        % Debug: Print matching results
        fprintf('Found %d common subjects\n', length(common_subjects));
        if length(common_subjects) > 0
            fprintf('Common subjects: %s\n', strjoin(common_subjects(1:min(5,end)), ', '));
        end
        
        % Report subject matching status
        original_n1 = length(subjects1);
        original_n2 = length(subjects2);
        matched_n = length(common_subjects);
        
        data1Table_matched = data1Table(idx1, :);
        data2Table_matched = data2Table(idx2, :);

        % Get data columns (excluding subject columns)
        data_cols1 = ~contains(lower(data1Table_matched.Properties.VariableNames), 'subject');
        data_cols2 = ~contains(lower(data2Table_matched.Properties.VariableNames), 'subject');
        
        % Extract data matrices for NaN checking
        data_matrix1 = data1Table_matched{:, data_cols1};
        data_matrix2 = data2Table_matched{:, data_cols2};
        
        % For correlation analysis, we use pairwise deletion (per channel)
        % Do NOT remove subjects here - NaN handling happens during correlation computation
        % This allows each channel to use all available subject pairs
        
        data1Table = data1Table_matched;
        data2Table = data2Table_matched;
        
        % Report NaN statistics for information
        n_nan_channels1 = sum(isnan(data_matrix1), 2);  % NaN count per subject in data1
        n_nan_channels2 = sum(isnan(data_matrix2), 2);  % NaN count per subject in data2
        total_channels = size(data_matrix1, 2);
        
        fprintf('NaN Statistics:\n');
        fprintf('  Data1: %.1f%% NaN channels per subject (mean), max=%d/%d channels\n', ...
            mean(n_nan_channels1)/total_channels*100, max(n_nan_channels1), total_channels);
        fprintf('  Data2: %.1f%% NaN channels per subject (mean), max=%d/%d channels\n', ...
            mean(n_nan_channels2)/total_channels*100, max(n_nan_channels2), total_channels);
        
        % Count valid pairs per channel (for reporting)
        valid_pairs_per_channel = zeros(1, total_channels);
        for ch = 1:total_channels
            valid_idx = ~isnan(data_matrix1(:, ch)) & ~isnan(data_matrix2(:, ch));
            valid_pairs_per_channel(ch) = sum(valid_idx);
        end
        fprintf('  Valid subject pairs per channel: min=%d, mean=%.1f, max=%d\n', ...
            min(valid_pairs_per_channel), mean(valid_pairs_per_channel), max(valid_pairs_per_channel));
        
        % Set final_n to matched_n (all matched subjects will be used)
        final_n = matched_n;
        nan_removed_n = 0;  % No subjects removed at this stage

        if matched_n < original_n1 || matched_n < original_n2
            missing1 = setdiff(subjects1, common_subjects);
            missing2 = setdiff(subjects2, common_subjects);
            
            warning_msg = sprintf(['Subject matching for correlation analysis:\n' ...
                                  'Original subjects in Data1: %d, Data2: %d\n' ...
                                  'Matched subjects: %d\n' ...
                                  'Missing from Data1: %s\n' ...
                                  'Missing from Data2: %s\n' ...
                                  'Note: NaN values handled per-channel during correlation'], ...
                                  original_n1, original_n2, matched_n, ...
                                  strjoin(missing1, ', '), strjoin(missing2, ', '));
            warning(warning_msg);
        end

        
        if final_n < 3
            error('Correlation analysis requires at least 3 matching subjects. Found: %d', final_n);
        end
        
        fprintf('Correlation analysis: Using %d matched subjects (NaN handled per-channel)\n', final_n);

        
        % Store subject matching info for results
        subject_match_info = struct();
        subject_match_info.original_n1 = original_n1;
        subject_match_info.original_n2 = original_n2;
        subject_match_info.matched_n = matched_n;
        subject_match_info.nan_removed_n = nan_removed_n;
        subject_match_info.final_n = final_n;
        subject_match_info.common_subjects = common_subjects;  % All common subjects (no filtering)
        
        % Match covariates to common subjects if using covariates
        if use_covariates && ~isempty(covariate_subjects)
            % Convert common_subjects to cell array of strings for matching
            common_subjects_str = cellstr(string(common_subjects));  % Use all common subjects
            
            [matched_subjects, idx_common, idx_cov] = intersect(common_subjects_str, covariate_subjects, 'stable');
            
            if length(matched_subjects) < length(common_subjects_str)
                warning('Not all subjects found in covariate file. %d subjects will be excluded.', ...
                    length(common_subjects_str) - length(matched_subjects));
                
                % Update data tables to only include subjects with covariates
                data1Table = data1Table(idx_common, :);
                data2Table = data2Table(idx_common, :);
                
                % Update subject match info
                subject_match_info.covariate_matched_n = length(matched_subjects);
                subject_match_info.covariate_excluded_n = length(common_subjects_str) - length(matched_subjects);
            else
                subject_match_info.covariate_matched_n = length(matched_subjects);
                subject_match_info.covariate_excluded_n = 0;
            end
            
            % Reorder covariates to match subjects
            covariates = covariate_data(idx_cov, :);
            
            fprintf('Matched %d subjects with covariates (%d excluded)\n', ...
                length(matched_subjects), length(common_subjects_str) - length(matched_subjects));
        end
    end


    % Remove subject columns 
    varToRemove1 = contains(lower(data1Table.Properties.VariableNames), lower('Subject'));
    data1Table(:, varToRemove1) = [];
    varToRemove2 = contains(lower(data2Table.Properties.VariableNames), lower('Subject'));
    data2Table(:, varToRemove2) = [];
    
    % Create output directory
    if ~exist(params.output_path, 'dir')
        mkdir(params.output_path);
    end
    
    % Channel setup (shared label-based montage selection; recording system
    % chosen via params.channels -> snpm_montage_registry).
    data_x = data1Table{:,:};
    data_y = data2Table{:,:};
    channel_labels = data1Table.Properties.VariableNames;
    chsel = snpm_setup_channels(params.channels, {data_x, data_y}, channel_labels);
    data_x = chsel.data{1};
    data_y = chsel.data{2};
    if onesample_vs_zero, data_y = zeros(size(data_x)); end   % one condition vs 0
    chanlocs        = chsel.chanlocs;
    neighbors       = chsel.neighbors;
    insidegoodch    = chsel.insidegoodch;
    channel_mapping = chsel.channel_mapping; %#ok<NASGU>
    select_mode     = chsel.select_mode;
    fprintf('%s montage: %d/%d data columns matched to channel locations\n', ...
        params.channels, numel(insidegoodch), numel(channel_labels));

    % Set analysis parameters from GUI
    E = 0.5;
    H = 2;
    alpha = 0.05;
    
    %Apply data type transformation
    switch params.datatype
        case 'absolute'
            % No transformation
        case 'logscale'
            data_x = log10(data_x);
            data_y = log10(data_y);
        case 'normalize'
            data_x = bsxfun(@minus,data_x,mean(data_x,2,'omitnan'));
            data_x = bsxfun(@rdivide,data_x,std(data_x,[],2,'omitnan'));
            data_y = bsxfun(@minus,data_y,mean(data_y,2,'omitnan'));
            data_y = bsxfun(@rdivide,data_y,std(data_y,[],2,'omitnan'));
    end
    
    
    % Calculate degrees of freedom
    g1subj_n = size(data1Table,1);
    g2subj_n = size(data2Table,1);
   
    
    if ismember(params.comparison,{'pairedT', 'unpairedT', 'onesampleT'})
        degrees_of_freedom = min(g1subj_n-1,g2subj_n-1);
        threshold = tinv(1 - alpha/2,degrees_of_freedom);
    elseif ismember(params.comparison,{'correlationP', 'correlationS'})
        n_subjects = min(g1subj_n, g2subj_n);
        
        % Adaptive threshold based on sample size for cluster analysis
        % This is used for identifying spatial clusters, not for statistical inference
        % (which is done via permutation testing)
        if n_subjects <= 15
            threshold = 0.50;  % Conservative for small samples
        elseif n_subjects <= 30
            threshold = 0.40;  % Moderate for medium samples
        elseif n_subjects <= 50
            threshold = 0.35;  % More liberal for larger samples
        else
            threshold = 0.30;  % Liberal for very large samples
        end
        
        fprintf('Correlation cluster threshold: n=%d subjects, r_threshold=%.2f (adaptive)\n', ...
            n_subjects, threshold);
        fprintf('Note: This threshold is for cluster identification only.\n');
        fprintf('      Statistical significance is determined by permutation testing.\n');
    end
    
    % Remove subjects with all NaN channels (must remove from BOTH datasets to keep them aligned)
    fprintf('Before NaN removal: data_x=%dx%d, data_y=%dx%d\n', size(data_x,1), size(data_x,2), size(data_y,1), size(data_y,2));
    
    % Find subjects with all NaN in either dataset
    all_nan_x = all(isnan(data_x), 2);
    all_nan_y = all(isnan(data_y), 2);
    subjects_to_remove = all_nan_x | all_nan_y;
    
    if any(subjects_to_remove)
        fprintf('Removing %d subjects with all NaN channels\n', sum(subjects_to_remove));
        data_x(subjects_to_remove, :) = [];
        data_y(subjects_to_remove, :) = [];
    end
    
    fprintf('After NaN removal: data_x=%dx%d, data_y=%dx%d\n', size(data_x,1), size(data_x,2), size(data_y,1), size(data_y,2));
    
    % Verify dimensions match for paired/correlation tests
    if contains(params.comparison, 'correlation') || contains(params.comparison, 'pairedT')
        if size(data_x, 1) ~= size(data_y, 1)
            error('Data dimension mismatch after NaN removal: data_x has %d rows, data_y has %d rows', ...
                size(data_x, 1), size(data_y, 1));
        end
    end
    
    % Global statistical analysis
    if ~contains(params.comparison, 'circ')
        res = global_stat_test(data_x, data_y, alpha, params.comparison, params.tail);
    end
    
    % Apply partial correlation if covariates are provided for correlation analyses
    if contains(params.comparison, 'correlation') && use_covariates && ~isempty(covariates)
        fprintf('Applying partial correlation with %d covariate(s)\n', size(covariates, 2));
        
        % Apply residualization to both data matrices
        % This removes the linear effects of covariates from the data
        data_x_residual = zeros(size(data_x));
        data_y_residual = zeros(size(data_y));
        
        for ch = 1:size(data_x, 2)
            % Get valid indices for this channel
            valid_idx = ~isnan(data_x(:, ch)) & ~isnan(data_y(:, ch));
            
            if sum(valid_idx) > size(covariates, 2) + 1  % Need enough subjects
                % Residualize data_x for this channel
                X_design = [ones(sum(valid_idx), 1), covariates(valid_idx, :)];
                beta_x = X_design \ data_x(valid_idx, ch);
                data_x_residual(valid_idx, ch) = data_x(valid_idx, ch) - X_design * beta_x;
                
                % Residualize data_y for this channel
                beta_y = X_design \ data_y(valid_idx, ch);
                data_y_residual(valid_idx, ch) = data_y(valid_idx, ch) - X_design * beta_y;
            else
                % Not enough subjects - keep as NaN
                data_x_residual(:, ch) = NaN;
                data_y_residual(:, ch) = NaN;
            end
        end
        
        % Use residualized data for permutation analysis
        data_x = data_x_residual;
        data_y = data_y_residual;
        
        fprintf('Residualization complete. Running permutation analysis on residuals.\n');
    end
    
    % Main permutation analysis
    if ~contains(params.comparison, 'circ')
        [T,p,~,~] = snpm_single_threshold_with_TFCE(data_x,data_y,neighbors,E,H,alpha,params.comparison,params.tail,params.permutations);
    else
        [T,p,~,~] = snpm_single_threshold_with_TFCE_circ(data_x,data_y,neighbors,E,H,alpha,params.comparison,params.tail,params.permutations);
    end
    
    % Cluster analysis
    if ~contains(params.comparison, 'circ')
        [Clusters] = snpm_cluster_analysis(data_x,data_y,threshold,neighbors,alpha,params.comparison,params.tail,params.permutations);
    else
        [Clusters] = snpm_cluster_analysis_circ(data_x,data_y,threshold,neighbors,alpha,params.comparison,params.tail,params.permutations);
    end
    
    % Find significant channels/clusters
    if contains(params.comparison, {'circ_wheeler_watson_Test', 'circ_WatsonsU2Test'})
        uncorrsigch = find(T.real_T>=threshold);
        correctTFCEsigch = find(p.correctedTFCE <=alpha);
    else
        uncorrsigch = find(p.real<=alpha);
        correctTFCEsigch = find(p.correctedTFCE <=alpha);
    end
    
    sigclusters = find([Clusters.p] <=alpha);
    SnPMsigch = [Clusters(sigclusters).channels];
    
    

    %% Plot figures
    % Extract basenames from file paths
    [~, data1_basename] = fileparts(params.data1_file);
    [~, data2_basename] = fileparts(params.data2_file);
    
    % Simple labels for plots
    data1_label = 'Data1';
    data2_label = 'Data2';
    
    data1_basename = sprintf('%s_%s',data1_basename, params.data1_sheet);
    data2_basename = sprintf('%s_%s',data2_basename, params.data2_sheet);
    if onesample_vs_zero, data2_basename = 'zero'; end   % single condition vs 0

    % Create filename with basenames and sheet names
    base_filename = sprintf('%s VS %s', ...
        data1_basename, data2_basename ); % use VS to separate two lines in figure title
    
    % Set savepath for plots
    savepath = params.output_path;
    
    % plot individual figure & between two groups/conditions
    
    % only plot individual 
    if onesample_vs_zero
        % Single condition vs 0: condition mean topography + the signed t-vs-0
        % map (T.real_T) as the effect map. The zero "condition B" card is
        % hidden in the report (results_struct.hide_condition_b below).
        avg_data_x = mean(data_x,1,'omitnan');
        avg_data_y = mean(data_y,1,'omitnan');   % all zeros; kept for data_summary
        tv = T.real_T; m = max(abs(tv(:)));
        if m == 0 || ~isfinite(m), m = 1; end
        try
            % skip2: the zero condition-B map is not plotted (it's hidden in the
            % report and constant data breaks topoplot); dat2 is a placeholder.
            plot_topoInd(avg_data_x, zeros(size(avg_data_x)), chanlocs, data1_basename, ...
                data2_basename, savepath, base_filename, 'diffdata', tv, ...
                'cmap_diff', 'jet', 'clim_diff', [-m m], 'skip2', true, ...
                'difftitle', 'Effect (t vs 0)');
        catch ME
            warning(ME.identifier, 'Could not generate one-sample topology plots: %s', ME.message);
        end
    elseif ismember(params.comparison,{'pairedT', 'unpairedT', 'onesampleT'})
        avg_data_x = mean(data_x,1,'omitnan');
        avg_data_y = mean(data_y,1,'omitnan');
        if (avg_data_x~=0) & (avg_data_y ~=0) % avoid the error from one sample t-test
            try
                plot_topoInd(avg_data_x,avg_data_y,chanlocs,data1_basename, data2_basename,savepath,base_filename);

                % Basic usage (uses defaults)
                %plot_topoInd(dat1, dat2, chanlocs, 'data1', 'data2', './output', 'difference');
                
                % Custom colormaps
                %plot_topoInd(dat1, dat2, chanlocs, 'data1', 'data2', './output', 'difference', ...
                %    'cmap1', 'viridis', 'cmap2', 'plasma', 'cmap_diff', 'RdBu');
                
                % Custom color limits
                %plot_topoInd(dat1, dat2, chanlocs, 'data1', 'data2', './output', 'difference', ...
                %    'clim1', [-2 2], 'clim2', [0 5], 'clim_diff', [-1 1]);


            catch ME
                warning(ME.identifier, 'Could not generate individual topology plots: %s', ME.message);
            end
        end
    elseif ismember(params.comparison,{'correlationP', 'correlationS'})
        % Correlation: condition-A/B mean topographies + the per-channel r-map
        % (r = T.real_T for correlation). The report's r-branch expects
        % <data1>_topo.png / <data2>_topo.png / <base>_topo.png; feed the r-map
        % as 'diffdata' so the 3rd map is r, not a meaningless mean difference.
        avg_data_x = mean(data_x,1,'omitnan');
        avg_data_y = mean(data_y,1,'omitnan');
        try
            plot_topoInd(avg_data_x, avg_data_y, chanlocs, data1_basename, data2_basename, ...
                savepath, base_filename, 'diffdata', T.real_T, 'cmap_diff', 'jet', 'clim_diff', [-1 1], ...
                'difftitle', 'Correlation r-map');
        catch ME
            warning(ME.identifier, 'Could not generate correlation topology plots: %s', ME.message);
        end
    elseif ismember(params.comparison,{'circ_wheeler_watson_Test', 'circ_WatsonsU2Test'})
        avg_data_x = circ_mean(data_x);
        avg_data_y = circ_mean(data_y);
        try
            plot_topoCircInd(rad2deg(avg_data_x),rad2deg(avg_data_y),chanlocs,data1_label, data2_label,savepath,base_filename);
        catch ME
            warning(ME.identifier, 'Could not generate circular individual topology plots: %s', E.message);
        end
    end 
    
    % plot t-map or corr-map
    try
        TopoplotSignificant_single(T.real_T,uncorrsigch,correctTFCEsigch,SnPMsigch,chanlocs,...
            insidegoodch,params.comparison, savepath,base_filename, select_mode);
    catch ME
        warning(ME.identifier, 'Could not generate significance topology plots: %s', ME.message);
    end
    
    % Calculate averages for display
    avg_data_x = mean(mean(data_x,'omitnan'),'omitnan');
    avg_data_y = mean(mean(data_y,'omitnan'),'omitnan');
    

    % Create results structure
    results_struct = struct();
    results_struct.T = T;
    results_struct.p = p;
    results_struct.Clusters = Clusters;
    results_struct.uncorrsigch = uncorrsigch;
    results_struct.correctTFCEsigch = correctTFCEsigch;
    results_struct.SnPMsigch = SnPMsigch;
    results_struct.chanlocs = chanlocs;
    results_struct.hide_condition_b = onesample_vs_zero;  % one condition vs 0 -> no B card
    results_struct.data_summary.data1_size = size(data_x);
    results_struct.data_summary.data2_size = size(data_y);
    results_struct.data_summary.data1_mean = avg_data_x;
    results_struct.data_summary.data2_mean = avg_data_y;
    results_struct.data_summary.data1_std = std(mean(data_x,2,'omitnan'));
    results_struct.data_summary.data2_std = std(mean(data_y,2,'omitnan'));
    
    % Save to file

    % Save results
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    %figtitle = [params.comparison '_' params.tail '_' timestamp];
    outputSname = fullfile(params.output_path, [base_filename, '_' timestamp]);

    save([outputSname '.mat'], 'results_struct');
    
    % Generate Excel output using func_genSnpmTable
    try
        func_genSnpmTable(results_struct, uncorrsigch, correctTFCEsigch, SnPMsigch, chanlocs, outputSname, select_mode);
        disp('Excel file generated successfully!');
    catch ME
        warning(ME.identifier, 'Could not generate Excel file: %s', ME.message);
    end

    % Generate text results
    results_text = {
        '=== SnPM Analysis Results ===',
        '',
        ['Data1 variable: ' params.data1_sheet],
        ['Data2 variable: ' params.data2_sheet],
        ['Analysis completed: ' datestr(now)],
        ['Comparison type: ' params.comparison],
        ['Data type: ' params.datatype],
        ['Tail: ' params.tail],
        ['Permutations: ' num2str(params.permutations)],
        '',
        '--- Data Summary ---',
        ['Data 1 size: ' mat2str(size(data_x))],
        ['Data 2 size: ' mat2str(size(data_y))],
        ['Data 1 mean: ' num2str(results_struct.data_summary.data1_mean, '%.4f')],
        ['Data 2 mean: ' num2str(results_struct.data_summary.data2_mean, '%.4f')],
        ['Data 1 std: ' num2str(results_struct.data_summary.data1_std, '%.4f')],
        ['Data 2 std: ' num2str(results_struct.data_summary.data2_std, '%.4f')]
        };

 

    if contains(params.comparison, 'correlation') && exist('subject_match_info', 'var')
        results_text{end+1} = '';
        results_text{end+1} = '--- Subject Matching (Correlation Analysis) ---';
        results_text{end+1} = ['Original subjects in Data1: ' num2str(subject_match_info.original_n1)];
        results_text{end+1} = ['Original subjects in Data2: ' num2str(subject_match_info.original_n2)];
        results_text{end+1} = ['Matched subjects used: ' num2str(subject_match_info.matched_n)];
        if subject_match_info.matched_n < max(subject_match_info.original_n1, subject_match_info.original_n2)
            results_text{end+1} = ['WARNING: Some subjects were excluded due to mismatch'];
        end
        
        % Store in results_struct too
        results_struct.subject_matching = subject_match_info;
        
        % Add covariate information if used
        if use_covariates && ~isempty(covariates)
            results_text{end+1} = '';
            results_text{end+1} = '--- Covariate Control ---';
            results_text{end+1} = ['Number of covariates: ' num2str(size(covariates, 2))];
            results_text{end+1} = ['Covariate names: ' strjoin(covariate_names, ', ')];
            if isfield(subject_match_info, 'covariate_matched_n')
                results_text{end+1} = ['Subjects with covariates: ' num2str(subject_match_info.covariate_matched_n)];
                if subject_match_info.covariate_excluded_n > 0
                    results_text{end+1} = ['Subjects excluded (no covariates): ' num2str(subject_match_info.covariate_excluded_n)];
                end
            end
            
            % Store in results_struct
            results_struct.covariates.names = covariate_names;
            results_struct.covariates.n_covariates = size(covariates, 2);
            results_struct.covariates.data = covariates;
        end
    end


    results_text{end+1} = '';
    results_text{end+1} = '--- Statistical Results ---';
    results_text{end+1} = ['Uncorrected significant channels: ' num2str(length(uncorrsigch))];
    results_text{end+1} = ['TFCE corrected significant channels: ' num2str(length(correctTFCEsigch))];
    results_text{end+1} = ['Significant clusters N = : ' num2str(length(sigclusters))];
    results_text{end+1} = ['SnPM significant channels N = : ' num2str(length(SnPMsigch))];
    results_text{end+1} = ['SnPM significant channels = : ' num2str(SnPMsigch)];
    results_text{end+1} = '';
    results_text{end+1} = ['Results saved to: ' outputSname '.mat'];
    results_text{end+1} = ['Excel file saved to: ' outputSname '.xlsx'];


    if exist('res', 'var') && ~isempty(res)
        results_text{end+1} = '';
        results_text{end+1} = '--- Global Statistics ---';
        results_text{end+1} = ['Global statistic: ' num2str(res.stat, '%.4f')];
        results_text{end+1} = ['Global p-value: ' num2str(res.pval, '%.4f')];
        
        % Save to results struct
        results_struct.global_stat = res.stat;
        results_struct.global_pval = res.pval;
    end

    %% Generate HTML Reports
    try
        fprintf('Generating HTML reports...\n');
        
        % Generate separate reports for TFCE and Cluster
        generateAnalysisReport(results_struct, params, base_filename, outputSname, uncorrsigch, correctTFCEsigch, SnPMsigch);
        
        % Add to results text
        results_text{end+1} = '';
        results_text{end+1} = ['Report saved to: ' outputSname '_report.html'];
        
        fprintf('HTML reports generated successfully!\n');
        
    catch ME
        warning(ME.identifier, 'Could not generate HTML reports: %s', ME.message);
    end

    % Save the GUI Result section as a CSV named to match the HTML report.
    try
        export_results_csv(results_text, outputSname, params);
        results_text{end+1} = ['Result section CSV saved to: ' outputSname '_report.csv'];
    catch ME
        warning(ME.identifier, 'Could not write result CSV: %s', ME.message);
    end
end




