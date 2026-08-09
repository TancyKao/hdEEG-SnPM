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

    % RETIRED circular comparison keys. These two named the pre-2026
    % Wheeler-Watson / Watson's U^2 engines (snpm_single_threshold_with_TFCE_circ
    % and snpm_cluster_analysis_circ), which have been deleted along with them.
    % Without this guard an old saved script carrying one of these keys falls
    % through EVERY route below into the legacy wide-format t-test flow and dies
    % deep inside a stat engine with a message that names neither the retired key
    % nor its replacement. Fail here instead, and say what to change.
    if ismember(params.comparison, {'circ_wheeler_watson_Test', 'circ_WatsonsU2Test'})
        error('core_snpm:retiredCircularComparison', ...
            ['params.comparison = ''%s'' is RETIRED: its engine has been deleted and it is ' ...
             'not a supported analysis. Replace it with one of the Tier-2 circular analyses:\n' ...
             '  ''circ_phase_group''      two independent groups, event-count-adjusted Hotelling ' ...
             'T^2 on the (cos,sin) embedding. PRIMARY replacement for a two-group phase ' ...
             'comparison (this is what circ_wheeler_watson_Test was used for).\n' ...
             '  ''circ_phase_group_u2''   same design, Watson''s U^2 omnibus. Replacement for ' ...
             'circ_WatsonsU2Test; note it CANNOT adjust for event count.\n' ...
             '  ''circ_corrAngLinear''    angle vs a linear measure, paired by subject.\n' ...
             'All three route to core_snpm_circ and REQUIRE params.circ_units and ' ...
             'params.circ_convention; the two phase-group analyses also require ' ...
             'params.count1_file / params.count2_file. See core_snpm_circ.m for the full ' ...
             'parameter list.'], params.comparison);
    end

    % Tier-2 circular (phase) analyses run on angles, which have no valid
    % datatype transform and no valid arithmetic channel mean, and they need a
    % per-channel event-count covariate. They route to their own pipeline; the
    % wide-format t-test / correlation flow below is untouched.
    if ismember(params.comparison, {'circ_phase_group', 'circ_phase_group_u2', ...
            'circ_corrAngLinear'})
        [results_struct, results_text] = core_snpm_circ(params);
        return
    end

    % Tier-1 guard: refuse an unadjusted analysis of the folded UNSIGNED
    % circular measure (snpm_circ_linearise). Runs before the GLM route so it
    % covers every preset, legacy and GLM alike.
    snpm_circ_tier1_guard(params);

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
        
        % NO SUBJECT IS REMOVED HERE, and that is NOT because the correlation
        % uses pairwise deletion -- it does not, and has not since 2026-08-02.
        % Missing data is handled later, by the COMPLETE-COLUMN rule: a CHANNEL
        % is tested only if it is usable in every analysed subject, the mask is
        % computed once before residualization and before the permutation loop,
        % and the excluded channels are reported (see the block headed
        % "CORRELATION ONLY -- evaluable channel mask" further down, the header
        % of dependencies/snpm_corr_columns.m, and the user-facing page
        % docs/explanation/missing-data-and-excluded-channels.md).
        %
        % Pairwise (per-channel) deletion was removed because it made the number
        % of complete pairs at a partially-missing channel a function of the
        % permutation, which mis-scales the null (measured family-wise error up
        % to 0.288 at nominal 0.05). The statistics printed below are DIAGNOSTIC
        % ONLY -- they describe the missingness, they do not describe how it is
        % handled.
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
        
        % Count complete pairs per channel. DIAGNOSTIC ONLY: this is the count a
        % pairwise-deletion analysis WOULD have used, and it is printed so the
        % analyst can see how much of the montage is at risk. Any channel whose
        % count here is below the matched-subject count is DROPPED ENTIRELY by
        % the complete-column rule; it is not analysed on a subset.
        valid_pairs_per_channel = zeros(1, total_channels);
        for ch = 1:total_channels
            valid_idx = ~isnan(data_matrix1(:, ch)) & ~isnan(data_matrix2(:, ch));
            valid_pairs_per_channel(ch) = sum(valid_idx);
        end
        fprintf(['  Complete subject pairs per channel (diagnostic, pre-mask): ' ...
            'min=%d, mean=%.1f, max=%d\n'], ...
            min(valid_pairs_per_channel), mean(valid_pairs_per_channel), max(valid_pairs_per_channel));
        fprintf(['  Channels short of %d complete pairs will be EXCLUDED, not analysed on a ' ...
            'subset: %d of %d\n'], matched_n, sum(valid_pairs_per_channel < matched_n), total_channels);

        % Provisional sample counts. They describe the SUBJECT MATCHING stage
        % only; the all-NaN row drop, the missing-covariate row drop and the
        % complete-column channel mask all happen later. Both fields are
        % refreshed from the analysed data before they are written to
        % results_struct.subject_matching (search for "Refresh the two fields").
        final_n = matched_n;
        nan_removed_n = 0;  % nothing dropped at THIS stage

        if matched_n < original_n1 || matched_n < original_n2
            missing1 = setdiff(subjects1, common_subjects);
            missing2 = setdiff(subjects2, common_subjects);
            
            warning_msg = sprintf(['Subject matching for correlation analysis:\n' ...
                                  'Original subjects in Data1: %d, Data2: %d\n' ...
                                  'Matched subjects: %d\n' ...
                                  'Missing from Data1: %s\n' ...
                                  'Missing from Data2: %s\n' ...
                                  'Note: remaining missing values are handled by the ' ...
                                  'complete-column rule -- a channel with any missing or ' ...
                                  'non-finite cell is excluded from the analysis entirely and ' ...
                                  'listed in the excluded-channel report, NOT analysed on the ' ...
                                  'subjects it happens to have.'], ...
                                  original_n1, original_n2, matched_n, ...
                                  strjoin(missing1, ', '), strjoin(missing2, ', '));
            warning(warning_msg);
        end

        
        if final_n < 3
            error('Correlation analysis requires at least 3 matching subjects. Found: %d', final_n);
        end
        
        fprintf(['Correlation analysis: %d matched subjects. Missing data is handled by the ' ...
            'COMPLETE-COLUMN rule (a channel must be usable in every analysed subject or it ' ...
            'is excluded and reported), not by pairwise deletion. See ' ...
            'docs/explanation/missing-data-and-excluded-channels.md\n'], matched_n);


        % Store subject matching info for results. final_n / nan_removed_n are
        % provisional here and are refreshed from the analysed data further down
        % (search for "Refresh the two fields") -- they used to be written out
        % unchanged, so the .mat claimed final_n = matched_n even for a run that
        % dropped rows or channels afterwards.
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
    is_source       = isfield(chsel, 'is_source') && chsel.is_source;
    fprintf('%s montage: %d/%d data columns matched to channel locations\n', ...
        params.channels, numel(insidegoodch), numel(channel_labels));

    % Source-space (2447 cortical voxels): enforce count + node-order invariants,
    % and reject signed current density on the group magnitude path. The stat
    % core is spatial-unit-agnostic; only ingestion/rendering differ.
    if is_source
        snpm_assert_source(chanlocs, neighbors, chsel.channel_mapping.channel_labels);
        % Non-negative magnitude guard: raw absolute band power must be >= 0.
        % Signed source current density (arbitrary per-subject sign) is invalid
        % for a group magnitude t-test -- cancellation destroys the effect.
        % data_x/data_y here are the RAW (untransformed) source inputs. For
        % 'absolute' a negative is a hard error. For logscale/normalize the
        % transform legitimately produces negatives, so we cannot inspect the
        % transformed values -- instead we check the RAW input and WARN (not
        % error, to avoid breaking legitimately-signed datatypes) if the raw
        % source magnitude is negative.
        neg_x = any(data_x(:) < 0);
        neg_y = ~onesample_vs_zero && any(data_y(:) < 0);
        if strcmp(params.datatype, 'absolute')
            if neg_x || neg_y
                error('snpm:source:signedMagnitude', ...
                    ['Source matrix contains negative values but non-negative ' ...
                     'band-power magnitudes are expected. Signed current density ' ...
                     'is invalid for a group magnitude t-test (arbitrary per-subject ' ...
                     'sign -> cancellation). Export voxel POWER (magnitude), or use a ' ...
                     'datatype that expects signed input.']);
            end
        elseif neg_x || neg_y
            warning('snpm:source:signedUnderTransform', ...
                ['Source input contains negative values with datatype=%s; if this ' ...
                 'is signed current density rather than a magnitude, the group test ' ...
                 'is invalid (arbitrary per-subject sign + cancellation).'], ...
                params.datatype);
        end
    end

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
    
    
    % ------------------------------------------------------------------
    % Drop subjects whose EVERY channel is NaN.
    %
    % Row alignment decides whether the drop is shared or independent:
    %   * pairedT / onesampleT / correlation -- the two matrices are ROW
    %     MATCHED (row i is the same subject in both), so a subject dropped
    %     from one must be dropped from the other: take the UNION.
    %   * unpairedT -- the two matrices are two INDEPENDENT groups with
    %     (possibly) different n. Their masks have different lengths, so they
    %     must be applied SEPARATELY. The old code OR'd them unconditionally,
    %     which threw MATLAB:sizeDimensionsMustMatch for every unpairedT with
    %     n1 ~= n2 whether or not any NaN was present -- i.e. every unequal
    %     between-group t-test in the toolbox was unreachable, including the
    %     spectral-folder unpairedT path where two cohort folders routinely
    %     differ in subject count.
    %
    % Dropping an all-NaN subject changes n1, which is fine and label
    % invariant: it is a property of the data, not of the group labels.
    % ------------------------------------------------------------------
    fprintf('Before NaN removal: data_x=%dx%d, data_y=%dx%d\n', size(data_x,1), size(data_x,2), size(data_y,1), size(data_y,2));

    rows_are_matched = ~strcmp(params.comparison, 'unpairedT');
    all_nan_x = all(isnan(data_x), 2);
    all_nan_y = all(isnan(data_y), 2);
    keep_rows_x = true(size(data_x,1),1);   % used below to keep covariates aligned

    if rows_are_matched
        if numel(all_nan_x) ~= numel(all_nan_y)
            error('snpm:rowMismatch', ...
                ['%s is a row-matched design but Data 1 has %d rows and Data 2 has %d. ' ...
                 'Subjects must be matched before analysis.'], ...
                params.comparison, size(data_x,1), size(data_y,1));
        end
        subjects_to_remove = all_nan_x | all_nan_y;
        if any(subjects_to_remove)
            fprintf('Removing %d subjects with all NaN channels (row-matched: dropped from both)\n', ...
                sum(subjects_to_remove));
            data_x(subjects_to_remove, :) = [];
            data_y(subjects_to_remove, :) = [];
            keep_rows_x = ~subjects_to_remove;
        end
    else
        if any(all_nan_x)
            fprintf('Removing %d Group A subjects with all NaN channels\n', sum(all_nan_x));
            data_x(all_nan_x, :) = [];
        end
        if any(all_nan_y)
            fprintf('Removing %d Group B subjects with all NaN channels\n', sum(all_nan_y));
            data_y(all_nan_y, :) = [];
        end
        keep_rows_x = ~all_nan_x;
    end

    fprintf('After NaN removal: data_x=%dx%d, data_y=%dx%d\n', size(data_x,1), size(data_x,2), size(data_y,1), size(data_y,2));

    if size(data_x,1) < 2 || size(data_y,1) < 2
        error('snpm:tooFewSubjects', ...
            'After dropping all-NaN subjects only %d and %d rows remain; not enough to test.', ...
            size(data_x,1), size(data_y,1));
    end

    % Keep the covariate matrix on the same rows as the data. (The old code
    % dropped rows from data_x/data_y only, silently misaligning covariates
    % against subjects in the partial-correlation residualization.)
    if ~isempty(covariates) && size(covariates,1) == numel(keep_rows_x)
        covariates = covariates(keep_rows_x, :);
    end

    % Calculate degrees of freedom (on the rows that actually remain)
    g1subj_n = size(data_x,1);
    g2subj_n = size(data_y,1);


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
    
    % Verify dimensions match for the ROW-MATCHED designs only. This used to be
    % written as contains(comparison,'pairedT'), which is also true for
    % 'unpairedT' -- so an unequal-group between-subjects t-test was rejected
    % here even after surviving the row-drop above. Two independent groups are
    % supposed to have different n.
    if rows_are_matched
        if size(data_x, 1) ~= size(data_y, 1)
            error('Data dimension mismatch after NaN removal: data_x has %d rows, data_y has %d rows', ...
                size(data_x, 1), size(data_y, 1));
        end
    end

    % ------------------------------------------------------------------
    % CORRELATION ONLY -- evaluable channel mask (permutation-invariant row set)
    %
    % A channel is evaluable if and only if it is COMPLETE across all analysed
    % rows. The mask is computed ONCE, here, before any residualization and
    % before the permutation loop, and both the observed map and every permuted
    % map are defined on exactly that channel set.
    %
    % WHY. The correlation null permutes the rows of data_x against a fixed
    % data_y. With pairwise (per-channel) deletion the number of complete pairs
    % at a partially-missing channel is a FUNCTION OF THE PERMUTATION, so the
    % observed r and the permuted r are computed at different sample sizes; raw
    % r has null spread ~1/sqrt(n-1), so the null is mis-scaled. Measured at a
    % channel with 25% missingness in a 30-subject design: observed n = 22 while
    % the permuted n ranged [14, 20] -- the observed value fell outside the
    % entire permuted range. Family-wise error at nominal 0.05, 2000 replicates:
    % 0.072 (25% disjoint missingness) and 0.288 (40%) for TFCE; with covariates
    % and the old zeros-imputation, 0.168 TFCE / 0.076 cluster. The rule below
    % removes both mechanisms at once, at the cost of dropping a whole channel
    % for one missing cell (the same trade the GLM tier already makes, see
    % dependencies/snpm_glm_fl_context.m).
    %
    % THE T-TESTS ARE DELIBERATELY NOT GIVEN THE COMPLETE-COLUMN RULE. What
    % separates them from correlation is not the effective n, it is the ANALYSED
    % SUBJECT SET:
    %
    %   1. Missingness travels with the row. A sign-flip swaps a subject's two
    %      conditions and a relabelling moves a whole subject between groups, so
    %      the SET of subjects that can be evaluated at a channel is the same for
    %      the observed labelling and for every permuted one. In the correlation
    %      null the rows of data_x are permuted against a FIXED data_y, which
    %      re-pairs x_i with y_j, so the set of complete pairs genuinely changes
    %      from permutation to permutation. That is the difference, and it is
    %      why one path needs the rule and the other does not.
    %
    %   2. What DOES vary under relabelling is the n1/n2 split of that fixed set,
    %      and with it the scale of the statistic: at 40% missingness the pooled
    %      standard error sqrt(1/n1 + 1/n2) moves by up to 34% across the
    %      reachable splits. (The pooled df n1+n2-2 does NOT move -- it is
    %      invariant by construction, so quoting its range, as an earlier version
    %      of this comment did, measured nothing.) That variation does not
    %      threaten validity: a permutation test is exact for ANY statistic, so
    %      long as the same function of the data is computed for the observed and
    %      the permuted labellings. Nichols & Holmes 2001 (Ref/NicholsHolmes.pdf,
    %      p. 7) say exactly this about a statistic whose variance differs across
    %      the map -- "The test will still be valid, but will be less sensitive".
    %      Measured family-wise error at nominal 0.05 across 20-40% missingness,
    %      including missingness confounded with group: 0.045-0.055.
    %
    %   3. What it costs is SENSITIVITY, and at the extreme the statistic stops
    %      being defined at all -- which is what the degeneracy guard below is
    %      for. That guard is a POWER guard, not a validity guard.
    % ------------------------------------------------------------------
    is_correlation = ismember(params.comparison, {'correlationP', 'correlationS'});
    is_ttest       = ismember(params.comparison, {'pairedT', 'unpairedT', 'onesampleT'});
    excluded_channels = [];
    exclusion_spatial = [];
    evaluable = [];                 % set by whichever channel rule applies below
    % n_available is assigned unconditionally a few lines down (every comparison
    % gets one), so it is deliberately NOT pre-initialised to [] here.

    % The comparison key the STAT ENGINES are given. It equals params.comparison
    % everywhere except partial Spearman, where the rank transform is applied to
    % the data before residualization and the engines then take a linear
    % (Pearson) correlation of the rank residuals -- see the residualization
    % block below. params.comparison itself is never rewritten: it is what the
    % filename, the topoplot and the report call the analysis.
    stat_comparison = params.comparison;

    if is_correlation
        % A subject with a missing COVARIATE cannot enter the partial
        % correlation at ANY channel (the design matrix row is undefined), so
        % it is dropped here rather than silently turning the whole map into
        % NaN, which is what the old code did.
        if use_covariates && ~isempty(covariates)
            if size(covariates,1) ~= size(data_x,1)
                error('snpm:covariateRowMismatch', ...
                    ['Covariate matrix has %d rows but %d subjects are analysed. Residualizing ' ...
                     'with misaligned rows would silently regress out the wrong subject.'], ...
                    size(covariates,1), size(data_x,1));
            end
            bad_cov = any(isnan(covariates), 2);
            if any(bad_cov)
                warning('snpm:covariateNaNRows', ...
                    ['%d of %d subjects have a missing covariate value and cannot enter the ' ...
                     'partial correlation at any channel; they are dropped from the analysis.'], ...
                    sum(bad_cov), numel(bad_cov));
                data_x(bad_cov, :) = [];
                data_y(bad_cov, :) = [];
                covariates(bad_cov, :) = [];
                if size(data_x,1) < 3
                    error('snpm:tooFewSubjects', ...
                        'Only %d subjects remain after dropping missing covariates; correlation needs >= 3.', ...
                        size(data_x,1));
                end
            end
        end

    end

    % ------------------------------------------------------------------
    % PER-CHANNEL AVAILABILITY -- computed for EVERY comparison, not just
    % correlation. It is what both channel rules below are built from, and it is
    % the data-quality signal a reader needs whichever analysis was run.
    %
    % USABLE CELL = FINITE, not merely non-NaN.
    %
    % isfinite, not ~isnan. With datatype 'logscale' a zero-power cell becomes
    % log10(0) = -Inf, and 'normalize' turns a zero-variance subject row into
    % Inf/NaN. -Inf passes ~isnan, so under an ~isnan test such a channel was
    % declared evaluable and never appeared in the exclusion report -- while
    % snpm_corr_columns returns NaN for a non-finite column regardless (its
    % r(~isfinite(r)) = NaN backstop). The channel therefore left the analysis
    % silently, which is exactly the failure mode the complete-column rule was
    % introduced to end. The GLM and LMM tiers already gate on isfinite
    % (dependencies/snpm_glm_fl_context.m, core_snpm_lmm.m); this makes the
    % mask, the exclusion report and the statistic agree across all three tiers.
    %
    % The count is PRE-mask, so it says WHY a channel was excluded and how close
    % the survivors came to being excluded. The unit it counts differs with the
    % design and that is deliberate: matched PAIRS for the row-matched analyses
    % (paired/one-sample/correlation), and total available SUBJECTS across both
    % groups for the unpaired test, where the two arms have different rows and a
    % "pair" does not exist.
    % ------------------------------------------------------------------
    if strcmp(params.comparison, 'unpairedT')
        n_available_x = sum(isfinite(data_x), 1);
        n_available_y = sum(isfinite(data_y), 1);
        n_available   = n_available_x + n_available_y;
    else
        usable = isfinite(data_x) & isfinite(data_y);
        n_available = sum(usable, 1);
    end

    if is_correlation
        evaluable = all(usable, 1);
        excluded_channels = snpm_excluded_channel_info(evaluable, chanlocs, params.comparison);
        excluded_channels.reason = ['missing or non-finite data (NaN, +/-Inf): channel not ' ...
            'usable in every analysed subject'];

        % A NON-FINITE, NON-NaN cell gets its own warning. It is not missing
        % data and the fix is different: -Inf here is log10(0) under
        % datatype 'logscale', i.e. a channel whose power really was zero and
        % which the transform sent to -Inf; +/-Inf or NaN under 'normalize' is a
        % subject whose channels were all identical (zero across-channel SD).
        % Folding it into the generic "missing data" message would send the
        % analyst looking for an absent recording that is not the problem.
        nonfinite_notnan = any((~isfinite(data_x) & ~isnan(data_x)) | ...
                               (~isfinite(data_y) & ~isnan(data_y)), 1);
        if any(nonfinite_notnan)
            nf_idx = find(nonfinite_notnan);
            nf_lab = cell(1, numel(nf_idx));
            for k = 1:numel(nf_idx)
                lab = '';
                if nf_idx(k) <= numel(chanlocs) && isfield(chanlocs, 'labels')
                    lab = chanlocs(nf_idx(k)).labels;
                end
                if isempty(lab), lab = sprintf('#%d', nf_idx(k)); end
                nf_lab{k} = lab;
            end
            warning('snpm:nonFiniteAfterTransform', ...
                ['%d channel(s) hold a non-finite value that is NOT missing data, after the ' ...
                 '''%s'' transform: %s. With ''logscale'' this is log10(0) -- a zero-power ' ...
                 'cell, not an absent one; with ''normalize'' it is a subject with zero ' ...
                 'across-channel variance. They are excluded from the analysis together with ' ...
                 'the genuinely missing channels. Fix the transform (floor the zeros, or use ' ...
                 '''absolute'') rather than treating these as missing recordings.'], ...
                numel(nf_idx), params.datatype, strjoin(nf_lab, ', '));
        end

        % WHERE the exclusions fell, not just how many. A survival count is a
        % power question; a survival PATTERN is a validity question -- losing
        % 12 scattered channels and losing 12 contiguous frontal channels have
        % the same count and completely different consequences for a frontal
        % hypothesis. Computed from the montage coordinates carried on chanlocs.
        exclusion_spatial = snpm_exclusion_spatial_profile(evaluable, chanlocs);
        if excluded_channels.n > 0 && exclusion_spatial.available
            fprintf('%s\n', exclusion_spatial.summary);
            if exclusion_spatial.lopsided
                warning('snpm:exclusionsLopsided', '%s', exclusion_spatial.summary);
            end
        end

        if ~any(evaluable)
            error('snpm:noEvaluableChannels', ...
                ['No channel is complete across all %d analysed subjects, so no channel can be ' ...
                 'tested with a permutation-invariant row set. Drop the subjects driving the ' ...
                 'missingness, or repair the data.'], size(data_x,1));
        end

        % Apply the mask BEFORE residualization so the covariate fit and the
        % permutation both see one channel set.
        data_x(:, ~evaluable) = NaN;
        data_y(:, ~evaluable) = NaN;

        if excluded_channels.n / numel(evaluable) > 0.10
            warning('snpm:manyChannelsExcluded', ...
                ['%d of %d channels (%.1f%%) were dropped for missing data. Above roughly 10%% ' ...
                 'it is usually better to drop the few subjects driving the missingness than to ' ...
                 'lose that much of the montage -- inspect the per-subject NaN counts printed above.'], ...
                excluded_channels.n, numel(evaluable), 100*excluded_channels.n/numel(evaluable));
        end

    elseif is_ttest
        % ------------------------------------------------------------------
        % DEGENERACY GUARD (t-tests). NOT the complete-column rule.
        %
        % The t-test maps stay otherwise unmasked -- a channel missing in a few
        % subjects is still tested, on the subjects it has, and that is valid
        % (see the three-point note above the correlation mask). What is masked
        % here is the much smaller set of channels where the statistic is
        % UNDEFINED OR DEGENERATE FOR SOME LABELLING IN THE PERMUTATION SET.
        % That is a power guard, not a validity guard, and the rule is derived
        % from the permutation set rather than tuned:
        %
        %   pairedT / onesampleT -- a sign flip swaps a subject's two conditions
        %       and never changes how many pairs a channel has, so the reachable
        %       n is the constant m. The statistic needs a variance, hence
        %           m >= 2 usable pairs.
        %       At m = 1 MATLAB's ttest returns df = 0 and t = -Inf, and
        %       ClusterEnhancement then builds hrange = dh:dh:max(abs(map))+dh
        %       from an infinite maximum and dies with MATLAB:pmaxsize
        %       ("Requested array exceeds the maximum possible variable size") --
        %       an opaque crash with no hint of which channel caused it.
        %
        %   unpairedT -- a relabelling moves whole subjects between the groups,
        %       so with m available subjects out of n1 + n2 the group-1 count
        %       reachable under some relabelling runs over
        %           max(0, m - n2) ... min(m, n1).
        %       Both groups hold >= 2 for EVERY relabelling iff
        %           m >= max(n1, n2) + 2.
        %       Below that some labelling puts one group at a single subject,
        %       where MATLAB's var of a scalar is 0 (not NaN), so ttest2 returns
        %       a FINITE t computed from a pooled SD that ignores that group and
        %       |t| explodes: measured max |t| = 7.4 over 20000 draws with
        %       n1 = 1, n2 = 10 on pure noise. That value enters the
        %       max-statistic null and sets the corrected threshold for the whole
        %       montage. Measured cost, 8x8 lattice, 800 reps: a planted,
        %       fully-observed 9-channel cluster is TFCE-detected 93.4% of the
        %       time with no degenerate channel present and 16.4% when four
        %       channels are available in only 3 of 20 subjects; masking those
        %       four restores it to 94.1%.
        %
        % On complete data both rules are no-ops (m = n >= 2 paired;
        % m = n1 + n2 >= max(n1,n2) + 2 whenever both groups have >= 2 rows,
        % which the earlier row guard already enforces).
        % ------------------------------------------------------------------
        n1 = size(data_x, 1);
        n2 = size(data_y, 1);
        if strcmp(params.comparison, 'unpairedT')
            min_units = max(n1, n2) + 2;
            unit_word = 'usable subjects (both groups pooled)';
            rule_word = sprintf(['degenerate under relabelling: fewer than %d of the %d ' ...
                'subjects are usable, so some group relabelling leaves one group with ' ...
                '< 2 subjects and ttest2 returns a finite t from a zero-variance group'], ...
                min_units, n1 + n2);
            why_word = sprintf(['They have fewer than %d usable subjects out of %d, so SOME ' ...
                'group relabelling in the permutation set puts one group at a single subject, ' ...
                'where MATLAB''s var of a scalar is 0 (not NaN) and ttest2 returns a finite t ' ...
                'from a pooled SD that ignores that group (measured max |t| = 7.4 on pure ' ...
                'noise with n1 = 1, n2 = 10). That value would set the max-statistic threshold ' ...
                'for the whole montage. This is the DEGENERACY guard, not the complete-column ' ...
                'rule: every other channel is still tested on the subjects it has.'], ...
                min_units, n1 + n2);
        else
            min_units = 2;
            unit_word = 'usable pairs';
            rule_word = ['degenerate under relabelling: fewer than 2 usable pairs, so the ' ...
                'paired t has df = 0 and returns +/-Inf for every sign-flip labelling'];
            why_word = ['They have fewer than 2 usable pairs, so the paired t has df = 0 and ' ...
                'returns +/-Inf under EVERY sign-flip labelling (sign flipping never changes ' ...
                'how many pairs a channel has). An infinite value also breaks the TFCE ' ...
                'threshold ladder outright (MATLAB:pmaxsize). This is the DEGENERACY guard, ' ...
                'not the complete-column rule: every other channel is still tested on the ' ...
                'pairs it has.'];
        end

        evaluable = n_available >= min_units;
        excluded_channels = snpm_excluded_channel_info(evaluable, chanlocs, ...
            params.comparison, rule_word, why_word);

        exclusion_spatial = snpm_exclusion_spatial_profile(evaluable, chanlocs);
        if excluded_channels.n > 0 && exclusion_spatial.available
            fprintf('%s\n', exclusion_spatial.summary);
            if exclusion_spatial.lopsided
                warning('snpm:exclusionsLopsided', '%s', exclusion_spatial.summary);
            end
        end

        if ~any(evaluable)
            error('snpm:noEvaluableChannels', ...
                ['No channel has the %d %s a %s needs under every labelling in the ' ...
                 'permutation set, so no channel can be tested. Drop the subjects driving ' ...
                 'the missingness, or repair the data.'], ...
                min_units, unit_word, params.comparison);
        end

        % Mask before the global test and before either permutation engine, so
        % the observed map, the null and the omnibus all see one channel set.
        data_x(:, ~evaluable) = NaN;
        data_y(:, ~evaluable) = NaN;
    end

    % TWO effective-N columns, deliberately not one.
    %
    %   n_analysed  post-mask. Under the whole-column rule this can only take
    %               two values -- the analysed subject count on every retained
    %               channel, and 0 on every excluded one -- so it is really a
    %               constant plus an exclusion flag. Reported per channel for
    %               the topoplot/table join, and the constant is stated once in
    %               the text header so a reader is not invited to look for
    %               variation that cannot exist.
    %   n_available pre-mask, genuinely per channel (old pairwise semantics).
    %               This is the diagnostic: it names the channels that were one
    %               subject short and the ones that were half empty.
    %
    % A single column called "per-channel n" conflated the two and looked as if
    % it varied meaningfully across retained channels when it cannot.
    %
    % FOR THE T-TESTS both columns exist too, but n_analysed genuinely VARIES
    % across retained channels, because only the degenerate channels are masked
    % and a channel missing in two subjects is still tested on the rest. There
    % is therefore no single constant to quote, so n_analysed_constant stays
    % empty for them and the "Effective N" text block (which states that
    % constant) does not fire.
    n_analysed = [];
    n_analysed_constant = [];
    if is_correlation
        % Same usability test as the mask above (isfinite, not ~isnan), so the
        % reported per-channel n cannot disagree with the channel set that was
        % actually analysed. Post-mask this is the analysed subject count on
        % every retained channel and 0 on every excluded one.
        n_analysed = sum(isfinite(data_x) & isfinite(data_y), 1);
        if any(n_analysed > 0)
            n_analysed_constant = max(n_analysed);
        end
    elseif is_ttest
        if strcmp(params.comparison, 'unpairedT')
            n_analysed = sum(isfinite(data_x), 1) + sum(isfinite(data_y), 1);
        else
            n_analysed = sum(isfinite(data_x) & isfinite(data_y), 1);
        end
    end

    % Apply partial correlation if covariates are provided for correlation analyses
    if contains(params.comparison, 'correlation') && use_covariates && ~isempty(covariates)
        fprintf('Applying partial correlation with %d covariate(s)\n', size(covariates, 2));

        % PARTIAL SPEARMAN: RANK FIRST, THEN RESIDUALIZE.
        %
        % This path used to residualize the RAW data and let the engine rank the
        % residuals afterwards. That is robust in its second step and NOT robust
        % in its first: a least-squares fit on the raw scale lets one outlier
        % distort the slope, and a distorted slope contaminates EVERY residual,
        % which ranking afterwards cannot undo. Ranking first bounds each
        % subject's influence before anything is fitted. It also matters for the
        % covariates this toolbox actually gets -- age, and Helfrich's sixteen
        % sleep-architecture variables -- whose relationships with the outcome
        % are monotone but frequently nonlinear, exactly the part a linear
        % residualization on the raw scale leaves in the residual.
        %
        % It also removes an internal inconsistency: zero-order correlationS
        % matched corr(...,'Type','Spearman') to 5.6e-17, and then the same named
        % analysis silently switched estimator as soon as a covariate was added.
        % Rank-then-residualize IS what partialcorr(x,y,Z,'Type','Spearman')
        % computes, so the map is now checkable against a builtin to machine
        % precision. Validated against partialcorr (map == partialcorr to
        % 1e-12) and T10 (the same equality end-to-end through this function),
        % and by the estimator note printed in T5.
        %
        % The permutation machinery is untouched: ranks are a fixed relabelling
        % of 1..n, so permuting rows of the rank residuals is the same operation
        % it always was, and the df used by the partial-correlation global p is
        % unchanged.
        %
        % Because the statistic on rank residuals is a LINEAR correlation (that
        % is the definition of the rank partial correlation -- you do not rank
        % twice), the engines are handed 'correlationP' for this run. The
        % analysis is still reported as correlationS.
        if strcmp(params.comparison, 'correlationS')
            data_x     = tiedrank(data_x);
            data_y     = tiedrank(data_y);
            covariates = tiedrank(covariates);
            stat_comparison = 'correlationP';
            fprintf(['Partial Spearman: ranking data and covariates before residualization, ' ...
                'then scoring a LINEAR correlation of the rank residuals in every engine ' ...
                '(the map equals partialcorr(x,y,Z,''Type'',''Spearman'')).\n']);
        end

        % Residualize both data matrices against [1, covariates], per channel.
        % Initialised to NaN, NOT zeros: with zeros a cell that was missing
        % survived residualization as a finite 0 -- exactly the residual mean --
        % which (a) left Spearman rho distorted by more than 0.05 in 29.2% of
        % channels over 2000 draws at 20% missingness, (b) reported a parametric
        % p on an inflated n, and (c) ERASED the missingness pattern, freezing
        % the pairwise-deletion denominators across permutations while the
        % numerator still lost mismatched pairs. Because the evaluable mask is
        % applied above, every retained channel is complete here, so the
        % sum(valid_idx) branch below can only fire on a channel the mask has
        % already dropped.
        data_x_residual = NaN(size(data_x));
        data_y_residual = NaN(size(data_y));

        for ch = 1:size(data_x, 2)
            % Get valid indices for this channel. isfinite, not ~isnan, to match
            % the evaluable mask above: on a retained channel every cell is
            % finite so the two agree, and on an excluded channel both are
            % all-false. Using ~isnan here would let a non-finite cell into the
            % design matrix if the mask above were ever bypassed.
            valid_idx = isfinite(data_x(:, ch)) & isfinite(data_y(:, ch));
            
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

    % Global statistical analysis (channel-averaged). Computed AFTER any
    % covariate residualization so the correlation global p reflects the partial
    % correlation; ncov shrinks its df to n-2-ncov (descriptive p only).
    ncov = 0;
    if contains(params.comparison, 'correlation') && use_covariates && ~isempty(covariates)
        ncov = size(covariates, 2);
    end

    % STAT_COMPARISON, NOT PARAMS.COMPARISON, FROM HERE DOWN.
    % The two differ on exactly one route: partial Spearman, where the block
    % above has already rank-transformed data_x, data_y and the covariates and
    % residualized the RANKS, so what remains for the engines is a LINEAR
    % correlation of rank residuals -- the definition of the rank partial
    % correlation, and what partialcorr(...,'Type','Spearman') returns.
    %
    % Handing the engines params.comparison ('correlationS') here made them
    % tiedrank the rank residuals a SECOND time, so the map was Spearman of the
    % residuals of ranks: a monotone functional nobody chose, equal to neither
    % the estimator named in the console message nor the pre-2026-08 behaviour
    % (measured n = 40, one covariate, seed 7: 0.1512 double-ranked vs 0.1808
    % partialcorr-Spearman vs 0.2077 residualize-raw-then-Spearman). It was not
    % an inference bug -- the same functional reached the observed and the
    % permuted maps, so the test stayed exact -- but the reported number was not
    % the one the run claimed to compute.
    %
    % All THREE consumers take it, deliberately: the omnibus, the TFCE engine
    % and the cluster engine must score the same functional or the global p and
    % the map describe different estimators. params.comparison is untouched and
    % still names the analysis in the filename, the topoplot and the report.
    res = global_stat_test(data_x, data_y, alpha, stat_comparison, params.tail, ncov);

    % Main permutation analysis. Circular comparisons never reach this point --
    % they are routed to core_snpm_circ by the early guard above, and the two
    % legacy _circ engines have been deleted.
    [T,p,~,~] = snpm_single_threshold_with_TFCE(data_x,data_y,neighbors,E,H,alpha,stat_comparison,params.tail,params.permutations);

    % Cluster analysis
    [Clusters] = snpm_cluster_analysis(data_x,data_y,threshold,neighbors,alpha,stat_comparison,params.tail,params.permutations);

    % Find significant channels/clusters. One rule for every analysis: the
    % uncorrected set is the channels whose uncorrected p is at or below alpha.
    % (The retired circular branch compared the raw statistic against a
    % cluster-forming threshold, which is not a p-value at all.)
    uncorrsigch = find(p.real<=alpha);
    correctTFCEsigch = find(p.correctedTFCE <=alpha);

    sigclusters = find([Clusters.p] <=alpha);
    SnPMsigch = [Clusters(sigclusters).channels];
    
    

    %% Plot figures
    % Extract basenames from file paths
    [~, data1_basename] = fileparts(params.data1_file);
    [~, data2_basename] = fileparts(params.data2_file);
    
    data1_basename = sprintf('%s_%s',data1_basename, params.data1_sheet);
    data2_basename = sprintf('%s_%s',data2_basename, params.data2_sheet);
    if onesample_vs_zero, data2_basename = 'zero'; end   % single condition vs 0

    % Create filename with basenames and sheet names
    base_filename = sprintf('%s VS %s', ...
        data1_basename, data2_basename ); % use VS to separate two lines in figure title
    
    % Set savepath for plots
    savepath = params.output_path;
    
    % plot individual figure & between two groups/conditions
    % Source-space systems have NO scalp layout: bypass every EEGLAB topoplot
    % (a 2-D scalp disc is meaningless for cortical voxels). The statistically
    % sufficient deliverable -- the significant-voxel list with coordinates and
    % cluster membership -- is emitted below via write_source_voxel_table.
    if ~is_source

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
    end

    % plot t-map or corr-map
    try
        TopoplotSignificant_single(T.real_T,uncorrsigch,correctTFCEsigch,SnPMsigch,chanlocs,...
            insidegoodch,params.comparison, savepath,base_filename, select_mode);
    catch ME
        warning(ME.identifier, 'Could not generate significance topology plots: %s', ME.message);
    end

    else
        % Source-space deliverable: significant-voxel list with coordinates and
        % cluster membership (topoplot skipped). Written as its own CSV next to
        % the .mat/.xlsx so a downstream 3-D brain viewer (deferred task) can
        % consume it directly.
        fprintf('Source-space system: scalp topoplot bypassed; emitting voxel list.\n');
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
    % Effective N -- the two columns described above the mask, wired to the
    % output contract here. They are absent entirely from the GLM/LMM structs and
    % every consumer is isfield/isempty guarded, so nothing else is affected.
    %
    % The t-tests now populate .per_channel_n and .n_available too (they used to
    % be [] for them). For a t-test .per_channel_n is NOT constant across
    % retained channels -- only the degenerate channels are masked -- so
    % .n_analysed_constant stays [] there, which is the flag the text header uses
    % to decide whether a single number is an honest summary.
    %
    %   .per_channel_n        POST-mask: the subject count the statistic was
    %                         actually computed at. Constant on every retained
    %                         channel, 0 on every excluded one. The name is kept
    %                         (rather than renamed to n_analysed) because it is
    %                         the published field: func_genSnpmTable's
    %                         'effectiveN' sheet, the effective-N
    %                         excluded==0 / retained==n assertion, and
    %                         core_snpm_circ all speak it.
    %   .n_available          PRE-mask matched pairs per channel: the genuinely
    %                         per-channel diagnostic that separates a channel
    %                         that was one subject short from one that was half
    %                         empty.
    %   .n_analysed_constant  that retained-channel constant, stated once so no
    %                         reader hunts for across-channel variation the
    %                         complete-column rule makes impossible.
    results_struct.per_channel_n       = n_analysed;
    results_struct.n_available         = n_available;
    results_struct.n_analysed_constant = n_analysed_constant;
    if ~isempty(excluded_channels)
        % Channels dropped for missing data (correlation only). Feeds the
        % 'excludedChannels' Excel sheet and the HTML banner, both of which
        % already exist for the GLM tier -- no new reporting surface.
        results_struct.excluded_channels = excluded_channels;
    end
    if ~isempty(exclusion_spatial)
        % WHERE those channels were, so the .mat carries the validity signal and
        % not just the count (see snpm_exclusion_spatial_profile).
        results_struct.exclusion_spatial = exclusion_spatial;
    end
    results_struct.hide_condition_b = onesample_vs_zero;  % one condition vs 0 -> no B card
    results_struct.is_source = is_source;           % source-space run (topoplot bypassed)
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

    % Source-space deliverable: significant-voxel list with coordinates +
    % cluster membership (replaces the scalp topoplot, which was bypassed above).
    if is_source
        try
            write_source_voxel_table(results_struct, chanlocs, outputSname);
        catch ME
            warning(ME.identifier, 'Could not write source voxel table: %s', ME.message);
        end
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
        results_text{end+1} = ['Subjects matched across the two files: ' num2str(subject_match_info.matched_n)];
        if subject_match_info.matched_n < max(subject_match_info.original_n1, subject_match_info.original_n2)
            results_text{end+1} = ['WARNING: Some subjects were excluded due to mismatch'];
        end

        % Refresh the two fields that are supposed to describe the FINAL
        % analysed sample. They are set at subject-matching time, BEFORE the
        % all-NaN row drop, the missing-covariate row drop and the
        % complete-column channel mask, and were previously written out
        % unchanged -- so the .mat and the results text could report an n larger
        % than the analysis actually used. Under the complete-column rule every
        % retained channel is analysed on the same rows, so a single number is
        % the honest summary.
        if ~isempty(n_analysed_constant)
            analysed_n = n_analysed_constant;
        else
            analysed_n = size(data_x, 1);
        end
        subject_match_info.final_n       = analysed_n;
        subject_match_info.nan_removed_n = subject_match_info.matched_n - analysed_n;
        if ~isempty(excluded_channels) && isstruct(excluded_channels)
            subject_match_info.n_channels_excluded = excluded_channels.n;
        end

        results_text{end+1} = ['Subjects entering every retained channel: ' num2str(analysed_n) ...
            ' (complete-column rule; identical row set for the observed map and every permutation)'];
        if subject_match_info.nan_removed_n > 0
            results_text{end+1} = ['Matched subjects dropped before analysis: ' ...
                num2str(subject_match_info.nan_removed_n) ' (all-NaN rows and/or missing covariates)'];
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


    if ~isempty(excluded_channels)
        results_text{end+1} = '';
        results_text{end+1} = '--- Channel Exclusions ---';
        if excluded_channels.n == 0
            results_text{end+1} = sprintf('Channels excluded: none (%d/%d evaluable)', ...
                excluded_channels.n_channels, excluded_channels.n_channels);
        else
            results_text{end+1} = sprintf('Channels excluded (%s): %d of %d - %s', ...
                excluded_channels.context, excluded_channels.n, excluded_channels.n_channels, ...
                strjoin(excluded_channels.labels, ', '));
            % The REASON is different in the two tiers and the distinction is the
            % whole point: correlation drops any incomplete channel (the complete
            % -column rule, a validity fix), a t-test drops only the channels that
            % are degenerate under some labelling (a power fix). Print it so a
            % reader is not left to assume the t-test map was masked for
            % missingness, which it was not.
            results_text{end+1} = ['Reason: ' excluded_channels.reason];
        end
        % WHERE the exclusions fell. A count is a power statement; the pattern is
        % a validity statement -- see snpm_exclusion_spatial_profile.
        if ~isempty(exclusion_spatial) && exclusion_spatial.available
            results_text{end+1} = exclusion_spatial.summary;
        end
        % For the t-tests the retained channels are NOT all analysed at the same
        % n, so state the range instead of a constant (the correlation tier
        % states its constant in the Effective N block below).
        if is_ttest && ~isempty(n_analysed)
            kept = n_analysed(evaluable);
            if strcmp(params.comparison, 'unpairedT')
                unit_txt = sprintf('subjects (of %d)', size(data_x,1) + size(data_y,1));
            else
                unit_txt = sprintf('pairs (of %d)', size(data_x,1));
            end
            results_text{end+1} = sprintf(['Retained channels are analysed on %d-%d %s ' ...
                '- the t-test map is deliberately NOT masked for ordinary missingness ' ...
                '(per-channel detail in the effectiveN sheet).'], ...
                min(kept), max(kept), unit_txt);
        end
    end

    % Effective N, stated once. Under the complete-column rule the analysed n is
    % the SAME on every retained channel, so quoting it per channel invites a
    % reader to look for variation that cannot exist; the per-channel detail that
    % does vary (what each channel HAD before the rule was applied) is the
    % pre-mask count, and it goes in the workbook's effectiveN sheet.
    if ~isempty(n_analysed_constant)
        results_text{end+1} = '';
        results_text{end+1} = '--- Effective N ---';
        results_text{end+1} = sprintf(['Subjects entering every retained channel: %d ' ...
            '(constant by construction - a channel is analysed only if it is complete).'], ...
            n_analysed_constant);
        if ~isempty(n_available)
            dropped = (n_analysed == 0);
            if any(dropped)
                results_text{end+1} = sprintf(['Excluded channels had %d-%d of %d subjects ' ...
                    'available before exclusion (per-channel detail in the effectiveN sheet).'], ...
                    min(n_available(dropped)), max(n_available(dropped)), n_analysed_constant);
            end
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

        % STATE K. The global is NOT a whole-head test in general: it is the mean
        % over the K channels that are usable in every analysed unit of BOTH arms
        % (see global_stat_test's common-channel rule). Averaging each unit over
        % its own available channels would make the arms measure different things
        % -- with a front-to-back power gradient that alone rejected 100% of the
        % time on null paired data with a different region missing per condition.
        % When K is much smaller than the montage this is a materially different
        % scientific claim from "whole head", so the count is reported, not
        % implied, and a large loss is called out.
        if isfield(res, 'n_channels') && isfield(res, 'n_channels_total')
            results_text{end+1} = sprintf(['Global scope: mean over %d of %d channels ' ...
                '(those complete in every analysed unit of both arms).'], ...
                res.n_channels, res.n_channels_total);
            if res.n_channels == 0
                results_text{end+1} = ['WARNING: no channel is complete in every analysed unit, ' ...
                    'so the global statistic is undefined (NaN). The channel-wise maps above ' ...
                    'are unaffected.'];
            elseif res.n_channels < 0.75 * res.n_channels_total
                results_text{end+1} = sprintf(['WARNING: the global rests on %.0f%% of the ' ...
                    'montage, so it is a regional average, not a whole-head result - read it ' ...
                    'with the exclusion pattern above.'], ...
                    100 * res.n_channels / res.n_channels_total);
            end
            results_struct.global_n_channels       = res.n_channels;
            results_struct.global_n_channels_total = res.n_channels_total;
        end

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




