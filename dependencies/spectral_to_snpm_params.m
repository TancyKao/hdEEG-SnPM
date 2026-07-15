function [params, n] = spectral_to_snpm_params(DATA, comparison, o)
%SPECTRAL_TO_SNPM_PARAMS  Slice multi-folder spectral data into core_snpm_analysis params.
% Shared by run_spectral_analysis.m (headless) and SnPMAnalysisGui (spectral mode).
%
% The spectral-folder source is multi-folder: each folder = one level of the
% design factor (load_spectral_dataset emits a 'level' column = folder label).
% This loads one (band,type,stage[s]) cell, slices it for the chosen comparison,
% writes temp CSV(s) under o.output_path, and returns a params struct for
% core_snpm_analysis. n is the sample size.
%
% Supported comparisons (folder = factor level):
%   pairedT / onesampleT  -> 2 folders, within-subject (subjects matched by sub-XX)
%   unpairedT             -> 2 folders, between-group (independent subjects)
%   anova1                -> >=2 folders, between-group (group_col = level)
%   rmanova               -> >=2 folders, within-subject (condition_col = level)
% (ancova/regression/correlation/mixed2way/mixedmodel need per-subject covariates
%  or trial-level data -> CSV input, not this folder source.)
%
% INPUTS
%   DATA : '' (folders carry their own paths) or a fallback root with condition-*
%   o    : struct
%       .folders (cellstr) one folder per level   .labels (cellstr) level labels
%       .band, .type, .stages (cellstr)
%       .output_path, .channels, .datatype, .permutations, .tail
%       .level_A, .level_B (optional) pick/order the two levels for 2-level tests

    o = defaults(o);
    lopts = struct('band',o.band, 'type',o.type, 'stages',{o.stages});
    if isfield(o,'folders') && ~isempty(o.folders), lopts.folders = o.folders; end
    if isfield(o,'labels')  && ~isempty(o.labels),  lopts.labels  = o.labels;  end
    if isfield(o,'subjects_csv') && ~isempty(o.subjects_csv), lopts.subjects_csv = o.subjects_csv; end
    [T, info] = load_spectral_dataset(DATA, lopts);

    params = struct('comparison',comparison, 'output_path',o.output_path, ...
        'channels',o.channels, 'datatype',o.datatype, 'permutations',o.permutations, 'tail',o.tail);

    lv = unique(T.level, 'stable');        % factor levels present (= folder labels)

    switch comparison
        case {'pairedT','onesampleT'}      % within-subject, 2 levels
            [lA,lB] = two_levels(o, lv, comparison);
            [f1,f2,n] = split_AB(T, info, lA, lB, o.output_path, true);
            params.data1_file=f1; params.data2_file=f2;
            params.data1_sheet='CSV File'; params.data2_sheet='CSV File';

        case 'unpairedT'                   % between-group, 2 levels
            [lA,lB] = two_levels(o, lv, comparison);
            [f1,f2,n] = split_AB(T, info, lA, lB, o.output_path, false);
            params.data1_file=f1; params.data2_file=f2;
            params.data1_sheet='CSV File'; params.data2_sheet='CSV File';

        case 'anova1'                      % between-group, >=2 levels
            params = write_glm(params, T, info, o.output_path);
            params.group_col = 'level';
            n = numel(unique(T.Subject));

        case 'rmanova'                     % within-subject, >=2 levels
            params = write_glm(params, T, info, o.output_path);
            params.condition_col = 'level'; params.subject_col = 'Subject';
            n = numel(unique(T.Subject));

        otherwise
            error('spectral_to_snpm_params:unsupported', ...
                ['Spectral-folder source supports pairedT/onesampleT/unpairedT/anova1/rmanova; ' ...
                 'got ''%s'' (covariate/predictor/trial-level designs need CSV input).'], comparison);
    end
end

% ---------------------------------------------------------------------------
function o = defaults(o)
    if ~isfield(o,'tail') || isempty(o.tail), o.tail = 'both'; end
end

function [lA, lB] = two_levels(o, lv, comparison)
    if isfield(o,'level_A') && ~isempty(o.level_A) && isfield(o,'level_B') && ~isempty(o.level_B)
        lA = char(o.level_A); lB = char(o.level_B);
    else
        assert(numel(lv)==2, '%s needs exactly 2 folders (levels); got %d (%s).', ...
            comparison, numel(lv), strjoin(cellstr(lv), ','));
        lA = char(lv(1)); lB = char(lv(2));
    end
end

function params = write_glm(params, T, info, outdir)
    gfile = fullfile(outdir, 'glm_input.csv'); writetable(T, gfile);
    params.data_file = gfile; params.data_sheet = 'CSV File';
    params.meta_cols = info.meta_cols;
end

function [f1, f2, ns] = split_AB(T, info, lA, lB, outdir, paired)
% Split the table into two channel-only wide CSVs by the two levels. paired=true
% keeps common subjects and orders both by Subject so core_snpm_analysis matches
% them row-by-row.
    sv = string(T.level);
    A = T(strcmpi(sv, lA), :);  B = T(strcmpi(sv, lB), :);
    assert(~isempty(A) && ~isempty(B), 'levels ''%s''/''%s'' not both present.', lA, lB);
    if paired
        common = intersect(string(A.Subject), string(B.Subject), 'stable');
        A = order_by_subject(A, common);  B = order_by_subject(B, common);
    end
    ns = height(A);
    chans = info.channel_cols;
    f1 = fullfile(outdir, sprintf('A_%s.csv', matlab.lang.makeValidName(char(lA))));
    f2 = fullfile(outdir, sprintf('B_%s.csv', matlab.lang.makeValidName(char(lB))));
    writetable(A(:, chans), f1);  writetable(B(:, chans), f2);
end

function T = order_by_subject(T, subj_order)
    [~, loc] = ismember(subj_order, string(T.Subject));
    T = T(loc(loc>0), :);
end
