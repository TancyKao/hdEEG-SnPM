function [results_struct, results_text] = core_snpm_lmm(params)
% Per-channel linear mixed model (LMM) analysis with permutation-based TFCE
% and cluster correction, following Stephan et al. 2021 (Current Biology).
%
% This is the trial-level (long-format) counterpart of core_snpm_analysis.
% It is invoked when params.comparison == 'mixedmodel'. Each channel is fit
% with fitlme; the effect of interest is mapped across the scalp and
% corrected for multiple comparisons with the existing neighbour/TFCE/cluster
% machinery (snpm_lmm_TFCE, snpm_lmm_cluster). It supports k>2 groups
% natively via a categorical group factor.
%
% INPUT: params struct
%   Data:
%     .data_file / .data1_file  - long-format file (one row per trial)
%     .data_sheet / .data1_sheet- sheet name, or 'CSV File' for csv
%     .output_path              - output directory
%     .channels                 - '164 channels' | '178 channels'
%     .datatype                 - 'absolute' | 'logscale' | 'normalize'
%   Model:
%     .lmm_dv          - dependent-variable column name (numeric)
%     .lmm_fixed       - fixed-effects formula, referencing 'POWER' and meta
%                        columns, e.g. 'POWER + group'
%     .lmm_random      - random-effects formula, e.g. '(1|Subject) + (1|time)'
%     .lmm_effect      - term whose effect is mapped ('POWER', 'group', ...)
%     .lmm_effect_type - 'continuous' (signed t) | 'factor' (F)
%     .lmm_meta_cols   - cellstr of non-channel columns (DV, Subject, group,
%                        time, ...). All other columns are treated as channels.
%     .lmm_perm        - 'within_subject' | 'group_label' (optional; default
%                        derived from effect_type)
%     .lmm_subject     - subject column name (optional, default 'Subject')
%     .lmm_group       - group column name (optional, default 'group')
%     .tail            - 'both' (default) | 'left' | 'right'
%     .permutations    - number of permutations
%
% OUTPUT: results_struct, results_text (same shape as core_snpm_analysis)

    %% Resolve inputs / defaults
    data_file  = getfield_default(params, 'data_file',  getfield_default(params, 'data1_file', ''));
    data_sheet = getfield_default(params, 'data_sheet', getfield_default(params, 'data1_sheet', 'CSV File'));
    tail       = getfield_default(params, 'tail', 'both');
    subj_col   = getfield_default(params, 'lmm_subject', 'Subject');
    grp_col    = getfield_default(params, 'lmm_group', 'group');

    if ~isfield(params, 'lmm_meta_cols') || isempty(params.lmm_meta_cols)
        error('core_snpm_lmm:noMetaCols', ...
            'params.lmm_meta_cols must list the non-channel columns (DV, Subject, group, time, ...).');
    end

    % Permutation scheme: derive from effect type if not given.
    if isfield(params, 'lmm_perm') && ~isempty(params.lmm_perm)
        perm_scheme = params.lmm_perm;
    elseif strcmpi(params.lmm_effect_type, 'factor')
        perm_scheme = 'group_label';
    else
        perm_scheme = 'within_subject';
    end

    %% Load long-format data
    if strcmp(data_sheet, 'CSV File')
        Tbl = readtable(data_file);
    else
        Tbl = readtable(data_file, 'Sheet', data_sheet);
    end

    varnames     = Tbl.Properties.VariableNames;
    meta_cols    = params.lmm_meta_cols;
    channel_cols = setdiff(varnames, meta_cols, 'stable');
    if isempty(channel_cols)
        error('core_snpm_lmm:noChannels', 'No channel columns left after removing meta columns.');
    end

    meta  = Tbl(:, meta_cols);
    power = Tbl{:, channel_cols};
    fprintf('LMM: %d trials, %d channel columns, %d meta columns\n', ...
        size(power, 1), numel(channel_cols), numel(meta_cols));

    %% Data type transform (applied to channel power, mirrors core_snpm_analysis)
    switch params.datatype
        case 'absolute'
            % no transformation (recommended: LMM uses absolute power)
        case 'logscale'
            power = log10(power);
        case 'normalize'
            power = bsxfun(@minus,  power, mean(power, 2, 'omitnan'));
            power = bsxfun(@rdivide, power, std(power, [], 2, 'omitnan'));
    end

    %% Channel setup (shared label-based montage selection)
    chsel = snpm_setup_channels(params.channels, {power}, channel_cols);
    power = chsel.data{1}; chanlocs = chsel.chanlocs;
    neighbors = chsel.neighbors; insidegoodch = chsel.insidegoodch;
    select_mode = chsel.select_mode;
    fprintf('%s montage: %d/%d channel columns matched to locations\n', ...
        params.channels, numel(insidegoodch), numel(channel_cols));

    % The TFCE/cluster machinery requires the channel-statistic vector to
    % align with the neighbour adjacency: size(power,2) must equal
    % size(neighbors,1). Fail fast with a clear message otherwise.
    if size(power, 2) ~= size(neighbors, 1)
        error('core_snpm_lmm:channelMismatch', ...
            ['Channel count (%d) does not match the neighbour matrix (%d rows). ' ...
             'For the ''%s'' montage the long-format file must have exactly %d channel columns ' ...
             '(same channel set as the existing wide-format tool).'], ...
            size(power, 2), size(neighbors, 1), params.channels, size(neighbors, 1));
    end

    %% Build the LMM spec
    spec = struct();
    spec.dv          = params.lmm_dv;
    spec.fixed       = params.lmm_fixed;
    spec.random      = params.lmm_random;
    spec.effect      = params.lmm_effect;
    spec.effect_type = params.lmm_effect_type;
    spec.subject     = subj_col;
    spec.group       = grp_col;
    spec.perm        = perm_scheme;
    spec.categorical = lmm_categorical_vars(meta, subj_col, grp_col, spec.random);

    alpha = 0.05;
    E = 0.5; H = 2;

    fprintf('LMM formula: %s ~ %s + %s\n', spec.dv, spec.fixed, spec.random);
    fprintf('Effect of interest: %s (%s); permutation scheme: %s; tail: %s\n', ...
        spec.effect, spec.effect_type, spec.perm, tail);

    %% Statistical analysis
    [T, p] = snpm_lmm_TFCE(power, meta, spec, neighbors, E, H, alpha, tail, params.permutations);
    Clusters = snpm_lmm_cluster(power, meta, spec, neighbors, alpha, params.permutations);

    %% Significant channels (same derivation as core_snpm_analysis:488-497)
    uncorrsigch      = find(p.real <= alpha);
    correctTFCEsigch = find(p.correctedTFCE <= alpha);
    sigclusters      = find([Clusters.p] <= alpha);
    SnPMsigch        = [Clusters(sigclusters).channels];

    %% Output naming
    [~, data_basename] = fileparts(data_file);
    base_filename = sprintf('%s_LMM_%s_on_%s', data_basename, spec.effect, spec.dv);
    base_filename = matlab.lang.makeValidName(base_filename);

    if ~exist(params.output_path, 'dir')
        mkdir(params.output_path);
    end
    timestamp   = datestr(now, 'yyyymmdd_HHMMSS');
    outputSname = fullfile(params.output_path, [base_filename, '_', timestamp]);
    savepath    = params.output_path;

    %% Plot significance topography (reused output function)
    try
        TopoplotSignificant_single(T.real_T, uncorrsigch, correctTFCEsigch, SnPMsigch, ...
            chanlocs, insidegoodch, params.comparison, savepath, base_filename, select_mode);
    catch ME
        warning(ME.identifier, 'Could not generate significance topology plots: %s', ME.message);
    end

    %% Results struct
    results_struct = struct();
    results_struct.T                = T;
    results_struct.p                = p;
    results_struct.Clusters         = Clusters;
    results_struct.uncorrsigch      = uncorrsigch;
    results_struct.correctTFCEsigch = correctTFCEsigch;
    results_struct.SnPMsigch        = SnPMsigch;
    results_struct.chanlocs         = chanlocs;

    nsubj = numel(unique(meta.(subj_col)));
    results_struct.data_summary.data1_size = size(power);
    results_struct.data_summary.data2_size = [nsubj, size(power, 1)];   % subjects, trials
    results_struct.data_summary.data1_mean = mean(power(:), 'omitnan');
    results_struct.data_summary.data2_mean = NaN;
    results_struct.data_summary.data1_std  = std(mean(power, 2, 'omitnan'), 'omitnan');
    results_struct.data_summary.data2_std  = NaN;

    results_struct.lmm.formula      = sprintf('%s ~ %s + %s', spec.dv, spec.fixed, spec.random);
    results_struct.lmm.effect       = spec.effect;
    results_struct.lmm.effect_type  = spec.effect_type;
    results_struct.lmm.perm_scheme  = spec.perm;
    results_struct.lmm.n_trials     = size(power, 1);
    results_struct.lmm.n_subjects   = nsubj;
    if ismember(grp_col, meta.Properties.VariableNames)
        results_struct.lmm.group_levels = categories(categorical(meta.(grp_col)));
    else
        results_struct.lmm.group_levels = {};
    end
    results_struct.lmm.dv     = spec.dv;
    results_struct.lmm.fixed  = spec.fixed;
    results_struct.lmm.random = spec.random;
    results_struct.lmm.tail   = tail;

    %% Effect-map + descriptive PNGs + global test for the LMM report
    is_factor = strcmpi(spec.effect_type, 'factor');
    descch = SnPMsigch; if isempty(descch), descch = correctTFCEsigch; end
    if isempty(descch), descch = 1:size(power, 2); end
    results_struct.lmm.effect_map_png       = '';
    results_struct.lmm.descriptive_png      = '';
    results_struct.lmm.descriptive_channels = numel(descch);
    results_struct.lmm.global_stat          = NaN;
    results_struct.lmm.global_pval          = NaN;

    % Effect map: jet + symmetric clim for signed t, hot for F.
    try
        tv = T.real_T; empng = [outputSname '_effectmap.png'];
        if is_factor
            save_single_topo(tv, chanlocs, 'Effect map (F)', empng, 'hot');
        else
            m = max(abs(tv(:))); if m == 0 || ~isfinite(m), m = 1; end
            save_single_topo(tv, chanlocs, 'Effect map (t)', empng, 'jet', [-m m]);
        end
        [~, fn, ext] = fileparts(empng); results_struct.lmm.effect_map_png = [fn ext];
    catch ME
        warning(ME.identifier, 'Could not generate LMM effect map: %s', ME.message);
    end

    % Descriptive: scatter (continuous) or group-means (factor), over descch.
    try
        dvcol = spec.dv; dpng = [outputSname '_descriptive.png'];
        cavg  = mean(power(:, descch), 2, 'omitnan');   % cluster-averaged power per trial
        gv    = lmm_effect_grouping(meta, spec);
        if is_factor
            if ismember(dvcol, meta.Properties.VariableNames) && isnumeric(meta.(dvcol))
                resp = double(meta.(dvcol));
            else
                resp = cavg;   % DV is POWER (the channel) -> summarise cluster-averaged power
            end
            if ~isempty(gv)
                plot_lmm_groupmeans(resp, gv, dvcol, dpng);
                [~, fn, ext] = fileparts(dpng); results_struct.lmm.descriptive_png = [fn ext];
            end
        else
            if ismember(dvcol, meta.Properties.VariableNames) && isnumeric(meta.(dvcol))
                plot_lmm_scatter(cavg, double(meta.(dvcol)), gv, spec.effect, dvcol, dpng);
                [~, fn, ext] = fileparts(dpng); results_struct.lmm.descriptive_png = [fn ext];
            end
        end
    catch ME
        warning(ME.identifier, 'Could not generate LMM descriptive plot: %s', ME.message);
    end

    % Global whole-head test: fit the effect on the channel-averaged signal.
    try
        [gstat, ~, gp] = snpm_lmm_fit(mean(power, 2, 'omitnan'), meta, spec);
        results_struct.lmm.global_stat = gstat(1);
        results_struct.lmm.global_pval = gp(1);
    catch ME
        warning(ME.identifier, 'Could not compute LMM global test: %s', ME.message);
    end

    %% Save .mat
    save([outputSname '.mat'], 'results_struct');

    %% Excel + HTML (reused output functions)
    try
        func_genSnpmTable(results_struct, uncorrsigch, correctTFCEsigch, SnPMsigch, chanlocs, outputSname, select_mode);
        disp('Excel file generated successfully!');
    catch ME
        warning(ME.identifier, 'Could not generate Excel file: %s', ME.message);
    end

    %% Text results
    results_text = {
        '=== SnPM Mixed Linear Model Results ===', '', ...
        ['Data file: ' data_file], ...
        ['Model: ' results_struct.lmm.formula], ...
        ['Effect of interest: ' spec.effect ' (' spec.effect_type ')'], ...
        ['Permutation scheme: ' spec.perm], ...
        ['Tail: ' tail], ...
        ['Permutations: ' num2str(params.permutations)], ...
        ['Trials: ' num2str(results_struct.lmm.n_trials) ', Subjects: ' num2str(nsubj)], ...
        '', '--- Statistical Results ---', ...
        ['Uncorrected significant channels: ' num2str(length(uncorrsigch))], ...
        ['TFCE corrected significant channels: ' num2str(length(correctTFCEsigch))], ...
        ['Significant clusters N = : ' num2str(length(sigclusters))], ...
        ['SnPM significant channels N = : ' num2str(length(SnPMsigch))], ...
        ['SnPM significant channels = : ' num2str(SnPMsigch)], ...
        '', ['Results saved to: ' outputSname '.mat'], ...
        ['Excel file saved to: ' outputSname '.xlsx'] ...
        }';

    %% HTML reports (reused output function)
    try
        fprintf('Generating HTML reports...\n');
        generateAnalysisReport_lmm(results_struct, params, base_filename, outputSname, uncorrsigch, correctTFCEsigch, SnPMsigch);
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

function v = getfield_default(s, f, default)
    if isfield(s, f) && ~isempty(s.(f))
        v = s.(f);
    else
        v = default;
    end
end

function gv = lmm_effect_grouping(meta, spec)
% Categorical grouping column implied by the mapped effect, for the descriptive
% plot: the meta column named in spec.effect (or, for an interaction a:b / a*b,
% the factor part that is a meta column, ignoring the POWER token and the DV).
% Returns [] when the effect has no meta grouping (a plain continuous slope).
    gv = [];
    toks = regexp(spec.effect, '[:*]', 'split');
    for i = 1:numel(toks)
        t = strtrim(toks{i});
        if isempty(t) || strcmpi(t, 'POWER') || strcmpi(t, spec.dv), continue; end
        if ismember(t, meta.Properties.VariableNames)
            gv = meta.(t); return;
        end
    end
end

function cats = lmm_categorical_vars(meta, subj_col, grp_col, random_formula)
% Variables to force categorical: subject, group, and any random-effect
% grouping variable (the token after '|' in the random formula).
    cats = {subj_col, grp_col};
    tok = regexp(random_formula, '\|\s*([A-Za-z]\w*)', 'tokens');
    for i = 1:numel(tok)
        cats{end+1} = tok{i}{1}; %#ok<AGROW>
    end
    cats = intersect(cats, meta.Properties.VariableNames, 'stable');
end
