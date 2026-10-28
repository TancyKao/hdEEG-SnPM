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
%     .tail            - 'both' (default; the only value accepted, see below)
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

    %% Up-front refusals: settings under which the permutation test cannot work
    % (checked before any file is read or any model fitted; same codes and
    % texts as the Python port).
    if ~strcmp(tail, 'both')
        error('core_snpm_lmm:oneSidedTail', ...
            ['mixedmodel requires tail = ''both'' (got ''%s''). The effect map is a two-sided t or ' ...
             'an F, and clusters are formed on the two-sided model p; a one-sided TFCE null of ' ...
             'signed maxima would not match it.'], char(tail));
    end
    % Cluster statistic ('sum' default | 'mean'), refused before any file is read.
    [cs, cs_line] = snpm_cluster_statistic(params);
    params.cluster_statistic = cs;
    eff_spec = struct('dv', params.lmm_dv, 'fixed', params.lmm_fixed, 'effect', params.lmm_effect, ...
        'effect_type', params.lmm_effect_type, 'perm', perm_scheme, 'group', grp_col);
    if strcmpi(perm_scheme, 'group_label')
        % coarse form here (the group levels are not known yet); the exact
        % '<group>_<level>' rule runs after the file is read, below
        snpm_lmm_check_effect([], eff_spec);
    else
        dv = char(params.lmm_dv);
        if strcmp(dv, 'POWER') || ~ismember(dv, params.lmm_meta_cols)
            if strcmp(dv, 'POWER'), why = '(it is the channel value).'; else, why = 'among lmm_meta_cols.'; end
            error('core_snpm_lmm:dvNotMetadata', ...
                ['The within_subject scheme shuffles the response within each subject, so the ' ...
                 'response must be a metadata column; ''%s'' is not %s Use group_label for a ' ...
                 'between-group effect on the channel value.'], dv, why);
        end
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

    %% Rows with a missing model variable: dropped ONCE, before anything else
    % (transform, channel setup, evaluable mask), so the observed fit and every
    % permutation use one row set; then the design is checked against the
    % permutation scheme (group constant within subject for group_label, >= 3
    % subjects, ...). See snpm_lmm_drop_missing_meta.
    drop_spec = struct('dv', params.lmm_dv, 'fixed', params.lmm_fixed, 'random', params.lmm_random, ...
        'subject', subj_col, 'group', grp_col, 'perm', perm_scheme);
    [meta, power, dropped_rows] = snpm_lmm_drop_missing_meta(meta, power, drop_spec);
    if dropped_rows.n > 0
        fprintf('LMM: %d trials left after dropping %d with a missing model variable\n', ...
            size(power, 1), dropped_rows.n);
    end

    % within_subject shuffles the DV inside each subject: a DV that is constant
    % within every subject leaves nothing to permute (every permuted map would
    % equal the observed one). Under group_label a metadata DV constant within
    % every subject makes every fit a perfect fit (the backstop would blank the
    % whole map): it is a subject-level analysis. Labelling-independent, so
    % decided once, here.
    if ~strcmp(char(params.lmm_dv), 'POWER') && ismember(char(params.lmm_dv), meta.Properties.VariableNames)
        w_dv = snpm_lmm_within_ratio(meta.(params.lmm_dv), meta.(subj_col));
        if w_dv <= 1e-6 && strcmpi(perm_scheme, 'within_subject')
            error('core_snpm_lmm:dvConstantWithinSubject', ...
                ['The within_subject scheme shuffles ''%s'' within each subject, but it is ' ...
                 'constant within each subject (within-subject share of its variation %.3g, ' ...
                 'threshold 1e-06), so the within-subject shuffle has nothing to permute. ' ...
                 'Use group_label for a between-group effect, or analyse subject means.'], ...
                char(params.lmm_dv), w_dv);
        elseif w_dv <= 1e-6
            error('core_snpm_lmm:dvConstantWithinSubject', ...
                ['The response ''%s'' is constant within each subject (within-subject share of ' ...
                 'its variation %.3g, threshold 1e-06), so every mixed-model fit is a perfect fit ' ...
                 'and no statistic is defined. Analyse subject means with anova1/unpairedT.'], ...
                char(params.lmm_dv), w_dv);
        end
    end

    % The effect must be testable (group_label: an actual level of the group
    % column) and must name exactly one coefficient (continuous) or a term
    % (factor) of the fixed part, and the fixed part must be full rank on the
    % metadata; otherwise every channel would come back NaN, or a different
    % hypothesis would be mapped. Decided from the metadata alone, before the
    % channel setup and any fit. The effect is returned as fitlme names it
    % ('POWER:group' -> 'group:POWER' for a factor; a two-level 'group' ->
    % 'group_B' for a continuous effect), and snpm_lmm_fit matches it exactly.
    eff_spec.categorical = lmm_categorical_vars(meta, subj_col, grp_col, params.lmm_random);
    effect_as_given = char(params.lmm_effect);
    params.lmm_effect = snpm_lmm_check_effect(meta, eff_spec);

    % within_subject shuffles the DV inside each subject, which leaves every
    % subject's mean of it in place: an effect whose variables are all constant
    % within each subject (group, group_B, a subject-level covariate) barely
    % moves under the shuffle, and the permuted statistics are the observed one
    % plus a systematic jitter -- there is no null for it (measured false-positive
    % rate 0.35 at a nominal 0.05 when POWER predicts the response). An effect
    % with POWER, or with any variable that varies within subject (group:time),
    % runs. See snpm_lmm_effect_within_share.
    if strcmpi(perm_scheme, 'within_subject')
        w_eff = snpm_lmm_effect_within_share(meta, effect_as_given, subj_col);
        if w_eff <= 1e-6
            error('core_snpm_lmm:effectConstantWithinSubject', ...
                ['The within_subject scheme shuffles ''%s'' within each subject, but the effect ' ...
                 '''%s'' is constant within each subject (largest within-subject share of ' ...
                 'variation among its variables %.3g, threshold 1e-06), so the shuffle leaves the ' ...
                 'effect in place and gives no null for it. Use group_label for a between-group ' ...
                 'effect, or test an effect that varies within subjects.'], ...
                char(params.lmm_dv), effect_as_given, w_eff);
        end
    end

    %% Data type transform (applied to channel power, mirrors core_snpm_analysis)
    % Capture raw negativity BEFORE transforming: logscale/normalize legitimately
    % produce negatives, so the source magnitude guard below must inspect the RAW
    % input, not the transformed power.
    raw_has_negative = any(power(:) < 0);
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
    is_source = isfield(chsel, 'is_source') && chsel.is_source;
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

    % Source-space (2447 cortical voxels): enforce count + node-order invariants
    % and reject signed current density (band power must be a non-negative
    % magnitude) when the DV is a source voxel column. This is Stephan et al.'s
    % source replication; the stat engine is unchanged.
    if is_source
        snpm_assert_source(chanlocs, neighbors, chsel.channel_mapping.channel_labels);
        if strcmp(params.datatype, 'absolute')
            if raw_has_negative
                warning('snpm:source:signedMagnitude', ...
                    ['Source voxel matrix contains negative values but non-negative ' ...
                     'band-power magnitudes are expected. Signed current density is ' ...
                     'invalid (arbitrary per-subject/trial sign -> cancellation). ' ...
                     'Confirm the DV columns are voxel POWER (magnitude).']);
            end
        elseif raw_has_negative
            warning('snpm:source:signedUnderTransform', ...
                ['Source input contains negative values with datatype=%s; if this ' ...
                 'is signed current density rather than a magnitude, the group test ' ...
                 'is invalid (arbitrary per-subject/trial sign + cancellation).'], ...
                params.datatype);
        end
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

    %% Evaluable-channel mask (null construction, decided ONCE)
    % A channel is evaluable if and only if it is complete across every
    % analysed row. snpm_lmm_fit deletes rows per channel, so an incomplete
    % channel is fitted on a subset while the permutation relabels/shuffles
    % across ALL rows — the permuted fits then see a different design (a group
    % factor can become constant on the subset -> anova returns FStat = NaN
    % without throwing) or a different DV distribution (when missingness is
    % associated with the DV) than the observed fit. Restricting to complete
    % channels makes `valid` all-true inside the fit, so the observed map and
    % every permuted map are computed on identical rows and identical channels.
    evaluable = all(isfinite(power), 1);
    % Shared exclusion surface (same struct the GLM/circular tiers emit, which
    % func_genSnpmTable and the report banners already know how to render).
    % Only the reason string is specialised: in the long format the rule is
    % completeness across analysed ROWS (trials), not across subjects -- a
    % channel missing one trial of an otherwise present subject is excluded.
    excluded_channels = snpm_excluded_channel_info(evaluable, chanlocs, ...
        sprintf('mixedmodel (%s)', spec.effect));
    excluded_channels.reason = ['missing data: channel not complete across all analysed ' ...
        'rows (trials), so the fitted row set would differ from the permuted one'];
    if ~any(evaluable)
        error('core_snpm_lmm:noEvaluableChannels', ...
            ['No channel is complete across all %d analysed rows, so no channel can be ' ...
             'tested against a permutation null. Drop the incomplete subjects/trials ' ...
             'instead of the channels.'], size(power, 1));
    end
    if excluded_channels.n > 0.10 * numel(evaluable)
        warning('core_snpm_lmm:manyChannelsExcluded', ...
            ['%d of %d channels (%.1f%%) are incomplete and were dropped. Above ~10%% of ' ...
             'the montage, dropping the few SUBJECTS/trials that carry the missingness ' ...
             'usually costs less information than dropping the channels — consider ' ...
             'filtering the input table instead.'], excluded_channels.n, numel(evaluable), ...
            100 * excluded_channels.n / numel(evaluable));
    end

    % Channels whose value is constant within each subject (within-subject share
    % of variation w <= 1e-6, snpm_lmm_within_ratio). Grouped by subject only,
    % so the property does not depend on the labelling: decided ONCE, here, and
    % held fixed for the observed map and every permutation.
    %   - channel is the response (dv = POWER, the group_label case): nothing
    %     within subject for the mixed model to estimate; fitlme cannot resolve
    %     such a fit (its subject/residual SD ratio runs to infinity).
    %   - channel is the predictor under within_subject: its slope is
    %     between-subject only, the shuffle cannot move it, and its t would sit
    %     in every permuted map unchanged, raising the max null for the rest.
    w_ch = nan(1, numel(evaluable));
    w_ch(evaluable) = snpm_lmm_within_ratio(power(:, evaluable), meta.(subj_col));
    const_ch = evaluable & w_ch <= 1e-6;
    if strcmp(spec.dv, 'POWER')
        excluded_channels = snpm_excluded_add(excluded_channels, const_ch, lmm_reason_const_response(), chanlocs);
        evaluable = evaluable & ~const_ch;
    elseif strcmpi(spec.perm, 'within_subject')
        excluded_channels = snpm_excluded_add(excluded_channels, const_ch, lmm_reason_const_predictor(), chanlocs);
        evaluable = evaluable & ~const_ch;
    else
        % group_label with a metadata response and POWER as predictor. Relabelling
        % whole subjects moves each subject's response and POWER together, so a
        % POWER constant within subject is a subject-level covariate and the
        % channel is KEPT (measured null false-positive rate nominal). Two
        % exceptions, both labelling-independent and decided once:
        %   - a channel with one value on every trial: POWER is collinear with
        %     the intercept and fitlme refuses the fit ("must be of full column
        %     rank") at every labelling, which used to leave it NaN unlisted;
        %   - a random slope on POWER by subject: its variance cannot be
        %     estimated where POWER is constant within subject (the fit sits on
        %     the variance boundary at every labelling).
        flat_ch = false(1, numel(evaluable));
        flat_ch(evaluable) = lmm_flat_columns(power(:, evaluable));
        excluded_channels = snpm_excluded_add(excluded_channels, flat_ch, lmm_reason_flat(), chanlocs);
        evaluable = evaluable & ~flat_ch;
        if lmm_power_slope_by_subject(spec.random, subj_col)
            slope_ch = evaluable & w_ch <= 1e-6;
            excluded_channels = snpm_excluded_add(excluded_channels, slope_ch, lmm_reason_const_slope(), chanlocs);
            evaluable = evaluable & ~slope_ch;
        end
    end
    if ~any(evaluable)
        error('core_snpm_lmm:noEvaluableChannels', ...
            'No channel is left to test: %d of %d excluded (%s).', ...
            excluded_channels.n, excluded_channels.n_channels, excluded_channels.reason);
    end

    % Parity recorder (scripts/export_parity_fixtures.m). Inert unless switched
    % on: it consumes no random number and nothing reads it back.
    if snpm_perm_tap('ison')
        snpm_perm_tap('record', 'lmm:engine_inputs', struct('power', power, 'meta', meta, ...
            'spec', spec, 'evaluable', evaluable, 'tail', tail, 'alpha', alpha, 'E', E, 'H', H, ...
            'keptcols', chsel.channel_mapping.original_data_columns, ...
            'channel_cols', {channel_cols}));
    end

    %% Statistical analysis
    [T, p] = snpm_lmm_TFCE(power, meta, spec, neighbors, E, H, alpha, tail, params.permutations, evaluable);
    Clusters = snpm_lmm_cluster(power, meta, spec, neighbors, alpha, params.permutations, evaluable, cs);
    % An observed channel whose fit had no usable residual variation
    % (snpm_lmm_degenerate, the per-labelling backstop) is NaN in the map; list it.
    und = false(1, numel(evaluable)); und(T.undefined_channels) = true;
    excluded_channels = snpm_excluded_add(excluded_channels, und, 'perfect fit, statistic undefined', chanlocs);

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
    % Source-space systems have no scalp layout -> skip the topoplot (a voxel
    % list with coordinates is emitted after the .mat/.xlsx below).
    if ~is_source
        try
            % Flat channel list (unnumbered cluster PNG) on purpose: the LMM tier
            % still renders through generateAnalysisReport_lmm, whose cluster table
            % numbers by ascending p with no positive/negative split, so a
            % sign-ordered PNG number would disagree with the table number. Pass the
            % struct('Clusters',...,'alpha',...) form here when the LMM tier moves to
            % the v2 report generator (phase 2).
            TopoplotSignificant_single(T.real_T, uncorrsigch, correctTFCEsigch, SnPMsigch, ...
                chanlocs, insidegoodch, params.comparison, savepath, base_filename, select_mode);
        catch ME
            warning(ME.identifier, 'Could not generate significance topology plots: %s', ME.message);
        end
    else
        fprintf('Source-space system: scalp topoplot bypassed; emitting voxel list.\n');
    end

    %% Results struct
    results_struct = struct();
    results_struct.T                = T;
    results_struct.p                = p;
    results_struct.Clusters         = Clusters;
    results_struct.uncorrsigch      = uncorrsigch;
    results_struct.correctTFCEsigch = correctTFCEsigch;
    results_struct.SnPMsigch        = SnPMsigch;
    % The significance level the maps, clusters and TFCE sets above were built
    % at. Published so the report numbers its clusters with the SAME threshold
    % the run used, rather than a second hard-coded .05 that could drift.
    results_struct.alpha = alpha;
    results_struct.cluster_statistic = cs;   % 'sum' | 'mean' (snpm_cluster_statistic)
    results_struct.chanlocs         = chanlocs;
    results_struct.is_source        = is_source;   % source-space run (topoplot bypassed)
    % Channel-exclusion surface shared with the GLM tier: .mat field ->
    % func_genSnpmTable 'excludedChannels' sheet -> HTML banner.
    results_struct.excluded_channels = excluded_channels;

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
    results_struct.lmm.dropped_rows = dropped_rows;   % rows dropped for a missing model variable
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
    if isempty(descch), descch = find(evaluable); end   % analysed set, not all columns
    results_struct.lmm.effect_map_png       = '';
    results_struct.lmm.descriptive_png      = '';
    results_struct.lmm.descriptive_channels = numel(descch);
    results_struct.lmm.global_stat          = NaN;
    results_struct.lmm.global_pval          = NaN;

    % Effect map: jet + symmetric clim for signed t, hot for F.
    % Source-space: no scalp layout -> skip (voxel list carries the per-voxel stat).
    if ~is_source
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

    % Global whole-head test: fit the effect on the channel-averaged signal,
    % averaged over the ANALYSED channels only (the same set the maps use).
    try
        [gstat, ~, gp] = snpm_lmm_fit(mean(power(:, evaluable), 2, 'omitnan'), meta, spec);
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
        snpm_write_settings_sheet(outputSname, cs, alpha);
        disp('Excel file generated successfully!');
    catch ME
        warning(ME.identifier, 'Could not generate Excel file: %s', ME.message);
    end

    % Source-space deliverable: significant-voxel list with coordinates +
    % cluster membership (replaces the scalp topoplot bypassed above).
    if is_source
        try
            write_source_voxel_table(results_struct, chanlocs, outputSname);
        catch ME
            warning(ME.identifier, 'Could not write source voxel table: %s', ME.message);
        end
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
        cs_line, ...
        ['Trials: ' num2str(results_struct.lmm.n_trials) ', Subjects: ' num2str(nsubj)], ...
        snpm_lmm_dropped_text(dropped_rows), ...
        '', '--- Statistical Results ---', ...
        ['Uncorrected significant channels: ' num2str(length(uncorrsigch))], ...
        ['TFCE corrected significant channels: ' num2str(length(correctTFCEsigch))], ...
        ['Significant clusters N = : ' num2str(length(sigclusters))], ...
        ['SnPM significant channels N = : ' num2str(length(SnPMsigch))], ...
        ['SnPM significant channels = : ' num2str(SnPMsigch)], ...
        lmm_excluded_text(excluded_channels), ...
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

% Reason strings are compared with the Python port; do not reword.
function s = lmm_reason_const_response()
    s = ['channel value constant within each subject: no within-subject variation for the mixed ' ...
         'model to estimate; analyse subject means with anova1/unpairedT'];
end

function s = lmm_reason_const_predictor()
    s = ['predictor constant within each subject: its effect is between-subject only, which the ' ...
         'within-subject shuffle cannot test (use group_label or regression on subject means)'];
end

function s = lmm_reason_const_slope()
    s = 'predictor constant within each subject: its random slope by subject cannot be estimated';
end

function s = lmm_reason_flat()
    s = 'same value in every subject';
end

function tf = lmm_flat_columns(V)
% one value on every row: total sum of squares <= 1e-20 x uncentred (the flat
% test of snpm_lmm_within_ratio)
    tf = sum((V - mean(V, 1)).^2, 1) <= 1e-20 * sum(V.^2, 1);
end

function tf = lmm_power_slope_by_subject(random_formula, subj_col)
% true when the random part has a term (... POWER ... | <subject>), e.g.
% (POWER|Subject) or (1 + POWER|Subject)
    tf = false;
    tok = regexp(char(random_formula), '\(([^|()]*)\|\s*(\w+)\s*\)', 'tokens');
    for i = 1:numel(tok)
        if strcmp(tok{i}{2}, subj_col) && ...
                any(strcmp(regexp(tok{i}{1}, '[A-Za-z]\w*', 'match'), 'POWER'))
            tf = true;
            return
        end
    end
end

function txt = lmm_excluded_text(excl)
% One results_text line (feeds the GUI Result box + the report CSV).
    if excl.n == 0
        txt = sprintf('Channels excluded for missing data: none (%d/%d evaluable)', ...
            excl.n_channels, excl.n_channels);
    elseif all(startsWith(excl.reasons, 'missing data'))   % .reasons: always set by snpm_excluded_add
        txt = sprintf('Channels excluded for missing data (%s): %d of %d — %s', ...
            excl.context, excl.n, excl.n_channels, strjoin(excl.labels, ', '));
    else
        % more than one rule can exclude: 'E1, E2 (reason A); E7 (reason B)'
        u = unique(excl.reasons, 'stable');
        parts = cell(1, numel(u));
        for k = 1:numel(u)
            parts{k} = sprintf('%s (%s)', strjoin(excl.labels(strcmp(excl.reasons, u{k})), ', '), u{k});
        end
        txt = sprintf('Channels excluded (%s): %d of %d — %s', ...
            excl.context, excl.n, excl.n_channels, strjoin(parts, '; '));
    end
end

function gv = lmm_effect_grouping(meta, spec)
% Categorical grouping column implied by the mapped effect, for the descriptive
% plot: the meta column named in spec.effect (or, for an interaction a:b / a*b,
% the factor part that is a meta column, ignoring the POWER token and the DV).
% Returns [] when the effect has no meta grouping (a plain continuous slope).
    gv = [];
    vn = meta.Properties.VariableNames;
    toks = regexp(spec.effect, '[:*]', 'split');
    for i = 1:numel(toks)
        t = strtrim(toks{i});
        if isempty(t) || strcmpi(t, 'POWER') || strcmpi(t, spec.dv), continue; end
        if ~ismember(t, vn)
            % a coefficient name (group_B) -> its column (the longest c with
            % the token starting with c_), as snpm_lmm_effect_within_share maps it
            c = vn(cellfun(@(v) startsWith(t, [v '_']), vn));
            if isempty(c), continue; end
            [~, m] = max(cellfun(@numel, c));
            t = c{m};
            if strcmpi(t, spec.dv), continue; end
        end
        gv = meta.(t); return;
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
