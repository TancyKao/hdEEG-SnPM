function [results_struct, results_text] = core_snpm_glm(params)
% General Linear Model group analyses with permutation TFCE + cluster
% correction. Trial- and subject-level designs that the legacy t-test/
% correlation engines cannot express: >2 groups (one-way ANOVA), ANCOVA,
% regression with nuisance covariates, repeated-measures ANOVA across >2
% conditions, and group x condition interaction.
%
% Invoked from core_snpm_analysis when params.comparison is one of:
%   'anova1' | 'ancova' | 'regression' | 'rmanova' | 'mixed2way'
% The user only maps columns; the statistic (t/F), contrast, and permutation
% scheme are chosen automatically by snpm_glm_design.
%
% INPUT: params struct
%   Data:  .data_file/.data1_file, .data_sheet/.data1_sheet, .output_path,
%          .channels ('164 channels'|'178 channels'), .datatype, .permutations
%          .meta_cols (cellstr) - non-channel columns; the rest are channels
%   Column mapping (only those the preset needs):
%          .group_col, .condition_col, .subject_col, .predictor_col,
%          .covariate_cols (cellstr), .effect ('interaction'|'condition' for mixed2way)
%   .posthoc (logical, default true) - pairwise post-hoc maps for omnibus F
%
% OUTPUT: results_struct, results_text (same shape as core_snpm_analysis)

    preset = params.comparison;

    data_file  = getfield_default(params, 'data_file',  getfield_default(params, 'data1_file', ''));
    data_sheet = getfield_default(params, 'data_sheet', getfield_default(params, 'data1_sheet', 'CSV File'));
    do_posthoc = getfield_default(params, 'posthoc', true);

    if ~isfield(params, 'meta_cols') || isempty(params.meta_cols)
        error('core_snpm_glm:noMetaCols', 'params.meta_cols must list the non-channel columns.');
    end
    % ANCOVA with no covariate is a one-way ANOVA that would be labelled
    % "covariate-adjusted". Refused before the file is read (text shared with
    % the Python port; do not reword).
    snpm_glm_require_covariate(preset, getfield_default(params, 'covariate_cols', {}));
    % Cluster statistic ('sum' default | 'mean'), validated before the file is
    % read and normalised in params so every effect below reads one value.
    [cs, cs_line] = snpm_cluster_statistic(params);
    params.cluster_statistic = cs;

    %% Load
    if strcmp(data_sheet, 'CSV File')
        Tbl = readtable(data_file);
    else
        Tbl = readtable(data_file, 'Sheet', data_sheet);
    end
    channel_cols = setdiff(Tbl.Properties.VariableNames, params.meta_cols, 'stable');
    if isempty(channel_cols)
        error('core_snpm_glm:noChannels', 'No channel columns after removing meta columns.');
    end
    meta  = Tbl(:, params.meta_cols);
    power = Tbl{:, channel_cols};
    fprintf('GLM (%s): %d observations, %d channel columns, %d meta columns\n', ...
        preset, size(power,1), numel(channel_cols), numel(params.meta_cols));

    %% Transform
    switch params.datatype
        case 'absolute'
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

    % Source-level (2447 cortical voxels) is NOT supported on the GLM presets.
    % These designs test between-/within-subject mean effects on SIGNED current
    % density, where arbitrary per-subject sign cancels the effect and yields an
    % invalid null; the GLM path also emits no spatial source deliverable. Reject
    % at the engine (load-bearing control for both the script and GUI doors)
    % before any statistics run. Route source data through the t-test
    % (core_snpm_analysis) or LMM (core_snpm_lmm) paths instead.
    if is_source
        error('snpm:source:notSupportedInGLM', ...
            ['Source-level analysis is not supported for GLM presets ' ...
             '(anova1/ancova/regression/rmanova/mixed2way); run source data ' ...
             'through the t-test (core_snpm_analysis) or LMM (core_snpm_lmm) ' ...
             'paths.']);
    end

    if size(power, 2) ~= size(neighbors, 1)
        error('core_snpm_glm:channelMismatch', ...
            ['Channel count (%d) does not match the neighbour matrix (%d rows). For the ''%s'' ' ...
             'montage the file must have exactly %d channel columns.'], ...
            size(power, 2), size(neighbors, 1), params.channels, size(neighbors, 1));
    end

    %% Build design from preset
    opts = struct();
    for f = {'group_col','condition_col','subject_col','predictor_col','covariate_cols','effect'}
        if isfield(params, f{1}), opts.(f{1}) = params.(f{1}); end
    end
    % blank / whitespace-only covariate names are dropped ({'age',''} -> {'age'})
    if isfield(opts, 'covariate_cols'), opts.covariate_cols = snpm_glm_covariate_list(opts.covariate_cols); end

    % Two-way mixed ANOVA needs a complete report: the group main effect
    % (between-subjects), the condition main effect (within-subject) and their
    % interaction (within-subject) each have a different error term and
    % permutation scheme, so they cannot come from one design. mixed2way routes
    % to a dedicated multi-effect pipeline; the single-effect flow below (and
    % every other preset) is untouched. params.effect is ignored here — the
    % report always covers all three effects.
    if strcmpi(preset, 'mixed2way')
        [results_struct, results_text] = run_mixed2way(preset, meta, power, opts, ...
            chanlocs, neighbors, insidegoodch, select_mode, params, data_file, channel_cols);
        return;
    end

    if strcmpi(preset, 'regression')
        snpm_glm_require_numeric_predictor(preset, meta, opts.predictor_col);
    end
    D = snpm_glm_design(preset, meta, opts);
    % every listed covariate constant on these rows -> dropped -> no covariate left
    snpm_glm_require_covariate(preset, getfield_default(opts, 'covariate_cols', {}), D);

    % Between-subject presets (anova1/ancova/regression -> 'free' permutation)
    % assume ONE row per subject: whole-subject relabelling is only a valid null
    % if subjects are exchangeable. Trial-level / repeated data (a subject in
    % several rows) has within-subject correlation that breaks exchangeability
    % and silently loses family-wise-error control. Guard against it here.
    % rmanova/mixed2way (perm_type 'within', repeated measures) are exempt.
    if strcmp(D.perm_type, 'free')
        guard_between_subject_rows(meta, opts, preset);
    end

    % Effect aliased with the nuisance (a covariate repeating the predictor or
    % the group): the contrast is not estimable and pinv would return a finite
    % but arbitrary statistic, so it is refused; a nuisance-only dependency
    % (two covariates spanning the same space) warns and runs. X is fixed
    % across Freedman-Lane permutations, so the check on it is exact.
    snpm_glm_estimability(D, 'enforce');

    % Saturated design (round-4 spec 3b, approved 2026-10-14): with no error
    % degrees of freedom every channel would hit zero-variance rule 2 and rule
    % 3 would mislabel the whole montage "perfect fit". That is a property of
    % the design, not the data, so it is refused up front. Rank, not columns.
    require_error_df(D.X);

    alpha = 0.05; E = 0.5; H = 2;
    fprintf('Effect: %s | statistic: %s | permutation: %s\n', ...
        D.effect_label, D.contrast_type, D.perm_type);

    %% Real statistic + permutation correction
    [stat, pm, ~, ~, zero_obs] = snpm_glm_stat(power, D.X, D.C);
    flctx = snpm_glm_fl_context(power, D);
    % Channels that missing data makes unevaluable in some permutations are
    % dropped from the observed map as well, so observed and permuted maps are
    % defined on one channel set. Computed once, before the permutation loop.
    % Zero-variance rule 1: channels with no variation left after the nuisance
    % model are dropped from the loop mask too (and reported with that reason).
    excluded_channels = excluded_channel_info(flctx.evaluable, chanlocs, preset);
    excluded_channels = snpm_excluded_add(excluded_channels, flctx.no_variation, ...
        reason_novar(), chanlocs);
    loop_ev = flctx.evaluable & ~flctx.no_variation;
    perm_stat_fn = @() permstat(flctx, D.X, D.C);
    snpm_perm_tap('record', ['loop:glm:main:' D.effect_label], []);
    [T, p, Clusters] = snpm_perm_correction(stat, pm, perm_stat_fn, neighbors, E, H, alpha, ...
        params.permutations, D.contrast_type, loop_ev, [], cs);
    % Zero-variance rule 3: an observed fit with a zero error term is reported
    % undefined (NaN) and listed; it stayed in every permuted map.
    [T, p, excluded_channels] = mark_undefined(T, p, zero_obs & loop_ev, ...
        excluded_channels, chanlocs);

    uncorrsigch      = find(p.real <= alpha);
    correctTFCEsigch = find(p.correctedTFCE <= alpha);
    sigclusters      = find([Clusters.p] <= alpha);
    SnPMsigch        = [Clusters(sigclusters).channels];

    %% Post-hoc pairwise maps (omnibus F presets only)
    % Each pair is its own two-level test on the rows of those two levels
    % (compute_posthoc / snpm_glm_pair), masked with the omnibus evaluable set.
    posthoc_results = empty_posthoc();
    posthoc_skipped = empty_skipped();
    if do_posthoc && strcmpi(D.contrast_type, 'F') && ~isempty(D.posthoc)
        [posthoc_results, posthoc_skipped] = compute_posthoc(power, meta, preset, opts, D, ...
            flctx.evaluable, neighbors, alpha, E, H, params.permutations, chanlocs, cs);
    end

    %% Output naming + plot
    [~, data_basename] = fileparts(data_file);
    base_filename = matlab.lang.makeValidName(sprintf('%s_%s', data_basename, preset));
    if ~exist(params.output_path, 'dir'), mkdir(params.output_path); end
    timestamp   = datestr(now, 'yyyymmdd_HHMMSS');
    outputSname = fullfile(params.output_path, [base_filename, '_', timestamp]);

    try
        % Clusters, not the flat SnPMsigch: the plotter numbers the clusters on
        % the Cluster image with the same order the report table uses. The run's
        % alpha rides along in the wrapper so both numberings key off one value;
        % the cell stops struct() expanding the struct array.
        TopoplotSignificant_single(T.real_T, uncorrsigch, correctTFCEsigch, ...
            struct('Clusters',{Clusters},'alpha',alpha), ...
            chanlocs, insidegoodch, preset, params.output_path, base_filename, select_mode);
    catch ME
        warning(ME.identifier, 'Could not generate significance topology plots: %s', ME.message);
    end

    %% Per-group mean topographies + standalone statistic map (omnibus F only)
    % For omnibus-F presets the report shows one mean topography per group plus
    % a standalone F-map (the two-means + diff layout used by t/r presets does
    % not apply). Group factor = the design's grouping column; levels come from
    % the same stable-unique logic snpm_glm_design's ref_dummy uses.
    group_labels    = {};
    group_mean_png  = {};
    fmap_png        = '';
    effect_map_only = false;   % regression: single slope map, no A/B mean cards
    if strcmpi(D.contrast_type, 'F')
        gv = [];
        if isfield(opts, 'group_col') && ~isempty(opts.group_col) ...
                && ismember(opts.group_col, meta.Properties.VariableNames)
            gv = meta.(opts.group_col);
        elseif isfield(opts, 'condition_col') && ~isempty(opts.condition_col) ...
                && ismember(opts.condition_col, meta.Properties.VariableNames)
            gv = meta.(opts.condition_col);
        end
        try
            if ~isempty(gv)
                [glevels, ~, gidx] = unique(gv, 'stable');
                for g = 1:numel(glevels)
                    rows = (gidx == g);
                    if ~any(rows), continue; end
                    gmean = mean(power(rows, :), 1, 'omitnan');
                    if all(isnan(gmean)), continue; end
                    lbl  = level_to_str(glevels, g);
                    safe = matlab.lang.makeValidName(lbl);
                    png  = [outputSname '_group-' safe '_mean.png'];
                    save_single_topo(gmean, chanlocs, sprintf('%s mean', lbl), png, 'hot');
                    [~, fn, ext] = fileparts(png);
                    group_labels{end+1}   = lbl;          %#ok<AGROW>
                    group_mean_png{end+1} = [fn ext];     %#ok<AGROW>
                end
            end
            % Standalone F-map from the per-channel F (T.real_T is the F path).
            fpng = [outputSname '_Fmap.png'];
            save_single_topo(T.real_T, chanlocs, 'F-map (omnibus)', fpng, 'hot');
            [~, fn, ext] = fileparts(fpng);
            fmap_png = [fn ext];
        catch ME
            warning(ME.identifier, 'Could not generate group / F-map topology plots: %s', ME.message);
        end
    else
        % t-contrast topographies. Two shapes:
        %  - 2-group factor (e.g. ANCOVA, 2 groups): condition-A/B mean topos +
        %    an effect map, via plot_topoInd -> <base>_Data1_topo.png /
        %    <base>_Data2_topo.png / <base>_topo.png (the names the report's
        %    t-branch expects), matching the legacy t-test layout.
        %  - no group factor (regression's continuous predictor): a single signed
        %    slope t-map only (A/B means are meaningless), -> <base>_topo.png,
        %    flagged so the report renders one card instead of three.
        try
            is_two_group = false;
            if isfield(opts, 'group_col') && ~isempty(opts.group_col) ...
                    && ismember(opts.group_col, meta.Properties.VariableNames)
                gv = meta.(opts.group_col);
                [glevels, ~, gidx] = unique(gv, 'stable');
                if numel(glevels) == 2
                    avgA = mean(power(gidx == 1, :), 1, 'omitnan');
                    avgB = mean(power(gidx == 2, :), 1, 'omitnan');
                    plot_topoInd(avgA, avgB, chanlocs, [base_filename '_Data1'], ...
                        [base_filename '_Data2'], params.output_path, base_filename);
                    is_two_group = true;
                end
            end
            if ~is_two_group
                tv = T.real_T; m = max(abs(tv(:)));
                if m == 0 || ~isfinite(m), m = 1; end
                save_single_topo(tv, chanlocs, [D.effect_label ' (t-map)'], ...
                    fullfile(params.output_path, [base_filename '_topo.png']), 'jet', [-m m]);
                effect_map_only = true;
            end
        catch ME
            warning(ME.identifier, 'Could not generate t-contrast topology plots: %s', ME.message);
        end
    end

    % Per-contrast significance maps for the post-hoc pairwise t-contrasts,
    % rendered exactly like the omnibus maps (black dots = uncorrected-sig,
    % white dots = correction-sig; one Cluster + one TFCE image that swap with
    % the report toggle). Reuses TopoplotSignificant_single; 'posthocT' gives a
    % signed diverging (jet) map. t = mean_a - mean_b, so warm = A > B.
    for ph = 1:numel(posthoc_results)
        safe = matlab.lang.makeValidName(posthoc_results(ph).label);
        map_base = [base_filename '_posthoc-' safe];
        try
            TopoplotSignificant_single(posthoc_results(ph).T.real_T, ...
                posthoc_results(ph).uncorrsigch, posthoc_results(ph).correctTFCEsigch, ...
                struct('Clusters',{posthoc_results(ph).Clusters},'alpha',alpha), ...
                chanlocs, insidegoodch, 'posthocT', ...
                params.output_path, map_base, select_mode);
            posthoc_results(ph).map_base = map_base;
        catch ME
            warning(ME.identifier, 'Could not generate post-hoc topology plots (%s): %s', ...
                posthoc_results(ph).label, ME.message);
        end
    end

    %% results_struct
    results_struct = struct();
    results_struct.T = T; results_struct.p = p; results_struct.Clusters = Clusters;
    results_struct.uncorrsigch = uncorrsigch;
    results_struct.correctTFCEsigch = correctTFCEsigch;
    results_struct.SnPMsigch = SnPMsigch;
    % The significance level the maps, clusters and TFCE sets above were built
    % at. Published so the report numbers its clusters with the SAME threshold
    % the run used, rather than a second hard-coded .05 that could drift.
    results_struct.alpha = alpha;
    % How each cluster was scored ('sum' | 'mean'; snpm_cluster_statistic).
    results_struct.cluster_statistic = cs;
    results_struct.chanlocs = chanlocs;
    results_struct.posthoc = posthoc_results;
    % Pairs that were NOT run, with the reason. A skipped pair has no map and no
    % entry in .posthoc, so without this list an output file could not tell a
    % reader that a comparison is absent; a console warning does not survive.
    results_struct.posthoc_skipped = posthoc_skipped;
    results_struct.excluded_channels = excluded_channels;
    results_struct.data_summary.data1_size = size(power);
    results_struct.data_summary.data2_size = [size(power,1), numel(channel_cols)];
    results_struct.data_summary.data1_mean = mean(power(:), 'omitnan');
    results_struct.data_summary.data2_mean = NaN;
    results_struct.data_summary.data1_std  = std(mean(power, 2, 'omitnan'), 'omitnan');
    results_struct.data_summary.data2_std  = NaN;
    results_struct.glm.preset = preset;
    results_struct.glm.effect_label = D.effect_label;
    results_struct.glm.contrast_type = D.contrast_type;
    results_struct.glm.effect_map_only = effect_map_only;  % regression -> single map
    results_struct.glm.perm_type = D.perm_type;
    results_struct.glm.n_obs = size(power, 1);
    % Observations, in the units the reader counts in. n_obs is ROWS; for a
    % repeated-measures design that is not the number of subjects, and for a
    % between-groups design "n = 32" hides whether it was 16 vs 16 or 28 vs 4.
    % Both are stated so the report's Observations row never has to infer one
    % from the other. group_sizes/group_size_labels are empty when the preset
    % has no grouping factor.
    [results_struct.glm.n_subjects, results_struct.glm.group_sizes, ...
        results_struct.glm.group_size_labels] = design_counts(meta, opts);
    results_struct.glm.group_labels   = group_labels;    % cellstr, parallel to ...
    results_struct.glm.group_mean_png = group_mean_png;  % ... mean-topo filenames
    results_struct.glm.fmap_png       = fmap_png;        % standalone F-map filename

    save([outputSname '.mat'], 'results_struct');

    try
        func_genSnpmTable(results_struct, uncorrsigch, correctTFCEsigch, SnPMsigch, chanlocs, outputSname, select_mode);
        snpm_write_settings_sheet(outputSname, cs, alpha);
        disp('Excel file generated successfully!');
    catch ME
        warning(ME.identifier, 'Could not generate Excel file: %s', ME.message);
    end

    %% Text
    results_text = {
        '=== SnPM GLM Group Analysis Results ===', '', ...
        ['Data file: ' data_file], ...
        ['Preset: ' preset], ...
        ['Effect: ' D.effect_label ' (' D.contrast_type ')'], ...
        ['Permutation: ' D.perm_type], ...
        ['Permutations: ' num2str(params.permutations)], ...
        cs_line, ...
        ['Observations: ' num2str(size(power,1))], ...
        excluded_text(excluded_channels), ...
        '', '--- Statistical Results ---', ...
        ['Uncorrected significant channels: ' num2str(numel(uncorrsigch))], ...
        ['TFCE corrected significant channels: ' num2str(numel(correctTFCEsigch))], ...
        ['Significant clusters N = : ' num2str(numel(sigclusters))], ...
        ['SnPM significant channels N = : ' num2str(numel(SnPMsigch))], ...
        '', ['Results saved to: ' outputSname '.mat'], ...
        ['Excel file saved to: ' outputSname '.xlsx'] }';
    results_text = [results_text; posthoc_text(posthoc_results, posthoc_skipped, '')];

    % params.tail is NOT set here. A GLM preset has no tail: an F contrast is a
    % non-negative omnibus statistic with no sidedness, and the 1-row t
    % contrasts are evaluated two-sided (|t| against a max-|t| null, see
    % snpm_perm_correction). This block used to park the contrast type ('t'/'F')
    % in params.tail purely to give generateAnalysisReport something to print,
    % which made the masthead of every GLM report read "one-tailed" -- the
    % report tests anything other than 'both' as one-sided. The contrast type
    % travels in its own field, results_struct.glm.contrast_type (set above and
    % read by stat_info), and the report now derives the tail wording from the
    % statistic instead of from an overloaded params.tail.
    try
        generateAnalysisReport(results_struct, params, base_filename, outputSname, uncorrsigch, correctTFCEsigch, SnPMsigch);
        results_text{end+1} = ['Report saved to: ' outputSname '_report.html'];
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

% ---------------------------------------------------------------------------
function [stat, pm] = permstat(flctx, X, C)
    Yp = snpm_glm_permute(flctx);
    [stat, pm] = snpm_glm_stat(Yp, X, C);
end

% ---- zero-variance rule (Method Spec 2026-10-14) -------------------------
% The reason strings are compared with the Python port; do not reword.
function s = reason_novar()
    s = 'no variation left after the nuisance model';
end

function s = reason_undefined()
    s = 'perfect fit, statistic undefined';
end

function s = reason_no_error_df()
    s = 'no error degrees of freedom left after fitting the design';
end

function tf = has_error_df(X)
% n_rows - rank(X) >= 1: at least one error degree of freedom.
    tf = size(X, 1) - rank(X) >= 1;
end

function require_error_df(X)
    if ~has_error_df(X)
        error('core_snpm_glm:noErrorDF', '%s', reason_no_error_df());
    end
end

function [T, p, excl] = mark_undefined(T, p, und, excl, chanlocs)
% Rule 3. `und` = channels in the loop mask whose OBSERVED fit hit rule 2
% (zero error term, statistic defined as 0). The 0 entered the observed TFCE
% integral and the channel stayed in every permuted map, which keeps the test
% exact; what is reported for it is NaN and an entry in the exclusion list,
% because a 0 there is a definition, not a measurement. It cannot be in an
% observed cluster: rule 2 also set its parametric p to 1.
    und = reshape(logical(und), 1, []);
    if ~any(und), return; end
    T.real_T(und)    = NaN;
    T.real_TFCE(und) = NaN;
    p.real(und)      = NaN;
    p.correctedTFCE(und) = NaN;
    T.excluded_channels = union(reshape(T.excluded_channels, 1, []), find(und));
    excl = snpm_excluded_add(excl, und, reason_undefined(), chanlocs);
end

function excl = loop_exclusions(finite_ev, no_var, und, chanlocs, context_label)
% One loop's exclusion report: missing data (the finiteness mask), then rule 1
% (no variation after the nuisance model), then rule 3 (perfect observed fit).
    excl = excluded_channel_info(finite_ev, chanlocs, context_label);
    excl = snpm_excluded_add(excl, no_var & finite_ev, reason_novar(), chanlocs);
    excl = snpm_excluded_add(excl, und, reason_undefined(), chanlocs);
end

function excl = excluded_channel_info(evaluable, chanlocs, context_label)
% Channels dropped because missing data makes them unevaluable in some
% permutations (see snpm_glm_fl_context). Reported loudly and returned so the
% caller can put them in the .mat / Excel / HTML — a channel must never leave an
% analysis silently.
    evaluable = reshape(logical(evaluable), 1, []);
    idx = find(~evaluable);
    labels = cell(1, numel(idx));
    for k = 1:numel(idx)
        lab = '';
        if idx(k) <= numel(chanlocs) && isfield(chanlocs, 'labels')
            lab = chanlocs(idx(k)).labels;
        end
        if isempty(lab), lab = sprintf('#%d', idx(k)); end
        labels{k} = lab;
    end
    excl = struct('n', numel(idx), 'n_channels', numel(evaluable), ...
        'index', idx, 'labels', {labels}, 'context', context_label, ...
        'reason', 'missing data: channel not evaluable in every permutation');
    if isempty(idx)
        return;
    end
    msg = sprintf(['%s: EXCLUDED %d of %d channels. Missing data means they cannot be ' ...
        'evaluated in every permutation, so they are NaN in the observed map as well as ' ...
        'in the null (a value scored against a null it could not enter would be invalid). ' ...
        'Channels: %s'], context_label, numel(idx), numel(evaluable), strjoin(labels, ', '));
    fprintf('%s\n', msg);
    warning('core_snpm_glm:channelsExcluded', '%s', msg);
end

function txt = excluded_text(excl)
% One results_text line per exclusion report (feeds the GUI Result box + CSV).
% Grouped by reason: a loop can exclude for missing data, for no variation
% left after the nuisance model, and for a perfect observed fit.
    if excl.n == 0
        txt = sprintf('Channels excluded: none (%d/%d evaluable)', ...
            excl.n_channels, excl.n_channels);
    else
        txt = sprintf('Channels excluded (%s): %d of %d — %s', ...
            excl.context, excl.n, excl.n_channels, excluded_by_reason(excl));
    end
end

function s = excluded_by_reason(excl)
% 'E1, E2 (reason A); E7 (reason B)'
    if isfield(excl, 'reasons') && numel(excl.reasons) == numel(excl.index)
        r = excl.reasons;
    else
        r = repmat({char(excl.reason)}, 1, numel(excl.index));
    end
    u = unique(r, 'stable');
    parts = cell(1, numel(u));
    for k = 1:numel(u)
        parts{k} = sprintf('%s (%s)', strjoin(excl.labels(strcmp(r, u{k})), ', '), u{k});
    end
    s = strjoin(parts, '; ');
end

function v = getfield_default(s, f, default)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = default; end
end

function guard_between_subject_rows(meta, opts, preset)
% Error if a between-subject ('free') preset is fed repeated-measures data (any
% subject appearing in more than one row). Uses the mapped subject column when
% present, otherwise any column whose name contains 'subject'. If no subject
% identifier is available the check is skipped (nothing to test against).
    subj_col = '';
    if isfield(opts, 'subject_col') && ~isempty(opts.subject_col) ...
            && ismember(opts.subject_col, meta.Properties.VariableNames)
        subj_col = opts.subject_col;
    else
        hit = find(strcmpi(meta.Properties.VariableNames, 'Subject'), 1);
        if isempty(hit)
            hit = find(contains(lower(meta.Properties.VariableNames), 'subject'), 1);
        end
        if ~isempty(hit), subj_col = meta.Properties.VariableNames{hit}; end
    end
    if isempty(subj_col), return; end   % no subject id -> cannot verify

    ids = meta.(subj_col);
    [uid, ~, gi] = unique(ids, 'stable');
    counts = accumarray(gi(:), 1);
    dup = find(counts > 1);
    if ~isempty(dup)
        error('core_snpm_glm:repeatedMeasures', ...
            ['Preset ''%s'' is a between-subject test (whole-subject permutation) and assumes ' ...
             'exactly one row per subject, but subject ''%s'' appears in %d rows (and %d other ' ...
             'subject(s) repeat). Repeated / trial-level measurements break exchangeability and ' ...
             'invalidate the permutation FWE control. Use ''rmanova'' for repeated conditions or ' ...
             '''mixedmodel'' for trial-level data.'], ...
            preset, level_to_str(uid, dup(1)), counts(dup(1)), numel(dup) - 1);
    end
end

function s = level_to_str(levels, i)
% Display label for the i-th level of a stable-unique group factor, robust to
% numeric / logical / categorical / string / cellstr columns.
    if iscell(levels), v = levels{i}; else, v = levels(i); end
    if isnumeric(v) || islogical(v)
        s = num2str(v);
    else
        s = char(string(v));
    end
end

% save_single_topo lives in dependencies/save_single_topo.m -- this file used to
% hold a byte-identical local copy, which shadowed it and had to be kept in sync
% by hand. Deleted; the shared function is on the path from the repo root.

% =========================================================================
% TWO-WAY MIXED ANOVA (mixed2way) — multi-effect pipeline
% =========================================================================
function [results_struct, results_text] = run_mixed2way(preset, meta, power, opts, ...
        chanlocs, neighbors, insidegoodch, select_mode, params, data_file, channel_cols)
% Computes the three tests of a two-way mixed ANOVA, each with its correct
% error term + permutation scheme, plus descriptive figures, a whole-head
% global test and interaction simple effects, and assembles one multi-effect
% results_struct that generateAnalysisReport renders as a combined report.
    alpha = 0.05; E = 0.5; H = 2;
    nperm = params.permutations;
    do_posthoc = getfield_default(params, 'posthoc', true);
    outpath = params.output_path;
    % validated and normalised by core_snpm_glm before the file was read
    [cs, cs_line] = snpm_cluster_statistic(params);

    for req = {'group_col','condition_col','subject_col'}
        if ~isfield(opts, req{1}) || isempty(opts.(req{1}))
            error('core_snpm_glm:mixed2wayCols', ...
                'mixed2way needs params.%s (group / condition / subject columns).', req{1});
        end
    end

    [~, data_basename] = fileparts(data_file);
    base_filename = matlab.lang.makeValidName(sprintf('%s_%s', data_basename, preset));
    if ~exist(outpath, 'dir'), mkdir(outpath); end
    timestamp   = datestr(now, 'yyyymmdd_HHMMSS');
    outputSname = fullfile(outpath, [base_filename '_' timestamp]);

    % ---- Balance (Method Spec B, approved 2026-10-14) ----
    % Every subject must have exactly one row per condition and one group.
    % Checked before anything is computed.
    require_mixed2way_balance(meta, opts);

    % ---- Error degrees of freedom (round-4 spec 3b) ----
    % Balanced, the group effect (subject means, anova1) has nS - nG error df
    % and the condition and interaction effects (full model) (nC-1)(nS - nG):
    % one of them is saturated exactly when all three are (one subject per
    % group). Nothing would be left to report, so the run is refused here,
    % with the main-analysis reason, rather than producing three skipped
    % effects.
    for eff_name = {'condition', 'interaction'}
        o = opts; o.effect = eff_name{1};
        Dchk = snpm_glm_design('mixed2way', meta, o);
        % the same design compute_effect is handed for this effect below
        snpm_glm_estimability(Dchk, 'enforce');
        require_error_df(Dchk.X);
    end

    % ---- ONE channel set for the whole report ----
    % A channel is analysed only if it is finite in EVERY row of the long
    % table. The same mask is applied to the group, condition and interaction
    % effects, the simple effects, every post-hoc, and the three whole-head
    % averages, so the excluded-channel list is identical for every effect.
    mask = all(isfinite(power), 1);

    fprintf('mixed2way: computing group, condition and interaction effects\n');

    % ---- (1) GROUP main effect: subject means, between-subjects ('free') ----
    snpm_perm_tap('record', 'scope:mixed2way:group', []);
    [eG, subj_power, subj_meta] = compute_group_effect(power, meta, opts, ...
        neighbors, alpha, E, H, nperm, do_posthoc, mask, chanlocs, cs);
    eG.key = 'group';
    eG = gen_effect_topos(eG, subj_power, subj_meta, opts.group_col, chanlocs, ...
        insidegoodch, select_mode, outpath, [base_filename '_group'], ...
        [outputSname '_group'], alpha);

    % ---- (2) CONDITION main effect: within-subject ----
    snpm_perm_tap('record', 'scope:mixed2way:condition', []);
    optsC = opts; optsC.effect = 'condition';
    Dc = snpm_glm_design('mixed2way', meta, optsC);
    eC = compute_effect(power, Dc, neighbors, alpha, E, H, nperm, do_posthoc, ...
        mask, meta, 'mixed2way', optsC, chanlocs, 'mixed2way condition', cs);
    eC.key = 'condition';
    eC = gen_effect_topos(eC, power, meta, opts.condition_col, chanlocs, ...
        insidegoodch, select_mode, outpath, [base_filename '_condition'], ...
        [outputSname '_condition'], alpha);

    % ---- (3) INTERACTION: within-subject. Map only here (the group x condition
    %          cell pattern is the descriptive grid); simple effects follow. ----
    snpm_perm_tap('record', 'scope:mixed2way:interaction', []);
    optsX = opts; optsX.effect = 'interaction';
    Dx = snpm_glm_design('mixed2way', meta, optsX);
    eX = compute_effect(power, Dx, neighbors, alpha, E, H, nperm, false, ...
        mask, meta, 'mixed2way', optsX, chanlocs, 'mixed2way interaction', cs);
    eX.key = 'interaction';
    eX = gen_effect_topos(eX, power, meta, '', chanlocs, ...
        insidegoodch, select_mode, outpath, [base_filename '_interaction'], ...
        [outputSname '_interaction'], alpha);
    [eX.simple_effects, eX.simple_skipped] = compute_simple_effects(power, meta, opts, neighbors, ...
        alpha, E, H, nperm, chanlocs, insidegoodch, select_mode, outpath, base_filename, outputSname, mask, cs);

    % ---- Global whole-head mixed ANOVA (parametric, channel-averaged) ----
    %
    % ONE AVERAGE SET, NOT THREE. All three effects are averaged over the same
    % channels: those finite in every row of the long table (the mask fixed at
    % the top of this function, which is also every effect's evaluable set).
    %
    % This comment used to say the opposite -- that each effect had its own
    % evaluable set and its own average, because the group effect was computed
    % on subject means taken over whatever rows a subject had
    % (mean(...,'omitnan')). That was the defect: a subject short of a
    % condition had its mean taken over fewer conditions, so a CONDITION effect
    % leaked into the GROUP test (with 7 of 20 subjects missing a condition the
    % group effect rejected 0.61 at nominal 0.05). The design is now required
    % to be balanced and the subject means are plain means over all conditions.
    %
    % Why one common set matters for a channel average at all: scalp power has
    % a strong spatial gradient, so the channel set fixes a per-unit offset.
    % Averaging units over different channel sets made a null paired contrast
    % reject 1.0000 of the time. Same rule as global_stat_test (legacy tier)
    % and core_snpm_lmm (LMM tier).
    %
    % The interaction SIMPLE EFFECTS (compute_simple_effects) deliberately have
    % no global test -- they emit maps only.
    %
    % ZERO-VARIANCE SPEC, GLM RULE 5: the average stays over the FINITENESS
    % mask (eff.finite_evaluable), not over each effect's loop mask. The loop
    % masks can differ by effect once rule 1 drops channels with no variation
    % left after that effect's nuisance model; averaging over them would give
    % the three whole-head tests three different channel sets again.
    Yg_G = global_average(subj_power, eG.finite_evaluable, 'mixed2way group');
    Yg_C = global_average(power,      eC.finite_evaluable, 'mixed2way condition');
    Yg_X = global_average(power,      eX.finite_evaluable, 'mixed2way interaction');
    [eG.global_stat, eG.global_pval, eG.global_df1, eG.global_df2, eG.global_reason] = compute_global_effect(Yg_G, eG.D);
    [eC.global_stat, eC.global_pval, eC.global_df1, eC.global_df2, eC.global_reason] = compute_global_effect(Yg_C, eC.D);
    [eX.global_stat, eX.global_pval, eX.global_df1, eX.global_df2, eX.global_reason] = compute_global_effect(Yg_X, eX.D);
    eG.global_n_channels = sum(eG.finite_evaluable);
    eC.global_n_channels = sum(eC.finite_evaluable);
    eX.global_n_channels = sum(eX.finite_evaluable);

    % ---- Descriptive figures: cell-mean grid + interaction line plot ----
    descriptive = struct('cellmean_png', {{}}, 'cellmean_labels', {{}}, ...
        'cell_rows', {{}}, 'cell_cols', {{}}, 'interaction_png', '', 'interaction_channels', []);
    try
        [cell_png, cell_lab, row_lab, col_lab] = plot_cellmean_grid(power, meta, ...
            opts.group_col, opts.condition_col, chanlocs, outputSname);
        descriptive.cellmean_png = cell_png; descriptive.cellmean_labels = cell_lab;
        descriptive.cell_rows = row_lab;     descriptive.cell_cols = col_lab;
    catch ME
        warning(ME.identifier, 'Could not build cell-mean grid: %s', ME.message);
    end
    ich = eX.SnPMsigch; if isempty(ich), ich = eX.correctTFCEsigch; end
    try
        [int_png, used_ch] = plot_interaction_lines(power, meta, opts.group_col, ...
            opts.condition_col, opts.subject_col, ich, [outputSname '_interaction_lines.png']);
        [~, fn, ext] = fileparts(int_png);
        descriptive.interaction_png = [fn ext];
        descriptive.interaction_channels = used_ch;
    catch ME
        warning(ME.identifier, 'Could not build interaction plot: %s', ME.message);
    end

    % ---- Assemble multi-effect results_struct (top level = interaction) ----
    effects = [normalize_effect(eG), normalize_effect(eC), normalize_effect(eX)];

    results_struct = struct();
    results_struct.T = eX.T; results_struct.p = eX.p; results_struct.Clusters = eX.Clusters;
    results_struct.uncorrsigch = eX.uncorrsigch;
    results_struct.correctTFCEsigch = eX.correctTFCEsigch;
    results_struct.SnPMsigch = eX.SnPMsigch;
    % The significance level the maps, clusters and TFCE sets above were built
    % at. Published so the report numbers its clusters with the SAME threshold
    % the run used, rather than a second hard-coded .05 that could drift.
    results_struct.alpha = alpha;
    % How each cluster was scored ('sum' | 'mean'), also on every effects(k)
    % and posthoc(k).
    results_struct.cluster_statistic = cs;
    results_struct.chanlocs = chanlocs;
    results_struct.posthoc = eX.posthoc;
    results_struct.posthoc_skipped = eX.posthoc_skipped;
    results_struct.excluded_channels = eX.excluded_channels;   % top level = interaction
    results_struct.data_summary.data1_size = size(power);
    results_struct.data_summary.data2_size = [size(power,1), numel(channel_cols)];
    results_struct.data_summary.data1_mean = mean(power(:), 'omitnan');
    results_struct.data_summary.data2_mean = NaN;
    results_struct.data_summary.data1_std  = std(mean(power, 2, 'omitnan'), 'omitnan');
    results_struct.data_summary.data2_std  = NaN;
    results_struct.glm.preset = preset;
    results_struct.glm.effect_label = eX.effect_label;
    results_struct.glm.contrast_type = eX.contrast_type;
    results_struct.glm.effect_map_only = false;
    results_struct.glm.perm_type = eX.perm_type;
    results_struct.glm.n_obs = size(power, 1);
    % Observations, in the units the reader counts in. n_obs is ROWS; for a
    % repeated-measures design that is not the number of subjects, and for a
    % between-groups design "n = 32" hides whether it was 16 vs 16 or 28 vs 4.
    % Both are stated so the report's Observations row never has to infer one
    % from the other. group_sizes/group_size_labels are empty when the preset
    % has no grouping factor.
    [results_struct.glm.n_subjects, results_struct.glm.group_sizes, ...
        results_struct.glm.group_size_labels] = design_counts(meta, opts);
    results_struct.glm.group_labels   = {};
    results_struct.glm.group_mean_png = {};
    results_struct.glm.fmap_png       = '';
    results_struct.effects = effects;
    results_struct.descriptive = descriptive;
    results_struct.global_stat = eX.global_stat;
    results_struct.global_pval = eX.global_pval;
    results_struct.global_reason = eX.global_reason;

    save([outputSname '.mat'], 'results_struct');

    try
        func_genSnpmTable(results_struct, eX.uncorrsigch, eX.correctTFCEsigch, ...
            eX.SnPMsigch, chanlocs, outputSname, select_mode);
        snpm_write_settings_sheet(outputSname, cs, alpha);
        disp('Excel file generated successfully!');
    catch ME
        warning(ME.identifier, 'Could not generate Excel file: %s', ME.message);
    end

    results_text = mixed2way_text(effects, data_file, preset, nperm, size(power,1), outputSname, cs_line);

    % No params.tail for mixed2way either -- see the note in the single-effect
    % path above. The three effects are reported with their own contrast types
    % in results_struct.effects(k).contrast_type.
    try
        generateAnalysisReport(results_struct, params, base_filename, outputSname, ...
            eX.uncorrsigch, eX.correctTFCEsigch, eX.SnPMsigch);
        results_text{end+1} = ['Report saved to: ' outputSname '_report.html'];
    catch ME
        warning(ME.identifier, 'Could not generate HTML report: %s', ME.message);
    end
    try
        export_results_csv(results_text, outputSname, params);
        results_text{end+1} = ['Result section CSV saved to: ' outputSname '_report.csv'];
    catch ME
        warning(ME.identifier, 'Could not write result CSV: %s', ME.message);
    end
end

% ---- one effect: real stat + FL permutation correction ------------------
function eff = compute_effect(Y, D, neighbors, alpha, E, H, nperm, do_posthoc, mask, meta, preset, opts, ...
        chanlocs, context_label, cs)
% mask  : 1 x nCh logical, the channel set this effect may test (fixed by the
%         caller before any permutation); ANDed with the design's own
%         Freedman-Lane evaluable set.
% meta, preset, opts : what D was built from, needed to rebuild each pairwise
%         post-hoc as its own two-level design (compute_posthoc).
% chanlocs, context_label : for this loop's own exclusion report.
% cs    : cluster statistic, 'sum' | 'mean' (snpm_cluster_statistic); also
%         handed to every post-hoc of this effect.
%
% Masks (zero-variance spec, GLM rules 1 and 4):
%   eff.finite_evaluable  finiteness mask (caller's mask AND the design's own
%                         Freedman-Lane finiteness). The whole-head average
%                         uses this one (rule 5), and so do the post-hocs as
%                         their starting mask.
%   eff.evaluable         the loop mask: finite_evaluable minus the channels
%                         with no variation left after the nuisance model.
    [stat, pm, df1, df2, zero_obs] = snpm_glm_stat(Y, D.X, D.C);
    flctx = snpm_glm_fl_context(Y, D);
    finite_ev = flctx.evaluable & reshape(logical(mask), 1, []);
    evaluable = finite_ev & ~flctx.no_variation;
    perm_fn = @() permstat(flctx, D.X, D.C);
    snpm_perm_tap('record', ['loop:glm:effect:' D.effect_label], []);
    [T, p, Clusters] = snpm_perm_correction(stat, pm, perm_fn, neighbors, E, H, alpha, ...
        nperm, D.contrast_type, evaluable, [], cs);
    excl = loop_exclusions(finite_ev, flctx.no_variation, false(size(finite_ev)), ...
        chanlocs, context_label);
    [T, p, excl] = mark_undefined(T, p, zero_obs & evaluable, excl, chanlocs);
    eff = struct();
    eff.evaluable = evaluable;
    eff.finite_evaluable = finite_ev;
    eff.excluded_channels = excl;
    eff.T = T; eff.p = p; eff.Clusters = Clusters;
    eff.uncorrsigch = find(p.real <= alpha);
    eff.correctTFCEsigch = find(p.correctedTFCE <= alpha);
    sigcl = find([Clusters.p] <= alpha);
    eff.SnPMsigch = [Clusters(sigcl).channels];
    eff.contrast_type = D.contrast_type;
    eff.effect_label = D.effect_label;
    eff.perm_type = D.perm_type;
    eff.df1 = df1; eff.df2 = df2;
    eff.cluster_statistic = cs;
    eff.D = D; eff.flctx = flctx;
    if do_posthoc && strcmpi(D.contrast_type, 'F') && ~isempty(D.posthoc)
        [eff.posthoc, eff.posthoc_skipped] = compute_posthoc(Y, meta, preset, opts, D, ...
            finite_ev, neighbors, alpha, E, H, nperm, chanlocs, cs);
    else
        eff.posthoc = empty_posthoc();
        eff.posthoc_skipped = empty_skipped();
    end
end

function ph = empty_posthoc()
% One entry per pairwise comparison that was RUN. Besides the maps it carries
% what the pair was run on, because a pair is tested on its own rows and that
% can be fewer subjects than the effect's header states:
%   .n_rows, .n_subjects      rows and subjects the pair used
%   .n_subjects_total         subjects in the effect's table (within designs)
%   .group_sizes              [n_a n_b] for between-subject pairs, [] otherwise
%   .dropped_subjects         cellstr: left out, only one of the two levels
%   .run_on                   the same in words, e.g. '23 of 25 subjects'
%   .excluded_channels        this pair's own exclusion report: missing data,
%                             no variation left after ITS nuisance model, and
%                             a perfect observed fit (zero-variance spec rule
%                             4: a pair can exclude a channel its omnibus
%                             tested)
%   .cluster_statistic        'sum' | 'mean', how this pair's clusters were scored
    ph = struct('label', {}, 'uncorrsigch', {}, 'correctTFCEsigch', {}, 'SnPMsigch', {}, ...
        'T', {}, 'p', {}, 'Clusters', {}, 'map_base', {}, 'n_rows', {}, 'n_subjects', {}, ...
        'n_subjects_total', {}, 'group_sizes', {}, 'dropped_subjects', {}, 'run_on', {}, ...
        'excluded_channels', {}, 'cluster_statistic', {});
end

function sk = empty_skipped()
% One entry per pairwise comparison that was NOT run: .label, .reason.
    sk = struct('label', {}, 'reason', {});
end

function lines = posthoc_text(posthoc, skipped, indent)
% Results-text lines for the post-hocs of one effect: what each was run on,
% who was left out, and every pair that was not run with its reason.
    lines = {};
    for k = 1:numel(posthoc)
        q = posthoc(k);
        s = sprintf('%sPost-hoc %s: %d TFCE-sig channels', indent, q.label, numel(q.correctTFCEsigch));
        if ~isempty(q.run_on)
            s = sprintf('%s. Tested on the rows of these two levels only: %s.', s, q.run_on);
        end
        if ~isempty(q.dropped_subjects)
            s = sprintf('%s Left out because they have only one of the two levels: %s.', s, ...
                strjoin(q.dropped_subjects, ', '));
        end
        if isfield(q, 'excluded_channels') && isstruct(q.excluded_channels) ...
                && q.excluded_channels.n > 0
            s = sprintf('%s Channels excluded for this pair: %s.', s, ...
                excluded_by_reason(q.excluded_channels));
        end
        lines{end+1, 1} = s; %#ok<AGROW>
    end
    for k = 1:numel(skipped)
        lines{end+1, 1} = sprintf('%sPost-hoc %s was not run: %s.', indent, ...
            skipped(k).label, skipped(k).reason); %#ok<AGROW>
    end
end

function [posthoc_results, skipped] = compute_posthoc(Y, meta, preset, opts, D, evaluable, neighbors, alpha, E, H, nperm, chanlocs, cs)
% Pairwise post-hoc t-maps for an omnibus-F effect (Method Spec A, approved
% 2026-10-14; rationale and rules in dependencies/snpm_glm_pair.m).
%
% Each "a vs b" is its OWN two-level test: the rows of levels a and b only,
% the preset's design rebuilt on those rows, Freedman-Lane with every column
% but the two-level dummy as nuisance, fresh permutation draws per pair. The
% recorded permutation index of a post-hoc loop therefore has the SUBSET's
% row count, not the omnibus row count.
%
% evaluable : the omnibus effect's FINITENESS mask, fixed before any
%             permutation. A post-hoc tests a channel only if it is in that set
%             AND finite in every row of the pair AND has variation left after
%             the pair's own nuisance model (zero-variance rule 1, applied per
%             loop). The omnibus's own zero-variance exclusions are not
%             inherited: rule 1 is a property of each loop's rows and nuisance
%             model, so a channel can be tested in the omnibus and excluded in
%             a pair (constant on that pair's rows), and each loop lists its own.
%
% A pair that cannot be tested is left out of posthoc_results, so that can
% hold fewer entries than D.posthoc; it is returned in `skipped` with the
% reason, and every output file states it (results text, CSV, workbook, HTML).
% Each pair that IS run carries what it was run on (rows, subjects, who was
% left out), for the same reason: the effect's header counts are not the
% post-hoc's.
    posthoc_results = empty_posthoc();
    skipped = empty_skipped();
    for ph = 1:numel(D.posthoc)
        lab = D.posthoc(ph).label;
        pr = snpm_glm_pair(preset, meta, opts, D.factor_col, D.posthoc(ph).a, D.posthoc(ph).b, lab);
        if ~pr.ok
            skipped(end+1) = struct('label', lab, 'reason', pr.reason); %#ok<AGROW>
            continue
        end
        Yp = Y(pr.rows, :);
        [phstat, phpm, ~, ~, zero_obs] = snpm_glm_stat(Yp, pr.D.X, pr.C);
        flph = snpm_glm_fl_context(Yp, pr.D);
        fin  = flph.evaluable & reshape(logical(evaluable), 1, []);
        ev   = fin & ~flph.no_variation;
        phfn = @() permstat(flph, pr.D.X, pr.C);
        snpm_perm_tap('record', ['loop:glm:posthoc:' lab], []);
        [Tph, pph, Cl_ph] = snpm_perm_correction(phstat, phpm, phfn, neighbors, E, H, alpha, ...
            nperm, 't', ev, [], cs);
        phex = loop_exclusions(fin, flph.no_variation, false(size(fin)), chanlocs, ...
            ['post-hoc ' lab]);
        [Tph, pph, phex] = mark_undefined(Tph, pph, zero_obs & ev, phex, chanlocs);
        k = numel(posthoc_results) + 1;
        posthoc_results(k).label = lab;
        posthoc_results(k).uncorrsigch = find(pph.real <= alpha);
        posthoc_results(k).correctTFCEsigch = find(pph.correctedTFCE <= alpha);
        posthoc_results(k).SnPMsigch = [Cl_ph(find([Cl_ph.p] <= alpha)).channels]; %#ok<FNDSB>
        % Full per-channel detail, so the report can present each pairwise
        % contrast as Cluster / TFCE / Uncorrected channel lists.
        posthoc_results(k).T = Tph;            % Tph.real_T = per-channel t, mean(a) - mean(b)
        posthoc_results(k).p = pph;            % pph.real, pph.correctedTFCE
        posthoc_results(k).Clusters = Cl_ph;   % cluster struct (.channels, .mass, .p)
        posthoc_results(k).map_base = '';
        posthoc_results(k).n_rows = pr.n_rows;
        posthoc_results(k).n_subjects = pr.n_subjects;
        posthoc_results(k).n_subjects_total = pr.n_subjects_total;
        posthoc_results(k).group_sizes = pr.group_sizes;
        posthoc_results(k).dropped_subjects = pr.dropped;
        posthoc_results(k).run_on = pr.run_on;
        posthoc_results(k).excluded_channels = phex;
        posthoc_results(k).cluster_statistic = cs;
        fprintf('  post-hoc %s: %d TFCE-sig channels (%s)\n', ...
            lab, numel(posthoc_results(k).correctTFCEsigch), pr.run_on);
    end
end

% ---- topoplots for one effect (significance maps + means/stat + post-hoc) --
function eff = gen_effect_topos(eff, Y, meta, level_col, chanlocs, insidegoodch, ...
        select_mode, out_path, map_base, outname_prefix, alpha)
% map_base    : base for the significance PNGs "<map_base> Cluster/TFCE.png".
% outname_prefix : prefix for mean / F / t map PNGs.
% level_col   : meta column for per-level mean topos (''=none -> stat map only).
% alpha       : the run's significance level, handed to the plotter so the
%               cluster numbers on the PNG match the report's cluster table.
    if strcmpi(eff.contrast_type, 'F'), topo_comp = 'mixed2way'; else, topo_comp = 'ancova'; end
    eff.map_base = map_base;
    try
        TopoplotSignificant_single(eff.T.real_T, eff.uncorrsigch, eff.correctTFCEsigch, ...
            struct('Clusters',{eff.Clusters},'alpha',alpha), ...
            chanlocs, insidegoodch, topo_comp, out_path, map_base, select_mode);
    catch ME
        warning(ME.identifier, 'Could not generate significance topology (%s): %s', map_base, ME.message);
    end

    eff.group_labels = {}; eff.group_mean_png = {}; eff.fmap_png = '';
    eff.data1_png = ''; eff.data2_png = ''; eff.diff_png = ''; eff.effect_map_only = false;

    have_levels = ~isempty(level_col) && ismember(level_col, meta.Properties.VariableNames);

    if strcmpi(eff.contrast_type, 'F')
        try
            if have_levels
                gv = meta.(level_col);
                [glevels, ~, gidx] = unique(gv, 'stable');
                for g = 1:numel(glevels)
                    rows = (gidx == g); if ~any(rows), continue; end
                    gmean = mean(Y(rows, :), 1, 'omitnan'); if all(isnan(gmean)), continue; end
                    lbl = level_to_str(glevels, g); safe = matlab.lang.makeValidName(lbl);
                    png = [outname_prefix '_level-' safe '_mean.png'];
                    save_single_topo(gmean, chanlocs, sprintf('%s mean', lbl), png, 'hot');
                    [~, fn, ext] = fileparts(png);
                    eff.group_labels{end+1} = lbl;          %#ok<AGROW>
                    eff.group_mean_png{end+1} = [fn ext];   %#ok<AGROW>
                end
            end
            fpng = [outname_prefix '_Fmap.png'];
            save_single_topo(eff.T.real_T, chanlocs, 'F-map', fpng, 'hot');
            [~, fn, ext] = fileparts(fpng); eff.fmap_png = [fn ext];
        catch ME
            warning(ME.identifier, 'Could not generate F / mean topology (%s): %s', map_base, ME.message);
        end
    else
        try
            is_two = false;
            if have_levels
                gv = meta.(level_col); [glevels, ~, gidx] = unique(gv, 'stable');
                if numel(glevels) == 2
                    avgA = mean(Y(gidx == 1, :), 1, 'omitnan');
                    avgB = mean(Y(gidx == 2, :), 1, 'omitnan');
                    plot_topoInd(avgA, avgB, chanlocs, [map_base '_Data1'], ...
                        [map_base '_Data2'], out_path, map_base);
                    eff.data1_png = [map_base '_Data1_topo.png'];
                    eff.data2_png = [map_base '_Data2_topo.png'];
                    eff.diff_png  = [map_base '_topo.png'];
                    eff.group_labels = {level_to_str(glevels,1), level_to_str(glevels,2)};
                    is_two = true;
                end
            end
            if ~is_two
                tv = eff.T.real_T; m = max(abs(tv(:))); if m == 0 || ~isfinite(m), m = 1; end
                save_single_topo(tv, chanlocs, [eff.effect_label ' (t-map)'], ...
                    fullfile(out_path, [map_base '_topo.png']), 'jet', [-m m]);
                eff.diff_png = [map_base '_topo.png'];
                eff.effect_map_only = true;
            end
        catch ME
            warning(ME.identifier, 'Could not generate t topology (%s): %s', map_base, ME.message);
        end
    end

    for ph = 1:numel(eff.posthoc)
        safe = matlab.lang.makeValidName(eff.posthoc(ph).label);
        pmb = [map_base '_posthoc-' safe];
        try
            TopoplotSignificant_single(eff.posthoc(ph).T.real_T, eff.posthoc(ph).uncorrsigch, ...
                eff.posthoc(ph).correctTFCEsigch, ...
                struct('Clusters',{eff.posthoc(ph).Clusters},'alpha',alpha), chanlocs, ...
                insidegoodch, 'posthocT', out_path, pmb, select_mode);
            eff.posthoc(ph).map_base = pmb;
        catch ME
            warning(ME.identifier, 'Could not generate post-hoc topology (%s): %s', ...
                eff.posthoc(ph).label, ME.message);
        end
    end
end

% ---- group main effect: subject means, between-subjects test -------------
function [eff, subj_power, subj_meta] = compute_group_effect(power, meta, opts, ...
        neighbors, alpha, E, H, nperm, do_posthoc, mask, chanlocs, cs)
% Subject means are PLAIN means over the subject's rows. run_mixed2way has
% already required one row per condition for every subject, so each mean is
% over the same conditions; a channel with any missing cell is outside `mask`
% and is NaN here as well (mean of a NaN), which is the intended result. The
% old mean(...,'omitnan') averaged a subject over whatever conditions it had
% and let a condition effect leak into the group test.
    [subjIDs, ~, si] = unique(meta.(opts.subject_col), 'stable');
    nS = numel(subjIDs); nCh = size(power, 2);
    subj_power = nan(nS, nCh);
    firstRow = zeros(nS, 1);
    for k = 1:nS
        rows = (si == k);
        subj_power(k, :) = mean(power(rows, :), 1);
        firstRow(k) = find(rows, 1);
    end
    gcolv = meta.(opts.group_col);
    subj_group = gcolv(firstRow);                    % subject's (constant) group, original type
    subj_meta  = table(subj_group, 'VariableNames', {opts.group_col});
    Dg = snpm_glm_design('anova1', subj_meta, struct('group_col', opts.group_col));
    snpm_glm_estimability(Dg, 'enforce');
    eff = compute_effect(subj_power, Dg, neighbors, alpha, E, H, nperm, do_posthoc, ...
        mask, subj_meta, 'anova1', struct('group_col', opts.group_col), chanlocs, 'mixed2way group', cs);
end

% ---- mixed2way balance requirement --------------------------------------
function require_mixed2way_balance(meta, opts)
% Every subject must have exactly one row per condition and belong to one
% group. A subject short of a condition would have its subject mean taken over
% fewer conditions, which biases the GROUP effect whenever conditions differ
% (measured: 7 of 20 subjects missing a condition, a condition effect present,
% no group effect -> the group test rejected 0.61 at nominal 0.05). A duplicated
% cell breaks the one-row-per-cell exchangeability the within-subject
% permutation assumes. Unbalanced data belong in the mixed model.
    [subj, ~, si] = unique(meta.(opts.subject_col), 'stable');
    [cond, ~, ci] = unique(meta.(opts.condition_col), 'stable');
    counts = accumarray([si(:), ci(:)], 1, [numel(subj), numel(cond)]);
    bad = find(any(counts ~= 1, 2));
    if ~isempty(bad)
        parts = cell(1, min(numel(bad), 10));
        for k = 1:numel(parts)
            s = bad(k);
            miss = find(counts(s, :) == 0);
            dup  = find(counts(s, :) > 1);
            bits = {};
            if ~isempty(miss)
                bits{end+1} = ['missing ' strjoin(arrayfun(@(j) level_to_str(cond, j), miss, ...
                    'UniformOutput', false), ', ')]; %#ok<AGROW>
            end
            if ~isempty(dup)
                bits{end+1} = ['duplicated ' strjoin(arrayfun(@(j) sprintf('%s (%d rows)', ...
                    level_to_str(cond, j), counts(s, j)), dup, 'UniformOutput', false), ', ')]; %#ok<AGROW>
            end
            parts{k} = sprintf('%s: %s', level_to_str(subj, s), strjoin(bits, '; '));
        end
        more = '';
        if numel(bad) > numel(parts), more = sprintf(' ... and %d more', numel(bad) - numel(parts)); end
        error('core_snpm_glm:mixed2wayUnbalanced', ...
            ['mixed2way needs every subject to have exactly one row per condition, and %d ' ...
             'subject(s) do not -- %s%s. A subject short of a condition would have its mean ' ...
             'taken over fewer conditions, which biases the group effect whenever the ' ...
             'conditions differ. Remove those subjects or complete their data; for ' ...
             'unbalanced or trial-level data use the mixed model (comparison ''mixedmodel'', ' ...
             'scripts/run_lmm_analysis.m).'], numel(bad), strjoin(parts, ' | '), more);
    end
    [~, ~, gi] = unique(meta.(opts.group_col), 'stable');
    ngrp = accumarray(si(:), gi(:), [numel(subj), 1], @(v) numel(unique(v)));
    multi = find(ngrp > 1);
    if ~isempty(multi)
        names = arrayfun(@(s) level_to_str(subj, s), multi(1:min(end, 10))', 'UniformOutput', false);
        error('core_snpm_glm:mixed2wayUnbalanced', ...
            ['mixed2way needs every subject to belong to exactly one group, and %d subject(s) ' ...
             'appear in more than one: %s. Group is the between-subject factor; if the same ' ...
             'people were measured under both "groups", that factor is within-subject and ' ...
             'this is not a mixed design -- use ''mixedmodel''.'], numel(multi), strjoin(names, ', '));
    end
end

% ---- interaction simple effects: condition within each group ------------
function [simples, skipped] = compute_simple_effects(power, meta, opts, neighbors, alpha, E, H, nperm, ...
        chanlocs, insidegoodch, select_mode, out_path, base_filename, outputSname, mask, cs)
% skipped : one entry (.label, .reason) per group whose simple effect was NOT
%           run, shaped like posthoc_skipped so every output file can say so.
%
% A simple effect is condition within ONE group, an rmanova on that group's
% rows. Its error term is the subject x condition residual, which needs at
% least 2 SUBJECTS. The guard used to count ROWS (sum(rows) < 2), so a group
% with one subject in 2+ conditions passed it, got an rmanova with zero error
% degrees of freedom, and was reported as "0 significant" -- a test that was
% never possible presented as a negative result. It now counts subjects and,
% below 2, skips the group with a reason that reaches the results text, the
% workbook and the HTML report.
    gcol = meta.(opts.group_col);
    [glev, ~, gi] = unique(gcol, 'stable');
    simples = normalize_effect(struct()); simples(:) = [];   % 0x0 typed struct
    skipped = empty_skipped();
    cnt = 0;
    for g = 1:numel(glev)
        rows = (gi == g);
        glab = sprintf('Condition within %s', level_to_str(glev, g));
        % Skip label and reason are compared with the Python port
        % (glm.plan_preset); do not reword.
        skip_label = sprintf('simple effect: condition within %s', level_to_str(glev, g));
        nsub = numel(unique(string(meta.(opts.subject_col)(rows))));
        if nsub < 2
            skipped(end+1) = struct('label', skip_label, 'reason', sprintf( ...
                'only %d subject in this group; a within-group test needs at least 2', ...
                nsub)); %#ok<AGROW>
            warning('core_snpm_glm:simpleEffectSkipped', 'Simple effect "%s" not run: %s.', ...
                skip_label, skipped(end).reason);
            continue;
        end
        % condition_levels: the WHOLE table's condition order, so each group's
        % simple effect lists and signs its pairs as the condition main effect
        % does even when that group's rows arrive in another order.
        sopts = struct('subject_col', opts.subject_col, 'condition_col', opts.condition_col);
        sopts.condition_levels = snpm_glm_levels(meta.(opts.condition_col));
        try
            Dg = snpm_glm_design('rmanova', meta(rows, :), sopts);
        catch ME
            warning(ME.identifier, 'simple-effect design failed (%s): %s', ...
                level_to_str(glev, g), ME.message);
            continue;
        end
        snpm_glm_estimability(Dg, 'enforce');
        % Round-4 spec 3b: a sub-loop with no error df is skipped with the
        % reason. Unreachable on a balanced design once the group has 2+
        % subjects ((n_g - 1)(nC - 1) >= 1), kept so the rule is local.
        if ~has_error_df(Dg.X)
            skipped(end+1) = struct('label', skip_label, 'reason', reason_no_error_df()); %#ok<AGROW>
            warning('core_snpm_glm:simpleEffectSkipped', 'Simple effect "%s" not run: %s.', ...
                skip_label, skipped(end).reason);
            continue;
        end
        snpm_perm_tap('record', ['scope:mixed2way:simple:' level_to_str(glev, g)], []);
        se = compute_effect(power(rows, :), Dg, neighbors, alpha, E, H, nperm, true, ...
            mask, meta(rows, :), 'rmanova', sopts, chanlocs, ...
            sprintf('mixed2way simple effect (%s)', level_to_str(glev, g)), cs);
        se.key = sprintf('condition@%s', level_to_str(glev, g));
        se.effect_label = glab;
        safe = matlab.lang.makeValidName(level_to_str(glev, g));
        se = gen_effect_topos(se, power(rows, :), meta(rows, :), opts.condition_col, ...
            chanlocs, insidegoodch, select_mode, out_path, ...
            [base_filename '_simple-' safe], [outputSname '_simple-' safe], alpha);
        cnt = cnt + 1;
        simples(cnt) = normalize_effect(se); %#ok<AGROW>
    end
end

% ---- observation counts for the report's model card ---------------------
function [nsub, gsizes, glabels] = design_counts(meta, opts)
% How many SUBJECTS the design has, and how they split across the grouping
% factor. Rows are not subjects in a repeated-measures design, and a total is
% not a split in a between-groups one; the report states both rather than
% guessing either. Everything is read from the meta table the design was built
% from, so it cannot disagree with the model that was fitted.
    nsub = NaN; gsizes = []; glabels = {};
    if ~istable(meta), return; end
    vn = meta.Properties.VariableNames;
    if isfield(opts,'subject_col') && ~isempty(opts.subject_col) && ismember(opts.subject_col, vn)
        nsub = numel(unique(string(meta.(opts.subject_col))));
    else
        nsub = height(meta);      % one row per subject when no subject column
    end
    if isfield(opts,'group_col') && ~isempty(opts.group_col) && ismember(opts.group_col, vn)
        gv = string(meta.(opts.group_col));
        if isfield(opts,'subject_col') && ~isempty(opts.subject_col) && ismember(opts.subject_col, vn)
            % Repeated measures: count each subject once, in the group of its
            % first row, so a 2-condition design does not report double the n.
            [~, first] = unique(string(meta.(opts.subject_col)), 'stable');
            gv = gv(first);
        end
        [lev, ~, li] = unique(gv, 'stable');
        gsizes  = accumarray(li(:), 1)';
        glabels = cellstr(lev(:)');
    end
end

% ---- global (whole-head, channel-averaged) parametric effect ------------
function [stat, pval, df1, df2, reason] = compute_global_effect(Yg, D)
% Parametric test on the channel average. Round-4 spec 3a: zero-variance rule
% 2 (statistic 0, p 1 on a zero error term) is a definition for permutation
% nulls only; a parametric p claims a place on the F/t distribution, and with
% a perfect fit the true statistic is infinite. So a perfectly fitted
% whole-head average is reported as undefined (NaN) with the GLM rule-3
% reason, and the report still shows the card with that reason.
    [stat, pval, df1, df2, zero_err] = snpm_glm_stat(Yg, D.X, D.C);
    reason = '';
    if any(zero_err)
        stat = NaN; pval = NaN;
        reason = reason_undefined();
    end
end

% ---- channel average over ONE effect's evaluable set --------------------
function Yg = global_average(Y, evaluable, label)
% Mean over the channels this effect can actually evaluate, so every analysed
% unit is averaged over the SAME channels (see the note at the call site).
% 'omitnan' is kept only as a belt-and-braces: inside the evaluable set every
% analysed row is finite by construction, so it never fires.
    evaluable = reshape(logical(evaluable), 1, []);
    if ~any(evaluable)
        warning('core_snpm_glm:noCommonChannel', ...
            ['%s global test: no channel is evaluable, so the whole-head statistic is ' ...
             'undefined (NaN). The channel-wise map for this effect is empty for the ' ...
             'same reason -- inspect the per-channel missingness.'], label);
        Yg = NaN(size(Y, 1), 1);
        return
    end
    Yg = mean(Y(:, evaluable), 2, 'omitnan');
end

% ---- canonical effect struct (fixed field set/order for concatenation) --
function e = normalize_effect(eff)
    e.key              = getfield_default(eff, 'key', '');
    e.effect_label     = getfield_default(eff, 'effect_label', '');
    e.contrast_type    = getfield_default(eff, 'contrast_type', 'F');
    e.perm_type        = getfield_default(eff, 'perm_type', '');
    e.T                = getfield_default(eff, 'T', []);
    e.p                = getfield_default(eff, 'p', []);
    e.Clusters         = getfield_default(eff, 'Clusters', []);
    e.uncorrsigch      = getfield_default(eff, 'uncorrsigch', []);
    e.correctTFCEsigch = getfield_default(eff, 'correctTFCEsigch', []);
    e.SnPMsigch        = getfield_default(eff, 'SnPMsigch', []);
    e.df1              = getfield_default(eff, 'df1', NaN);
    e.df2              = getfield_default(eff, 'df2', NaN);
    e.cluster_statistic = getfield_default(eff, 'cluster_statistic', 'sum');
    e.map_base         = getfield_default(eff, 'map_base', '');
    e.group_labels     = getfield_default(eff, 'group_labels', {});
    e.group_mean_png   = getfield_default(eff, 'group_mean_png', {});
    e.fmap_png         = getfield_default(eff, 'fmap_png', '');
    e.data1_png        = getfield_default(eff, 'data1_png', '');
    e.data2_png        = getfield_default(eff, 'data2_png', '');
    e.diff_png         = getfield_default(eff, 'diff_png', '');
    e.effect_map_only  = getfield_default(eff, 'effect_map_only', false);
    e.posthoc          = getfield_default(eff, 'posthoc', empty_posthoc());
    e.posthoc_skipped  = getfield_default(eff, 'posthoc_skipped', empty_skipped());
    e.simple_effects   = getfield_default(eff, 'simple_effects', struct([]));
    % Simple effects that were NOT run, with the reason (.label, .reason).
    e.simple_skipped   = getfield_default(eff, 'simple_skipped', empty_skipped());
    e.excluded_channels = getfield_default(eff, 'excluded_channels', ...
        struct('n', 0, 'n_channels', 0, 'index', [], 'labels', {{}}, ...
               'context', '', 'reason', ''));
    e.global_stat      = getfield_default(eff, 'global_stat', NaN);
    e.global_pval      = getfield_default(eff, 'global_pval', NaN);
    e.global_df1       = getfield_default(eff, 'global_df1', NaN);
    e.global_df2       = getfield_default(eff, 'global_df2', NaN);
    % Why the global is NaN when the engine declared it undefined (round-4
    % spec 3a: 'perfect fit, statistic undefined'); '' otherwise.
    e.global_reason    = getfield_default(eff, 'global_reason', '');
    % How many channels the global average was taken over. NOT cosmetic: the
    % global is only a whole-head claim when this is the whole montage. In a
    % mixed2way report every effect uses the one channel set fixed at the top
    % of run_mixed2way (finite in every row of the long table), so the three
    % globals rest on the SAME count; it is stored per effect because the
    % report reads it per effect.
    e.global_n_channels = getfield_default(eff, 'global_n_channels', NaN);
end

function txt = mixed2way_text(effects, data_file, preset, nperm, nobs, outputSname, cs_line)
    txt = {
        '=== SnPM Two-way Mixed ANOVA (group / condition / interaction) ===', '', ...
        ['Data file: ' data_file], ...
        ['Preset: ' preset], ...
        ['Permutations: ' num2str(nperm)], ...
        cs_line, ...
        ['Observations: ' num2str(nobs)], '', '--- Effects ---'}';
    for k = 1:numel(effects)
        e = effects(k);
        txt{end+1} = sprintf('%s (%s, %s): %d TFCE-sig, %d cluster-sig channels', ...
            e.key, e.contrast_type, e.perm_type, numel(e.correctTFCEsigch), numel(e.SnPMsigch)); %#ok<AGROW>
        if ~isnan(e.global_stat)
            if isnan(e.global_n_channels)
                kstr = '';
            else
                kstr = sprintf(' (mean over %d channels complete in every analysed unit)', ...
                    e.global_n_channels);
            end
            txt{end+1} = sprintf('    global %s = %.3f, p = %.4g%s', ...
                e.contrast_type, e.global_stat, e.global_pval, kstr); %#ok<AGROW>
        elseif ~isempty(e.global_reason)
            txt{end+1} = sprintf('    global %s: %s', e.contrast_type, e.global_reason); %#ok<AGROW>
        end
        if e.excluded_channels.n > 0
            txt{end+1} = ['    ' excluded_text(e.excluded_channels)]; %#ok<AGROW>
        end
        txt = [txt; posthoc_text(e.posthoc, e.posthoc_skipped, '    ')]; %#ok<AGROW>
        for s = 1:numel(e.simple_effects)
            se = e.simple_effects(s);
            txt{end+1} = sprintf('  %s: %d TFCE-sig, %d cluster-sig channels', ...
                se.effect_label, numel(se.correctTFCEsigch), numel(se.SnPMsigch)); %#ok<AGROW>
            if se.excluded_channels.n > 0
                txt{end+1} = ['    ' excluded_text(se.excluded_channels)]; %#ok<AGROW>
            end
            txt = [txt; posthoc_text(se.posthoc, se.posthoc_skipped, '    ')]; %#ok<AGROW>
        end
        for s = 1:numel(e.simple_skipped)
            txt{end+1} = sprintf('  %s was not run: %s.', e.simple_skipped(s).label, ...
                e.simple_skipped(s).reason); %#ok<AGROW>
        end
    end
    txt{end+1} = ['Results saved to: ' outputSname '.mat'];
    txt{end+1} = ['Excel file saved to: ' outputSname '.xlsx'];
end

