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

    D = snpm_glm_design(preset, meta, opts);

    % Between-subject presets (anova1/ancova/regression -> 'free' permutation)
    % assume ONE row per subject: whole-subject relabelling is only a valid null
    % if subjects are exchangeable. Trial-level / repeated data (a subject in
    % several rows) has within-subject correlation that breaks exchangeability
    % and silently loses family-wise-error control. Guard against it here.
    % rmanova/mixed2way (perm_type 'within', repeated measures) are exempt.
    if strcmp(D.perm_type, 'free')
        guard_between_subject_rows(meta, opts, preset);
    end

    alpha = 0.05; E = 0.5; H = 2;
    fprintf('Effect: %s | statistic: %s | permutation: %s\n', ...
        D.effect_label, D.contrast_type, D.perm_type);

    %% Real statistic + permutation correction
    [stat, pm] = snpm_glm_stat(power, D.X, D.C);
    flctx = make_fl_context(power, D);
    perm_stat_fn = @() permstat(flctx, D.X, D.C);
    [T, p, Clusters] = snpm_perm_correction(stat, pm, perm_stat_fn, neighbors, E, H, alpha, ...
        params.permutations, D.contrast_type);

    uncorrsigch      = find(p.real <= alpha);
    correctTFCEsigch = find(p.correctedTFCE <= alpha);
    sigclusters      = find([Clusters.p] <= alpha);
    SnPMsigch        = [Clusters(sigclusters).channels];

    %% Post-hoc pairwise maps (omnibus F presets only)
    posthoc_results = struct('label', {}, 'uncorrsigch', {}, 'correctTFCEsigch', {}, 'SnPMsigch', {}, ...
        'T', {}, 'p', {}, 'Clusters', {}, 'map_base', {});
    if do_posthoc && strcmpi(D.contrast_type, 'F') && ~isempty(D.posthoc)
        for ph = 1:numel(D.posthoc)
            Cph = pad_contrast(D.posthoc(ph).C, size(D.X, 2));
            [phstat, phpm] = snpm_glm_stat(power, D.X, Cph);
            % nuisance for a pairwise t = everything except this contrast's columns
            flph = flctx;   % same FL residualization (nuisance = original design nuisance)
            phfn = @() permstat(flph, D.X, Cph);
            [Tph, pph, Cl_ph] = snpm_perm_correction(phstat, phpm, phfn, neighbors, E, H, alpha, ...
                params.permutations, 't');
            posthoc_results(ph).label = D.posthoc(ph).label;
            posthoc_results(ph).uncorrsigch = find(pph.real <= alpha);
            posthoc_results(ph).correctTFCEsigch = find(pph.correctedTFCE <= alpha);
            posthoc_results(ph).SnPMsigch = [Cl_ph(find([Cl_ph.p] <= alpha)).channels]; %#ok<FNDSB>
            % Carry the full per-channel detail so the report can present each
            % pairwise contrast as Cluster / TFCE / Uncorrected channel lists
            % (statistic t, uncorrected p, TFCE p, cluster-level p), not counts.
            posthoc_results(ph).T = Tph;            % Tph.real_T = per-channel t
            posthoc_results(ph).p = pph;            % pph.real, pph.correctedTFCE
            posthoc_results(ph).Clusters = Cl_ph;   % cluster struct (.channels, .p)
            fprintf('  post-hoc %s: %d TFCE-sig channels\n', ...
                D.posthoc(ph).label, numel(posthoc_results(ph).correctTFCEsigch));
        end
    end

    %% Output naming + plot
    [~, data_basename] = fileparts(data_file);
    base_filename = matlab.lang.makeValidName(sprintf('%s_%s', data_basename, preset));
    if ~exist(params.output_path, 'dir'), mkdir(params.output_path); end
    timestamp   = datestr(now, 'yyyymmdd_HHMMSS');
    outputSname = fullfile(params.output_path, [base_filename, '_', timestamp]);

    try
        TopoplotSignificant_single(T.real_T, uncorrsigch, correctTFCEsigch, SnPMsigch, ...
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
                posthoc_results(ph).SnPMsigch, chanlocs, insidegoodch, 'posthocT', ...
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
    results_struct.chanlocs = chanlocs;
    results_struct.posthoc = posthoc_results;
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
    results_struct.glm.group_labels   = group_labels;    % cellstr, parallel to ...
    results_struct.glm.group_mean_png = group_mean_png;  % ... mean-topo filenames
    results_struct.glm.fmap_png       = fmap_png;        % standalone F-map filename

    save([outputSname '.mat'], 'results_struct');

    try
        func_genSnpmTable(results_struct, uncorrsigch, correctTFCEsigch, SnPMsigch, chanlocs, outputSname, select_mode);
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
        ['Observations: ' num2str(size(power,1))], ...
        '', '--- Statistical Results ---', ...
        ['Uncorrected significant channels: ' num2str(numel(uncorrsigch))], ...
        ['TFCE corrected significant channels: ' num2str(numel(correctTFCEsigch))], ...
        ['Significant clusters N = : ' num2str(numel(sigclusters))], ...
        ['SnPM significant channels N = : ' num2str(numel(SnPMsigch))], ...
        '', ['Results saved to: ' outputSname '.mat'], ...
        ['Excel file saved to: ' outputSname '.xlsx'] }';
    for ph = 1:numel(posthoc_results)
        results_text{end+1} = sprintf('Post-hoc %s: %d TFCE-sig channels', ...
            posthoc_results(ph).label, numel(posthoc_results(ph).correctTFCEsigch)); %#ok<AGROW>
    end

    % generateAnalysisReport reads params.tail (display only); GLM presets
    % have no tail, so provide a descriptive value.
    if ~isfield(params, 'tail') || isempty(params.tail)
        params.tail = D.contrast_type;
    end
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

function flctx = make_fl_context(Y, D)
    Z = D.X(:, D.nuisance_idx);
    if isempty(Z), Z = ones(size(Y,1), 1); end
    betaZ = Z \ Y;
    Zfit  = Z * betaZ;
    flctx = struct('Zfit', Zfit, 'R', Y - Zfit, 'eb', D.eb, 'perm_type', D.perm_type);
end

function Cpad = pad_contrast(C, p)
    Cpad = zeros(size(C,1), p);
    Cpad(:, 1:size(C,2)) = C;
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

function save_single_topo(data, chanlocs, title_str, out_png, cmap, clim)
% Single-map scalp topography saved to OUT_PNG. Mirrors create_topo_plot in
% dependencies/plot_topoInd.m (settings, sizing, dpi) so the per-group means
% and F-map match the look of the legacy difference maps.
%
% CLIM (optional) sets explicit colour limits. Signed maps (e.g. post-hoc
% pairwise t) must pass symmetric limits [-m m] so 0 sits at the colormap
% centre; without it 'minmax'/'auto' scaling is asymmetric and a mixed +/-
% map can read as all-one-colour. Omit CLIM (or pass []) to keep the legacy
% data-driven 'minmax' scaling used by the 'hot' F-map / group means.
    if nargin < 6, clim = []; end
    if isempty(clim), maplim = 'minmax'; else, maplim = clim; end
    plot_settings = {'headrad', 0.6, 'style', 'map', 'electrodes', 'on', ...
        'maplimits', maplim, 'whitebk', 'on'};
    fig = figure('Position', [50 50 450 450], 'Color', 'w', 'PaperPositionMode', 'auto');
    cl = onCleanup(@() close(fig)); %#ok<NASGU>
    topoplot(data, chanlocs, plot_settings{:});
    if ischar(cmap) || iscell(cmap) || isnumeric(cmap)
        if iscell(cmap), colormap(cmap{1}); else, colormap(cmap); end
    end
    if isempty(clim), caxis('auto'); else, caxis(clim); end
    axis off;
    set(gca, 'XLim', [-0.55 0.55], 'YLim', [-0.59 0.59]);
    delete(findobj(gca, 'Marker', '.'));
    title(title_str, 'Interpreter', 'none');
    colorbar;
    print(fig, '-dpng', '-r300', out_png);
end

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

    fprintf('mixed2way: computing group, condition and interaction effects\n');

    % ---- (1) GROUP main effect: subject means, between-subjects ('free') ----
    [eG, subj_power, subj_meta] = compute_group_effect(power, meta, opts, ...
        neighbors, alpha, E, H, nperm, do_posthoc);
    eG.key = 'group';
    eG = gen_effect_topos(eG, subj_power, subj_meta, opts.group_col, chanlocs, ...
        insidegoodch, select_mode, outpath, [base_filename '_group'], [outputSname '_group']);

    % ---- (2) CONDITION main effect: within-subject ----
    optsC = opts; optsC.effect = 'condition';
    Dc = snpm_glm_design('mixed2way', meta, optsC);
    eC = compute_effect(power, Dc, neighbors, alpha, E, H, nperm, do_posthoc);
    eC.key = 'condition';
    eC = gen_effect_topos(eC, power, meta, opts.condition_col, chanlocs, ...
        insidegoodch, select_mode, outpath, [base_filename '_condition'], [outputSname '_condition']);

    % ---- (3) INTERACTION: within-subject. Map only here (the group x condition
    %          cell pattern is the descriptive grid); simple effects follow. ----
    optsX = opts; optsX.effect = 'interaction';
    Dx = snpm_glm_design('mixed2way', meta, optsX);
    eX = compute_effect(power, Dx, neighbors, alpha, E, H, nperm, false);
    eX.key = 'interaction';
    eX = gen_effect_topos(eX, power, meta, '', chanlocs, ...
        insidegoodch, select_mode, outpath, [base_filename '_interaction'], [outputSname '_interaction']);
    eX.simple_effects = compute_simple_effects(power, meta, opts, neighbors, ...
        alpha, E, H, nperm, chanlocs, insidegoodch, select_mode, outpath, base_filename, outputSname);

    % ---- Global whole-head mixed ANOVA (parametric, channel-averaged) ----
    Yg      = mean(power, 2, 'omitnan');
    Yg_subj = mean(subj_power, 2, 'omitnan');
    [eG.global_stat, eG.global_pval, eG.global_df1, eG.global_df2] = compute_global_effect(Yg_subj, eG.D);
    [eC.global_stat, eC.global_pval, eC.global_df1, eC.global_df2] = compute_global_effect(Yg,      eC.D);
    [eX.global_stat, eX.global_pval, eX.global_df1, eX.global_df2] = compute_global_effect(Yg,      eX.D);

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
    results_struct.chanlocs = chanlocs;
    results_struct.posthoc = eX.posthoc;
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
    results_struct.glm.group_labels   = {};
    results_struct.glm.group_mean_png = {};
    results_struct.glm.fmap_png       = '';
    results_struct.effects = effects;
    results_struct.descriptive = descriptive;
    results_struct.global_stat = eX.global_stat;
    results_struct.global_pval = eX.global_pval;

    save([outputSname '.mat'], 'results_struct');

    try
        func_genSnpmTable(results_struct, eX.uncorrsigch, eX.correctTFCEsigch, ...
            eX.SnPMsigch, chanlocs, outputSname, select_mode);
        disp('Excel file generated successfully!');
    catch ME
        warning(ME.identifier, 'Could not generate Excel file: %s', ME.message);
    end

    results_text = mixed2way_text(effects, data_file, preset, nperm, size(power,1), outputSname);

    if ~isfield(params, 'tail') || isempty(params.tail), params.tail = 'F'; end
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
function eff = compute_effect(Y, D, neighbors, alpha, E, H, nperm, do_posthoc)
    [stat, pm, df1, df2] = snpm_glm_stat(Y, D.X, D.C);
    flctx = make_fl_context(Y, D);
    perm_fn = @() permstat(flctx, D.X, D.C);
    [T, p, Clusters] = snpm_perm_correction(stat, pm, perm_fn, neighbors, E, H, alpha, ...
        nperm, D.contrast_type);
    eff = struct();
    eff.T = T; eff.p = p; eff.Clusters = Clusters;
    eff.uncorrsigch = find(p.real <= alpha);
    eff.correctTFCEsigch = find(p.correctedTFCE <= alpha);
    sigcl = find([Clusters.p] <= alpha);
    eff.SnPMsigch = [Clusters(sigcl).channels];
    eff.contrast_type = D.contrast_type;
    eff.effect_label = D.effect_label;
    eff.perm_type = D.perm_type;
    eff.df1 = df1; eff.df2 = df2;
    eff.D = D; eff.flctx = flctx;
    if do_posthoc && strcmpi(D.contrast_type, 'F') && ~isempty(D.posthoc)
        eff.posthoc = compute_posthoc(Y, D, flctx, neighbors, alpha, E, H, nperm);
    else
        eff.posthoc = empty_posthoc();
    end
end

function ph = empty_posthoc()
    ph = struct('label', {}, 'uncorrsigch', {}, 'correctTFCEsigch', {}, 'SnPMsigch', {}, ...
        'T', {}, 'p', {}, 'Clusters', {}, 'map_base', {});
end

function posthoc_results = compute_posthoc(Y, D, flctx, neighbors, alpha, E, H, nperm)
% Pairwise post-hoc t-maps for an omnibus-F effect (same body as the legacy
% single-effect path). Nuisance = the omnibus design's nuisance (flctx reused).
    posthoc_results = empty_posthoc();
    for ph = 1:numel(D.posthoc)
        Cph = pad_contrast(D.posthoc(ph).C, size(D.X, 2));
        [phstat, phpm] = snpm_glm_stat(Y, D.X, Cph);
        phfn = @() permstat(flctx, D.X, Cph);
        [Tph, pph, Cl_ph] = snpm_perm_correction(phstat, phpm, phfn, neighbors, E, H, alpha, ...
            nperm, 't');
        posthoc_results(ph).label = D.posthoc(ph).label;
        posthoc_results(ph).uncorrsigch = find(pph.real <= alpha);
        posthoc_results(ph).correctTFCEsigch = find(pph.correctedTFCE <= alpha);
        posthoc_results(ph).SnPMsigch = [Cl_ph(find([Cl_ph.p] <= alpha)).channels]; %#ok<FNDSB>
        posthoc_results(ph).T = Tph;
        posthoc_results(ph).p = pph;
        posthoc_results(ph).Clusters = Cl_ph;
        posthoc_results(ph).map_base = '';
        fprintf('  post-hoc %s: %d TFCE-sig channels\n', ...
            D.posthoc(ph).label, numel(posthoc_results(ph).correctTFCEsigch));
    end
end

% ---- topoplots for one effect (significance maps + means/stat + post-hoc) --
function eff = gen_effect_topos(eff, Y, meta, level_col, chanlocs, insidegoodch, ...
        select_mode, out_path, map_base, outname_prefix)
% map_base    : base for the significance PNGs "<map_base> Cluster/TFCE.png".
% outname_prefix : prefix for mean / F / t map PNGs.
% level_col   : meta column for per-level mean topos (''=none -> stat map only).
    if strcmpi(eff.contrast_type, 'F'), topo_comp = 'mixed2way'; else, topo_comp = 'ancova'; end
    eff.map_base = map_base;
    try
        TopoplotSignificant_single(eff.T.real_T, eff.uncorrsigch, eff.correctTFCEsigch, ...
            eff.SnPMsigch, chanlocs, insidegoodch, topo_comp, out_path, map_base, select_mode);
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
                eff.posthoc(ph).correctTFCEsigch, eff.posthoc(ph).SnPMsigch, chanlocs, ...
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
        neighbors, alpha, E, H, nperm, do_posthoc)
    [subjIDs, ~, si] = unique(meta.(opts.subject_col), 'stable');
    nS = numel(subjIDs); nCh = size(power, 2);
    subj_power = nan(nS, nCh);
    firstRow = zeros(nS, 1);
    for k = 1:nS
        rows = (si == k);
        subj_power(k, :) = mean(power(rows, :), 1, 'omitnan');
        firstRow(k) = find(rows, 1);
    end
    gcolv = meta.(opts.group_col);
    subj_group = gcolv(firstRow);                    % subject's (constant) group, original type
    subj_meta  = table(subj_group, 'VariableNames', {opts.group_col});
    Dg = snpm_glm_design('anova1', subj_meta, struct('group_col', opts.group_col));
    eff = compute_effect(subj_power, Dg, neighbors, alpha, E, H, nperm, do_posthoc);
end

% ---- interaction simple effects: condition within each group ------------
function simples = compute_simple_effects(power, meta, opts, neighbors, alpha, E, H, nperm, ...
        chanlocs, insidegoodch, select_mode, out_path, base_filename, outputSname)
    gcol = meta.(opts.group_col);
    [glev, ~, gi] = unique(gcol, 'stable');
    simples = normalize_effect(struct()); simples(:) = [];   % 0x0 typed struct
    cnt = 0;
    for g = 1:numel(glev)
        rows = (gi == g);
        if sum(rows) < 2, continue; end
        try
            Dg = snpm_glm_design('rmanova', meta(rows, :), ...
                struct('subject_col', opts.subject_col, 'condition_col', opts.condition_col));
        catch ME
            warning(ME.identifier, 'simple-effect design failed (%s): %s', ...
                level_to_str(glev, g), ME.message);
            continue;
        end
        se = compute_effect(power(rows, :), Dg, neighbors, alpha, E, H, nperm, true);
        se.key = sprintf('condition@%s', level_to_str(glev, g));
        se.effect_label = sprintf('Condition within %s', level_to_str(glev, g));
        safe = matlab.lang.makeValidName(level_to_str(glev, g));
        se = gen_effect_topos(se, power(rows, :), meta(rows, :), opts.condition_col, ...
            chanlocs, insidegoodch, select_mode, out_path, ...
            [base_filename '_simple-' safe], [outputSname '_simple-' safe]);
        cnt = cnt + 1;
        simples(cnt) = normalize_effect(se); %#ok<AGROW>
    end
end

% ---- global (whole-head, channel-averaged) parametric effect ------------
function [stat, pval, df1, df2] = compute_global_effect(Yg, D)
    [stat, pval, df1, df2] = snpm_glm_stat(Yg, D.X, D.C);
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
    e.map_base         = getfield_default(eff, 'map_base', '');
    e.group_labels     = getfield_default(eff, 'group_labels', {});
    e.group_mean_png   = getfield_default(eff, 'group_mean_png', {});
    e.fmap_png         = getfield_default(eff, 'fmap_png', '');
    e.data1_png        = getfield_default(eff, 'data1_png', '');
    e.data2_png        = getfield_default(eff, 'data2_png', '');
    e.diff_png         = getfield_default(eff, 'diff_png', '');
    e.effect_map_only  = getfield_default(eff, 'effect_map_only', false);
    e.posthoc          = getfield_default(eff, 'posthoc', empty_posthoc());
    e.simple_effects   = getfield_default(eff, 'simple_effects', struct([]));
    e.global_stat      = getfield_default(eff, 'global_stat', NaN);
    e.global_pval      = getfield_default(eff, 'global_pval', NaN);
    e.global_df1       = getfield_default(eff, 'global_df1', NaN);
    e.global_df2       = getfield_default(eff, 'global_df2', NaN);
end

function txt = mixed2way_text(effects, data_file, preset, nperm, nobs, outputSname)
    txt = {
        '=== SnPM Two-way Mixed ANOVA (group / condition / interaction) ===', '', ...
        ['Data file: ' data_file], ...
        ['Preset: ' preset], ...
        ['Permutations: ' num2str(nperm)], ...
        ['Observations: ' num2str(nobs)], '', '--- Effects ---'}';
    for k = 1:numel(effects)
        e = effects(k);
        txt{end+1} = sprintf('%s (%s, %s): %d TFCE-sig, %d cluster-sig channels', ...
            e.key, e.contrast_type, e.perm_type, numel(e.correctTFCEsigch), numel(e.SnPMsigch)); %#ok<AGROW>
        if ~isnan(e.global_stat)
            txt{end+1} = sprintf('    global %s = %.3f, p = %.4g', ...
                e.contrast_type, e.global_stat, e.global_pval); %#ok<AGROW>
        end
    end
    txt{end+1} = ['Results saved to: ' outputSname '.mat'];
    txt{end+1} = ['Excel file saved to: ' outputSname '.xlsx'];
end

