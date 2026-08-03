function [results_struct, results_text] = core_snpm_circ(params)
%CORE_SNPM_CIRC  Tier-2 circular (phase) analyses with TFCE + cluster correction.
%
% Routed from core_snpm_analysis by an early guard, exactly like core_snpm_glm
% and core_snpm_lmm. Supported params.comparison:
%
%   'circ_phase_group'     PRIMARY. Covariate-adjusted Hotelling T^2 on the
%                          (cos,sin) embedding, two independent groups.
%   'circ_phase_group_u2'  SECONDARY. Watson's U^2, same design. Omnibus over
%                          the whole distribution, but it CANNOT adjust for
%                          event count.
%   'circ_corrAngLinear'   Circular-linear association between the angle and a
%                          linear measure, paired by subject.
%
% (circ_wheeler_watson_Test and circ_WatsonsU2Test are RETIRED, together with
% the two _circ engines that implemented them.)
%
% ------------------------------------------------------------------------
% THE FOUR THINGS THAT MAKE THIS CORRECT RATHER THAN MERELY RUNNABLE
%
% 1. THE ANGLE MATRIX IS NEVER PERMUTED. Freedman-Lane permutes the RESPONSE
%    through snpm_glm_permute, which applies ONE permutation index to every
%    column of a stacked [cos, sin] matrix. Shuffling angles invites the bug
%    where cos and sin move with different indices and the embedding silently
%    desynchronises; passing the response makes that unrepresentable. Under FL
%    the permuted response leaves the unit circle (only the residuals move, the
%    nuisance fit is held fixed), which is why the permuted statistic is
%    evaluated by snpm_circ_hotelling_perm rather than snpm_circ_hotelling.
%
% 2. THE EVENT-COUNT COVARIATE IS PER CHANNEL, not a per-subject scalar. Event
%    count varies strongly across the montage and the confound acts channel-wise,
%    so a per-subject scalar under-adjusts exactly at the posterior channels
%    where the problem is worst. Each channel therefore gets its own
%    Freedman-Lane nuisance fit, built ONCE before the permutation loop.
%    params.circ_count_covariate_mode = 'per_channel' (default) | 'per_subject'.
%
% 3. THE ESTIMABILITY GATE BELONGS HERE, NOT TO THE STATISTIC. min_per_group
%    defaults to OFF inside snpm_circ_hotelling / snpm_circ_watsons_u2 on
%    purpose: a non-zero default would be recomputed on every permuted call, so
%    the analysed channel set would vary across the null and the family-wise
%    correction would no longer be calibrated. The >=8-per-group mask is
%    computed ONCE from the observed labels here and held FIXED, and that is
%    asserted at run time (see assert_gate_fixed below).
%
% 4. THERE IS NO GLOBAL / OMNIBUS TEST, deliberately. global_stat_test averages
%    channels arithmetically, which is meaningless for angles; and a circularly
%    correct whole-head mean would average over a real 25-50 deg
%    anterior-posterior phase gradient. results_struct.global_stat/global_pval
%    are never set, so generateAnalysisReport simply omits the Global test card.
%    A DESCRIPTIVE panel (per-group mean direction, circular SD, Rayleigh) is
%    emitted instead and is labelled no-inference.
%
% Cluster statistic is cluster MASS (snpm_perm_correction's), which differs from
% Helfrich's FieldTrip 'maxsize'; disclosed in the results text.
%
% INPUT params (beyond the usual data/output/channels/permutations):
%   .circ_units             REQUIRED 'rad'|'deg'          (no default)
%   .circ_convention        REQUIRED, see snpm_circ_import (no default)
%   .circ_zero_offset_deg   required for convention 'custom'
%   .count1_file/.count2_file  REQUIRED for the two phase-group analyses
%   .circ_count_covariate_mode 'per_channel' (default) | 'per_subject'
%   .circ_collinearity_max  |point-biserial| above which log(count) is treated
%                           as collinear with group. DEFAULT 0.8.
%   .measure_col            circ_corrAngLinear: name of the linear measure column
%
% See also SNPM_CIRC_HOTELLING, SNPM_CIRC_HOTELLING_PERM, SNPM_CIRC_WATSONS_U2,
%          SNPM_CIRC_CORR_ANGLINEAR, SNPM_CIRC_IMPORT, SNPM_PERM_CORRECTION.

    comparison = params.comparison;
    alpha = 0.05; E = 0.5; H = 2;
    nperm = params.permutations;
    is_u2   = strcmp(comparison, 'circ_phase_group_u2');
    is_grp  = ismember(comparison, {'circ_phase_group', 'circ_phase_group_u2'});
    is_corr = strcmp(comparison, 'circ_corrAngLinear');
    if ~(is_grp || is_corr)
        error('core_snpm:circUnknownComparison', ...
            'core_snpm_circ does not handle comparison ''%s''.', comparison);
    end

    % Tails: the circular statistics are non-negative omnibus quantities. A
    % one-sided request is not a weaker version of them, it is meaningless, so
    % it is refused rather than quietly ignored.
    if isfield(params, 'tail') && ~isempty(params.tail) && ~strcmpi(params.tail, 'both')
        error('core_snpm:circTailNotSupported', ...
            ['params.tail = ''%s'' is not defined for a circular analysis. The ' ...
             'Hotelling F, Watson U^2 and circular-linear F are all non-negative ' ...
             'omnibus statistics with a single upper tail; use tail = ''both''.'], ...
            params.tail);
    end
    params.tail = 'both';
    if isfield(params, 'datatype') && ~isempty(params.datatype) && ...
            ~strcmpi(params.datatype, 'absolute')
        error('core_snpm:circDatatypeNotSupported', ...
            ['params.datatype = ''%s'' is undefined on a circular domain (log10 of ' ...
             'a radian goes complex; a cross-channel z-score of an angle is ' ...
             'meaningless). Use ''absolute''.'], params.datatype);
    end
    params.datatype = 'absolute';

    if ~exist(params.output_path, 'dir'), mkdir(params.output_path); end

    %% ---- load angles -------------------------------------------------
    D1 = read_wide(params.data1_file, getf(params, 'data1_sheet', 'CSV File'));
    D2 = read_wide(params.data2_file, getf(params, 'data2_sheet', 'CSV File'));

    res_gate = 'warn';
    if is_u2, res_gate = 'error'; end
    iopt = struct('units', getf(params, 'circ_units', ''), ...
        'convention', getf(params, 'circ_convention', ''), ...
        'zero_offset_deg', getf(params, 'circ_zero_offset_deg', []), ...
        'resolution_gate', res_gate);

    iopt.label = 'Angles file 1';
    [A, imp1] = snpm_circ_import(D1.M, iopt);
    yv = [];

    if is_grp
        % Group B angles must live on the same channels; reorder by name.
        D2 = align_columns(D2, D1.cols, 'Group B angle file');
        iopt.label = 'Angles file 2';
        A2 = snpm_circ_import(D2.M, iopt);
        n1g = size(A,1); n2g = size(A2,1);
        A  = [A; A2];
        gi = [ones(n1g,1); zeros(n2g,1)];   % 1 = group A
    else
        % Circular-linear: match subjects between the angle file and the
        % single-column linear measure file.
        [A, yv] = match_measure(A, D1, D2, params);
        gi = [];
        n1g = size(A,1); n2g = 0;
    end
    nSubj = size(A, 1);
    fprintf(['Circular import: units=%s convention=%s -> rotation %+.1f deg applied; ' ...
        'angular resolution %.4g deg\n'], imp1.units, imp1.convention, ...
        imp1.rotation_deg, imp1.resolution_deg);

    %% ---- channel setup (shared label-based montage selection) --------
    chsel        = snpm_setup_channels(params.channels, {A}, D1.cols);
    A            = chsel.data{1};
    chanlocs     = chsel.chanlocs;
    neighbors    = chsel.neighbors;
    insidegoodch = chsel.insidegoodch;
    select_mode  = chsel.select_mode;
    keptcols     = chsel.channel_mapping.original_data_columns;
    nCh          = size(A, 2);
    fprintf('%s montage: %d/%d data columns matched to channel locations\n', ...
        params.channels, nCh, numel(D1.cols));
    if size(neighbors, 1) ~= nCh
        error('core_snpm:circChannelMismatch', ...
            'Selected %d channels but the neighbour matrix has %d rows.', ...
            nCh, size(neighbors, 1));
    end

    %% ---- event counts -> per-channel log covariate --------------------
    logC = []; usecov = false(1, nCh); collinear = false(1, nCh);
    cov_report = struct('mode', 'none', 'n_dropped_constant', 0, ...
        'n_collinear', 0, 'pb_max', NaN, 'pb_median', NaN, ...
        'n_cells_no_count', 0, 'frac_cells_no_count', 0);
    if is_grp
        Cnt = load_counts(params, D1, D2, keptcols);
        [A, logC, usecov, collinear, cov_report] = snpm_circ_count_covariate( ...
            A, gi, Cnt, params, chanlocs);
    end

    %% ---- convention / TurtleWave-inversion check ----------------------
    % Pooled over BOTH groups, so it is label-invariant and cannot contaminate
    % the inference; restricted to frontal channels that pass the gate.
    gate0 = gate_mask(A, gi, is_grp, 8) & ~collinear;
    snpm_circ_check_inversion(A, gate0, chanlocs, imp1);

    %% ---- estimability gate (computed ONCE, held FIXED) ----------------
    gate = gate0;
    if is_grp
        [~, ~, ~, ~, exg] = snpm_circ_hotelling(A, gi);      % the documented idiom
        gate = gate & (exg.n1 >= 8) & (exg.n2 >= 8);
    else
        gate = gate & (sum(isfinite(A) & isfinite(yv), 1) >= 8);
    end
    % Channels using a per-channel covariate must be COMPLETE: Freedman-Lane
    % fits the nuisance over whole columns (same rule as snpm_glm_fl_context).
    complete = all(isfinite(A), 1);
    evaluable = gate & (~usecov | complete);
    if ~any(evaluable)
        error('core_snpm:circNoEvaluableChannels', ...
            ['No channel is estimable: every channel fails the >=8-per-group gate, ' ...
             'the completeness requirement, or the log(count) collinearity check. ' ...
             '%d/%d channels were collinear with group.'], sum(collinear), nCh);
    end
    excluded_channels = excluded_info(evaluable, chanlocs, comparison);
    fprintf('Estimability: %d/%d channels analysed (%d excluded)\n', ...
        sum(evaluable), nCh, excluded_channels.n);

    %% ---- observed statistic + permutation correction ------------------
    tfce_dh = [];   % [] = ClusterEnhancement's own default (0.1); see U^2 below
    if is_grp && ~is_u2
        hctx = snpm_circ_hotelling_perm(gi, logC, usecov, evaluable);
        [stat, pm] = hctx.eval(cos(A), sin(A));
        flctx = snpm_circ_fl_context(A, logC, usecov, evaluable);
        perm_fn = make_hotelling_perm(flctx, hctx, nCh, evaluable, stat);
    elseif is_u2
        % Watson U^2 has no covariate slot: the count confound is NOT adjusted.
        if any(usecov)
            warning('core_snpm:circU2Unadjusted', ...
                ['Watson U^2 cannot take the event-count covariate, so this map is ' ...
                 'UNADJUSTED for event count (max |point-biserial(log count, group)| ' ...
                 '= %.2f). Use comparison = ''circ_phase_group'' if that association ' ...
                 'is material.'], cov_report.pb_max);
        end
        sopt = struct('valid_channels', evaluable);
        [stat, pm] = snpm_circ_watsons_u2(A, gi, sopt);
        report_u2_ties(A, evaluable);
        % TFCE integration step. U^2 lives on ~[0, 0.5], so the default
        % dh = 0.1 leaves ~5 integration levels and enhances 50-79 of 178
        % channels to exactly 0, which makes the integral meaningless. Pass
        % dh = 0.005 (~91 levels) straight through to ClusterEnhancement; it is
        % applied identically to the observed and to every permuted map.
        tfce_dh = 0.005;
        perm_fn = make_u2_perm(A, gi, sopt, evaluable, stat);
    else
        sopt = struct('valid_channels', evaluable);
        [stat, pm, ~, ~, rho] = snpm_circ_corr_anglinear(yv, A, [], sopt);
        flctx = struct('Zfit', repmat(mean(yv, 'omitnan'), nSubj, 1), ...
            'R', yv - mean(yv, 'omitnan'), 'eb', ones(nSubj,1), 'perm_type', 'free');
        perm_fn = make_corr_perm(flctx, A, sopt, evaluable, stat);
    end

    [T, p, Clusters] = snpm_perm_correction(stat, pm, perm_fn, neighbors, E, H, ...
        alpha, nperm, 'F', evaluable, tfce_dh);

    uncorrsigch      = find(p.real <= alpha);
    correctTFCEsigch = find(p.correctedTFCE <= alpha);
    sigclusters      = find([Clusters.p] <= alpha);
    SnPMsigch        = [Clusters(sigclusters).channels];

    %% ---- descriptive circular panel (NO inference) --------------------
    desc = circ_descriptive(A, gi, is_grp, evaluable);

    %% ---- output naming + figures --------------------------------------
    [~, b1] = fileparts(params.data1_file);
    [~, b2] = fileparts(params.data2_file);
    b1 = sprintf('%s_%s', b1, getf(params, 'data1_sheet', 'CSV File'));
    b2 = sprintf('%s_%s', b2, getf(params, 'data2_sheet', 'CSV File'));
    base_filename = sprintf('%s VS %s', b1, b2);
    timestamp   = datestr(now, 'yyyymmdd_HHMMSS');
    outputSname = fullfile(params.output_path, [base_filename, '_', timestamp]);

    if is_grp
        plot_circ_topos(A, gi, chanlocs, b1, b2, params.output_path, base_filename);
    else
        % Circular-linear: one mean-direction map + the statistic map. The
        % report's hide_condition_b layout asks for exactly these two files.
        try
            plot_topoInd(circ_mean_deg(A), [], chanlocs, b1, b2, params.output_path, ...
                base_filename, 'cmap1', 'hsv', 'clim1', [0 360], 'skip2', true, ...
                'diffdata', T.real_T, 'cmap_diff', 'hot', ...
                'difftitle', 'Circular-linear F');
        catch ME
            warning(ME.identifier, 'Could not generate circular-linear topology plots: %s', ME.message);
        end
    end
    try
        TopoplotSignificant_single(T.real_T, uncorrsigch, correctTFCEsigch, SnPMsigch, ...
            chanlocs, insidegoodch, comparison, params.output_path, base_filename, select_mode);
    catch ME
        warning(ME.identifier, 'Could not generate significance topology plots: %s', ME.message);
    end

    %% ---- results struct (same contract the table/report/topo chain reads)
    results_struct = struct();
    results_struct.T = T; results_struct.p = p; results_struct.Clusters = Clusters;
    results_struct.uncorrsigch = uncorrsigch;
    results_struct.correctTFCEsigch = correctTFCEsigch;
    results_struct.SnPMsigch = SnPMsigch;
    results_struct.chanlocs = chanlocs;
    results_struct.excluded_channels = excluded_channels;
    results_struct.hide_condition_b  = is_corr;   % no "group B" topography
    results_struct.per_channel_n = sum(isfinite(A), 1);
    results_struct.circ = struct('comparison', comparison, 'import', imp1, ...
        'rotation_deg', imp1.rotation_deg, 'convention', imp1.convention, ...
        'units', imp1.units, 'covariate', cov_report, 'gate_min_per_group', 8, ...
        'evaluable', evaluable, 'cluster_statistic', 'mass (sum stat / n channels)', ...
        'n_group1', n1g, 'n_group2', n2g);
    if is_corr, results_struct.circ.rho = rho; end
    results_struct.circ_descriptive = desc;
    % Recorded in the report header via the existing effect-label slot. NOTE:
    % glm.contrast_type is deliberately NOT set, so stat_info keeps resolving
    % the statistic from the comparison key ('circ').
    results_struct.glm.effect_label = sprintf( ...
        '%s | zero convention %s (rotation %+.1f deg) | count covariate: %s', ...
        stat_label(comparison), imp1.convention, imp1.rotation_deg, cov_report.mode);
    results_struct.data_summary.data1_size = [n1g, nCh];
    results_struct.data_summary.data2_size = [max(n2g,0), nCh];
    results_struct.data_summary.data1_mean = desc.group(1).mean_deg;
    if numel(desc.group) > 1
        results_struct.data_summary.data2_mean = desc.group(2).mean_deg;
        results_struct.data_summary.data2_std  = desc.group(2).sd_deg;
    else
        results_struct.data_summary.data2_mean = NaN;
        results_struct.data_summary.data2_std  = NaN;
    end
    results_struct.data_summary.data1_std = desc.group(1).sd_deg;

    save([outputSname '.mat'], 'results_struct');
    try
        func_genSnpmTable(results_struct, uncorrsigch, correctTFCEsigch, SnPMsigch, ...
            chanlocs, outputSname, select_mode);
        disp('Excel file generated successfully!');
    catch ME
        warning(ME.identifier, 'Could not generate Excel file: %s', ME.message);
    end

    results_text = circ_text(params, comparison, imp1, cov_report, desc, ...
        excluded_channels, uncorrsigch, correctTFCEsigch, sigclusters, SnPMsigch, ...
        outputSname, nSubj, n1g, n2g);

    try
        generateAnalysisReport(results_struct, params, base_filename, outputSname, ...
            uncorrsigch, correctTFCEsigch, SnPMsigch);
        results_text{end+1} = ['Report saved to: ' outputSname '_report.html'];
    catch ME
        warning(ME.identifier, 'Could not generate HTML reports: %s', ME.message);
    end
    try
        export_results_csv(results_text, outputSname, params);
        results_text{end+1} = ['Result section CSV saved to: ' outputSname '_report.csv'];
    catch ME
        warning(ME.identifier, 'Could not write result CSV: %s', ME.message);
    end
end

% =========================================================================
% PERMUTATION CLOSURES
% Each asserts, on its FIRST call, that the analysed channel set is byte-
% identical between the observed map and the permuted map. A gate that is
% recomputed inside the loop decalibrates the max-statistic null silently;
% this makes that failure loud.
% =========================================================================
function fn = make_hotelling_perm(flctx, hctx, nCh, evaluable, obs)
    first = true;
    fn = @perm;
    function [s, pp] = perm()
        Yp = snpm_glm_permute(flctx);
        [s, pp] = hctx.eval(Yp(:, 1:nCh), Yp(:, nCh+1:end));
        if first, assert_gate_fixed(hctx.valid, evaluable, obs, s); first = false; end
    end
end

function fn = make_u2_perm(A, gi, sopt, evaluable, obs)
    first = true; n = numel(gi);
    fn = @perm;
    function [s, pp] = perm()
        [s, pp] = snpm_circ_watsons_u2(A, gi(randperm(n)), sopt);
        if first, assert_gate_fixed(sopt.valid_channels, evaluable, obs, s); first = false; end
    end
end

function fn = make_corr_perm(flctx, A, sopt, evaluable, obs)
    first = true;
    fn = @perm;
    function [s, pp] = perm()
        yp = snpm_glm_permute(flctx);
        [s, pp] = snpm_circ_corr_anglinear(yp, A, [], sopt);
        if first, assert_gate_fixed(sopt.valid_channels, evaluable, obs, s); first = false; end
    end
end

function assert_gate_fixed(used_mask, evaluable, obs, perm_stat)
    assert(isequal(logical(used_mask(:)'), logical(evaluable(:)')), ...
        'core_snpm:circGateDrift', ...
        ['The estimability mask handed to the statistic is not the one computed ' ...
         'from the observed labels. The gate must be fixed across permutations.']);
    a = isnan(obs); b = isnan(perm_stat);
    assert(isequal(a(evaluable), b(evaluable)), 'core_snpm:circGateDrift', ...
        ['The observed and permuted statistic maps are defined on different ' ...
         'channel sets (%d observed vs %d permuted finite channels inside the ' ...
         'gate). The max-statistic null would not be calibrated.'], ...
        sum(~a(evaluable)), sum(~b(evaluable)));
end

% =========================================================================
% GATES, DIAGNOSTICS, DESCRIPTIVES
% =========================================================================
function m = gate_mask(A, gi, is_grp, minn)
    if is_grp
        m = (sum(isfinite(A) & (gi == 1), 1) >= minn) & ...
            (sum(isfinite(A) & (gi == 0), 1) >= minn);
    else
        m = sum(isfinite(A), 1) >= minn;
    end
end

function report_u2_ties(A, evaluable)
% Non-gating tie diagnostic for the ECDF-based statistic: how many DISTINCT
% angles each analysed channel actually has, against its n.
    ch = find(evaluable);
    d = zeros(1, numel(ch)); nn = zeros(1, numel(ch));
    for k = 1:numel(ch)
        v = A(:, ch(k)); v = v(isfinite(v));
        nn(k) = numel(v); d(k) = numel(unique(v));
    end
    ratio = d ./ max(nn, 1);
    fprintf(['Watson U^2 tie diagnostic: distinct angles / n per channel -- ' ...
        'min %.2f, median %.2f, max %.2f (%d channels with ties)\n'], ...
        min(ratio), median(ratio), max(ratio), sum(d < nn));
end

function desc = circ_descriptive(A, gi, is_grp, evaluable)
% DESCRIPTIVE ONLY -- there is deliberately no whole-head inferential test on
% the circular tier (see the header). Rayleigh z/p are computed in closed form
% from the vendored circ_r; circ_rtest is not part of the vendored subset.
    if is_grp, groups = {gi == 1, gi == 0}; labels = {'Group A', 'Group B'};
    else,      groups = {true(size(A,1),1)}; labels = {'All subjects'};
    end
    ch = find(evaluable);
    G = struct('label', {}, 'n', {}, 'mean_deg', {}, 'sd_deg', {}, 'R', {}, ...
        'rayleigh_z', {}, 'rayleigh_p', {});
    for k = 1:numel(groups)
        a = A(groups{k}, ch); a = a(isfinite(a));
        n = numel(a);
        if n == 0
            G(k) = struct('label', labels{k}, 'n', 0, 'mean_deg', NaN, 'sd_deg', NaN, ...
                'R', NaN, 'rayleigh_z', NaN, 'rayleigh_p', NaN);
            continue
        end
        r  = abs(sum(exp(1i * a))) / n;
        mu = mod(rad2deg(angle(sum(exp(1i * a)))), 360);
        sd = rad2deg(sqrt(-2 * log(max(r, realmin))));      % circular SD
        z  = n * r^2;
        pv = exp(sqrt(1 + 4*n + 4*(n^2 - (n*r)^2)) - (1 + 2*n));
        G(k) = struct('label', labels{k}, 'n', n, 'mean_deg', mu, 'sd_deg', sd, ...
            'R', r, 'rayleigh_z', z, 'rayleigh_p', min(max(pv, 0), 1));
    end
    desc = struct('group', G, 'n_channels', numel(ch), 'note', ...
        ['Descriptive only. No whole-head inferential test is reported for the ' ...
         'circular tier: an arithmetic channel mean is invalid for angles, and a ' ...
         'circularly correct whole-head mean averages over a real 25-50 degree ' ...
         'anterior-posterior phase gradient. Inference is the channel-wise ' ...
         'permutation map above.']);
end

% =========================================================================
% FIGURES
% =========================================================================
function plot_circ_topos(A, gi, chanlocs, b1, b2, savepath, base_filename)
% Mean direction per group + the CIRCULAR difference. A linear colormap on a
% circular quantity renders 0 and 359 degrees as maximally different, so the
% colormap is cyclic (hsv) and the limits are fixed; the difference uses
% circ_dist, never the arithmetic mA - mB.
    mA = circ_mean_deg(A(gi == 1, :));
    mB = circ_mean_deg(A(gi == 0, :));
    d  = rad2deg(circ_dist(deg2rad(mA), deg2rad(mB)));
    try
        plot_topoInd(mA, mB, chanlocs, b1, b2, savepath, base_filename, ...
            'cmap1', 'hsv', 'clim1', [0 360], 'cmap2', 'hsv', 'clim2', [0 360], ...
            'diffdata', d, 'cmap_diff', 'hsv', 'clim_diff', [-180 180], ...
            'difftitle', 'Mean-direction difference (circ\_dist, deg)');
    catch ME
        warning(ME.identifier, 'Could not generate circular topology plots: %s', ME.message);
    end
end

function m = circ_mean_deg(A)
    z = exp(1i * A); z(isnan(A)) = 0;
    m = mod(rad2deg(angle(sum(z, 1))), 360);
end

% =========================================================================
% IO HELPERS
% =========================================================================
function S = read_wide(file, sheet)
% Wide table -> struct with .subj (normalised ids or []), .cols, .M.
    if isempty(sheet) || strcmpi(sheet, 'CSV File')
        Tb = readtable(file);
    else
        Tb = readtable(file, 'Sheet', sheet);
    end
    vn = Tb.Properties.VariableNames;
    issubj = contains(lower(vn), 'subject');
    if any(issubj)
        S.subj = normalise_ids(Tb{:, find(issubj, 1)});
    else
        S.subj = [];
    end
    S.cols = vn(~issubj);
    S.M    = Tb{:, ~issubj};
    if ~isnumeric(S.M)
        error('core_snpm:circNonNumeric', ...
            'Non-numeric data columns in %s.', file);
    end
end

function ids = normalise_ids(raw)
    if isnumeric(raw)
        ids = arrayfun(@(x) sprintf('sub%03d', x), raw, 'UniformOutput', false);
    else
        ids = cellstr(string(raw));
    end
    ids = ids(:);
end

function Cnt = load_counts(params, D1, D2, keptcols)
% Event counts are a REQUIRED input for the two phase-group analyses: one wide
% file per group, subjects x channels, with the SAME subject-identifier column
% and the SAME channel columns as the corresponding angle file. Channel names
% must be an exact SET match (order free, reordered here by name) and subject
% ids an exact set match after sub%03d normalisation; both are hard errors that
% name the offending entries.
    if ~isfield(params,'count1_file') || isempty(params.count1_file) || ...
       ~isfield(params,'count2_file') || isempty(params.count2_file)
        error('core_snpm:circCountsRequired', ...
            ['params.count1_file and params.count2_file are REQUIRED for ''%s''. ' ...
             'A subject''s preferred phase is estimated from their detected ' ...
             'events, so its precision -- and therefore the resultant length the ' ...
             'Hotelling test reads -- depends on how many events there were. ' ...
             'Without the count there is no way to tell a genuine phase ' ...
             'difference from a detection-rate difference. Supply one wide file ' ...
             'per group: subjects x channels, same subject column and same ' ...
             'channel columns as the corresponding angle file.'], params.comparison);
    end
    C1 = read_wide(params.count1_file, getf(params, 'count1_sheet', 'CSV File'));
    C2 = read_wide(params.count2_file, getf(params, 'count2_sheet', 'CSV File'));
    C1 = align_columns(C1, D1.cols, 'Group A count file');
    C2 = align_columns(C2, D1.cols, 'Group B count file');
    C1 = align_subjects(C1, D1.subj, 'Group A count file');
    C2 = align_subjects(C2, D2.subj, 'Group B count file');
    Cnt = [C1.M; C2.M];
    Cnt = Cnt(:, keptcols);
end

function S = align_columns(S, ref_cols, what)
% Exact SET match on channel names (order free); reorder into ref order.
    [tf, loc] = ismember(ref_cols, S.cols);
    missing = ref_cols(~tf);
    extra   = setdiff(S.cols, ref_cols);
    if ~isempty(missing) || ~isempty(extra)
        error('core_snpm:circColumnMismatch', ...
            ['%s must have exactly the same channel columns as the reference ' ...
             'angle file (order does not matter). Missing: {%s}. Unexpected: {%s}.'], ...
            what, strjoin(missing, ', '), strjoin(extra, ', '));
    end
    S.M    = S.M(:, loc);
    S.cols = ref_cols;
end

function S = align_subjects(S, ref_ids, what)
% Exact SET match on subject ids after sub%03d normalisation; reorder.
    if isempty(ref_ids)
        error('core_snpm:circSubjectColumn', ...
            ['The angle file has no Subject column, so %s cannot be matched to it. ' ...
             'Both files need a subject identifier column.'], what);
    end
    if isempty(S.subj)
        error('core_snpm:circSubjectColumn', '%s has no Subject column.', what);
    end
    [tf, loc] = ismember(ref_ids, S.subj);
    missing = ref_ids(~tf);
    extra   = setdiff(S.subj, ref_ids);
    if ~isempty(missing) || ~isempty(extra)
        error('core_snpm:circSubjectMismatch', ...
            ['%s must contain exactly the same subjects as its angle file. ' ...
             'Missing: {%s}. Unexpected: {%s}.'], what, ...
            strjoin(missing(:)', ', '), strjoin(cellstr(string(extra(:)))', ', '));
    end
    S.M    = S.M(loc, :);
    S.subj = ref_ids;
end

function [A, y] = match_measure(A1, D1, D2, params)
% circ_corrAngLinear: angles (wide) + a single linear measure column, paired by
% subject.
    if isempty(D1.subj) || isempty(D2.subj)
        error('core_snpm:circSubjectColumn', ...
            'circ_corrAngLinear needs a Subject column in both the angle and the measure file.');
    end
    col = getf(params, 'measure_col', '');
    if isempty(col)
        if numel(D2.cols) ~= 1
            error('core_snpm:circMeasureColumn', ...
                ['The measure file has %d non-subject columns {%s}; name the one to ' ...
                 'use in params.measure_col.'], numel(D2.cols), strjoin(D2.cols, ', '));
        end
        col = D2.cols{1};
    end
    k = find(strcmp(D2.cols, col), 1);
    if isempty(k)
        error('core_snpm:circMeasureColumn', ...
            'params.measure_col ''%s'' is not a column of the measure file.', col);
    end
    [subj, i1, i2] = intersect(D1.subj, D2.subj, 'stable');
    if numel(subj) < 8
        error('core_snpm:circTooFewSubjects', ...
            'Only %d subjects match between the angle and measure files (need >= 8).', ...
            numel(subj));
    end
    A = A1(i1, :);
    y = D2.M(i2, k);
end

% =========================================================================
% SMALL UTILITIES
% =========================================================================
function v = getf(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function s = stat_label(comparison)
    switch comparison
        case 'circ_phase_group',    s = 'Hotelling T^2 on (cos,sin), two groups';
        case 'circ_phase_group_u2', s = 'Watson U^2, two groups';
        otherwise,                  s = 'Circular-linear association (F)';
    end
end

function excl = excluded_info(evaluable, chanlocs, context_label)
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
        'reason', ['not estimable: fewer than 8 subjects per group, incomplete ' ...
                   'column under a per-channel covariate, or log(count) collinear with group']);
end

function txt = circ_text(params, comparison, imp, cov, desc, excl, unc, tfce, ...
        sigclusters, snpmch, outputSname, nSubj, n1g, n2g)
    txt = {
        '=== SnPM Circular (Tier 2) Analysis Results ===', '', ...
        ['Comparison: ' comparison ' (' stat_label(comparison) ')'], ...
        ['Angle units declared: ' imp.units], ...
        ['Zero convention: ' imp.convention], ...
        sprintf('Rotation applied: %+.1f deg (0 = slow-oscillation up-state peak)', imp.rotation_deg), ...
        sprintf('Angular resolution of the export: %.4g deg', imp.resolution_deg), ...
        ['Event-count covariate: ' cov.mode], ...
        ['Permutations: ' num2str(params.permutations)], ...
        sprintf('Subjects: %d (group A %d, group B %d)', nSubj, n1g, n2g) }';
    if ~strcmp(cov.mode, 'none')
        txt{end+1} = sprintf(['Count/group confound: point-biserial r median %.2f, ' ...
            'max |r| %.2f; covariate dropped at %d channels (constant), ' ...
            '%d channels non-estimable (collinear)'], cov.pb_median, cov.pb_max, ...
            cov.n_dropped_constant, cov.n_collinear);
    end
    txt{end+1} = sprintf('Excluded channels: %d of %d', excl.n, excl.n_channels);
    txt{end+1} = '';
    txt{end+1} = '--- Descriptive circular summary (NO inference) ---';
    for k = 1:numel(desc.group)
        g = desc.group(k);
        txt{end+1} = sprintf(['%s: n=%d angles, mean direction %.1f deg, circular ' ...
            'SD %.1f deg, R=%.3f, Rayleigh z=%.2f p=%.3g'], g.label, g.n, ...
            g.mean_deg, g.sd_deg, g.R, g.rayleigh_z, g.rayleigh_p); %#ok<AGROW>
    end
    txt{end+1} = desc.note;
    txt{end+1} = '';
    txt{end+1} = '--- Statistical Results ---';
    txt{end+1} = ['Uncorrected significant channels: ' num2str(numel(unc))];
    txt{end+1} = ['TFCE corrected significant channels: ' num2str(numel(tfce))];
    txt{end+1} = ['Significant clusters N = : ' num2str(numel(sigclusters))];
    txt{end+1} = ['SnPM significant channels N = : ' num2str(numel(snpmch))];
    txt{end+1} = ['Cluster statistic: MASS (sum of the statistic / n channels). ' ...
        'This differs from the FieldTrip maxsize used by Helfrich et al.'];
    txt{end+1} = ['Note: the parametric per-channel p is a cluster-forming and ' ...
        'descriptive quantity only. For Watson U^2 at n=8 per group the exact ' ...
        'enumerated null gives 0.1553 / 0.1787 / 0.2412 at alpha .10 / .05 / .01 ' ...
        'against the asymptotic 0.152 / 0.187 / 0.268 -- good agreement at .10 ' ...
        'and .05, conservative in the far tail. The inferential p is always the ' ...
        'permutation p.'];
    txt{end+1} = '';
    txt{end+1} = ['Results saved to: ' outputSname '.mat'];
    txt{end+1} = ['Excel file saved to: ' outputSname '.xlsx'];
end
