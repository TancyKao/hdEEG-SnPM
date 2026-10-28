function D = snpm_glm_design(preset, meta, opts)
% Translate a friendly preset + column mapping into the GLM objects the
% engine needs: design matrix X, contrast C (1 row -> t, >1 -> F), the
% nuisance columns for Freedman-Lane permutation, exchangeability blocks,
% and the permutation type. This is where all the "no statistics knowledge
% required" logic lives.
%
% INPUTS
%   preset : 'anova1' | 'ancova' | 'regression' | 'rmanova' | 'mixed2way'
%   meta   : table, one row per observation (subject for between-subject
%            presets; subject x condition for rmanova/mixed2way)
%   opts   : struct of column names -
%       .group_col, .condition_col, .subject_col, .predictor_col,
%       .covariate_cols (cellstr), .effect ('interaction'|'group'|'condition'
%       for mixed2way; default 'interaction'),
%       .condition_levels (optional; rmanova/mixed2way) the condition level
%       order to use instead of first appearance in these rows
%
% OUTPUT D with fields:
%   X            : nObs x p full design
%   C            : q x p contrast (q==1 -> t, q>1 -> F)
%   contrast_type: 't' | 'F'
%   nuisance_idx : logical 1 x p, the columns held fixed under Freedman-Lane
%   eb           : nObs x 1 exchangeability-block id (subject for 'within')
%   perm_type    : 'free' (between) | 'within' (repeated measures)
%   effect_label : human-readable effect name
%   posthoc      : struct array (label, a, b) of pairwise comparisons that
%                  follow an effect with three or more levels: mean(a) - mean(b),
%                  a and b being level indices in order of first appearance in
%                  the factor column. Each is run as its OWN two-level test on
%                  the rows of a and b only -- see snpm_glm_pair. (Until
%                  2026-10-14 this carried a contrast row C in the omnibus
%                  model; that scheme left the other levels' effects in the
%                  permuted residuals and is retired.)
%   factor_col   : name of the column whose levels .posthoc compares ('' if none)
%   col_names    : 1 x p cellstr, one name per column of X, for messages:
%                  '(Intercept)'; 'predictor <col>'; 'covariate <col>' (numeric
%                  or logical) or 'covariate <col>_<level>' (dummy, reference
%                  level dropped); '<group_col>_<level>', '<condition_col>_<level>',
%                  '<subject_col>_<level>'; '<group_col>_<l>:<condition_col>_<l>'.
%                  The predictor/covariate prefixes tell the two apart when the
%                  same column is on both sides (snpm_glm_estimability).

    D = struct('X',[], 'C',[], 'contrast_type','t', 'nuisance_idx',[], ...
        'eb',[], 'perm_type','free', 'effect_label','', ...
        'posthoc',struct('label',{},'a',{},'b',{}), 'factor_col','', 'col_names',{{}});

    n = height(meta);
    [cov, cov_names] = get_covariates(meta, opts);   % nObs x nCov, mean-centered (may be empty)

    switch lower(preset)
        case 'anova1'
            [Dg, gl] = ref_dummy(meta.(opts.group_col));
            X = [ones(n,1), Dg];
            D.col_names = ['(Intercept)', level_names(opts.group_col, gl)];
            ofint = false(1, size(X,2)); ofint(2:end) = true;
            D.X = X;
            D.C = contrast_from_idx(ofint);
            D.nuisance_idx = ~ofint;                 % intercept only
            D.perm_type = 'free';  D.eb = ones(n,1);
            D.effect_label = sprintf('group (%d levels)', numel(gl));
            D.posthoc = pairwise_levels(gl);  D.factor_col = opts.group_col;

        case 'ancova'
            [Dg, gl] = ref_dummy(meta.(opts.group_col));
            X = [ones(n,1), Dg, cov];
            D.col_names = ['(Intercept)', level_names(opts.group_col, gl), cov_names];
            ng = size(Dg,2);
            ofint = false(1, size(X,2)); ofint(1+(1:ng)) = true;
            D.X = X;
            D.C = contrast_from_idx(ofint);
            D.nuisance_idx = ~ofint;                 % intercept + covariates
            D.perm_type = 'free';  D.eb = ones(n,1);
            D.effect_label = sprintf('group (%d levels), covariate-adjusted', numel(gl));
            D.posthoc = pairwise_levels(gl);  D.factor_col = opts.group_col;

        case 'regression'
            pred = meta.(opts.predictor_col);
            pred = pred - mean(pred, 'omitnan');
            X = [ones(n,1), pred, cov];
            D.col_names = ['(Intercept)', {['predictor ' opts.predictor_col]}, cov_names];
            ofint = false(1, size(X,2)); ofint(2) = true;
            D.X = X;
            D.C = contrast_from_idx(ofint);
            D.contrast_type = 't';
            D.nuisance_idx = ~ofint;                 % intercept + covariates
            D.perm_type = 'free';  D.eb = ones(n,1);
            D.effect_label = sprintf('slope of %s', opts.predictor_col);

        case 'rmanova'
            [Ds, sl] = ref_dummy(meta.(opts.subject_col));
            [Dc, cl] = ref_dummy(meta.(opts.condition_col), condition_order(opts));
            X = [ones(n,1), Ds, Dc];
            D.col_names = ['(Intercept)', level_names(opts.subject_col, sl), ...
                level_names(opts.condition_col, cl)];
            nc = size(Dc,2);
            ofint = false(1, size(X,2)); ofint(end-nc+1:end) = true;
            D.X = X;
            D.C = contrast_from_idx(ofint);
            D.nuisance_idx = ~ofint;                 % intercept + subject means
            D.perm_type = 'within';
            D.eb = block_ids(meta.(opts.subject_col));
            D.effect_label = sprintf('condition (%d levels), within-subject', numel(cl));
            % Pairwise condition post-hoc: each pair is its own two-level
            % rmanova on the rows of those two conditions (snpm_glm_pair).
            D.posthoc = pairwise_levels(cl);  D.factor_col = opts.condition_col;

        case 'mixed2way'
            effect = 'interaction';
            if isfield(opts,'effect') && ~isempty(opts.effect), effect = lower(opts.effect); end
            % Subject intercepts absorb the subject and (nested) group main
            % effect, so the group main effect is not estimable here — it is a
            % between-subjects test on subject means (core_snpm_glm handles it).
            [Ds, sl] = ref_dummy(meta.(opts.subject_col));
            [Dc, cl] = ref_dummy(meta.(opts.condition_col), condition_order(opts));
            ns = size(Ds,2); nc = size(Dc,2);
            base_names = ['(Intercept)', level_names(opts.subject_col, sl), ...
                level_names(opts.condition_col, cl)];
            switch effect
                case 'interaction'
                    % Full model; test the group x condition product columns.
                    [Dg, gl] = ref_dummy(meta.(opts.group_col));
                    Dx = interaction_cols(Dg, Dc); nx = size(Dx,2);
                    X = [ones(n,1), Ds, Dc, Dx];
                    D.col_names = [base_names, interaction_names(opts, gl, cl)];
                    ofint = false(1, size(X,2)); ofint(end-nx+1:end) = true;
                    D.effect_label = 'two-way mixed ANOVA (group x condition)';
                case 'condition'
                    % FULL model (Method Spec C, approved 2026-10-14):
                    %     X = [1, Ds, Dc, Dxc]
                    % with the group-by-condition columns Dxc built from a
                    % group code that makes the Dc coefficients the condition
                    % MAIN effect, so a real interaction sits in the nuisance
                    % model instead of in the residuals that are permuted.
                    %
                    % The model this replaces was the reduced [1, Ds, Dc]. It
                    % gave the pooled condition effect, but left any real
                    % interaction in the Freedman-Lane residuals: with a
                    % crossover interaction a planted condition effect was
                    % detected 0.000 of the time. Plain reference coding in
                    % the full model would not do either: Dc would then be the
                    % condition effect in the reference group only.
                    %
                    % The condition F from this model is the standard mixed-
                    % ANOVA F (ranova, and Type III anovan): 8 + 14 subjects x
                    % 3 conditions gives F(2,40) = 30.9004 in all three.
                    [Dg, gl] = ref_dummy(meta.(opts.group_col));
                    Eg  = main_effect_group_code(Dg);
                    Dxc = interaction_cols(Eg, Dc);      % group outer, condition inner
                    X = [ones(n,1), Ds, Dc, Dxc];
                    D.col_names = [base_names, interaction_names(opts, gl, cl)];
                    ofint = false(1, size(X,2)); ofint(1+ns+(1:nc)) = true;
                    D.effect_label = sprintf('condition main effect (%d levels)', numel(cl));
                    % Pairwise condition post-hoc: each pair is this same
                    % model on the rows of those two conditions (snpm_glm_pair).
                    D.posthoc = pairwise_levels(cl);  D.factor_col = opts.condition_col;
                otherwise
                    error('snpm_glm_design:mixed2wayGroup', ...
                        ['For mixed2way the between-group main effect is a subject-level test; ' ...
                         'use anova1/ancova on subject means for that. Supported effects here: ' ...
                         '''interaction'' (default) or ''condition''.']);
            end
            D.X = X;
            D.C = contrast_from_idx(ofint);
            D.nuisance_idx = ~ofint;                 % everything except the effect
            D.perm_type = 'within';
            D.eb = block_ids(meta.(opts.subject_col));

        otherwise
            error('snpm_glm_design:unknownPreset', 'Unknown preset: %s', preset);
    end

    % contrast type from number of rows
    if size(D.C, 1) == 1
        D.contrast_type = 't';
    else
        D.contrast_type = 'F';
    end
end

% ---------------------------------------------------------------------------
function Eg = main_effect_group_code(Dg)
% =========================================================================
% GROUP CODING FOR THE mixed2way CONDITION MAIN EFFECT -- THE ONE PLACE.
%
% SUM-TO-ZERO: the reference-coded group dummies with the first group's rows
% set to -1 in every column. The condition coefficients are then the condition
% effect averaged over groups with EQUAL weight per group (unweighted marginal
% means), which is what a standard mixed ANOVA (ranova, Type III) reports.
%
% The alternative is CENTRED coding (Dg - mean(Dg)), which weights each group
% by its number of subjects. The two agree when groups are the same size and
% differ otherwise. Which one the toolbox should report is still to be
% confirmed by Tancy (open as of 2026-10-14); to switch, change this function
% only -- nothing else depends on the choice.
% =========================================================================
    first_group = (sum(Dg, 2) == 0);
    Eg = Dg - first_group;
end

% ---------------------------------------------------------------------------
function order = condition_order(opts)
% opts.condition_levels (optional): the condition level order fixed from
% outside, so that the condition effect inside one group lists and signs its
% pairs exactly as the whole table does. See snpm_glm_levels.
    order = [];
    if isfield(opts, 'condition_levels') && ~isempty(opts.condition_levels)
        order = opts.condition_levels;
    end
end

function [Dd, levels] = ref_dummy(x, order)
% Reference-coded dummies (first level = reference): n x (k-1). Levels are in
% order of first appearance unless `order` fixes them (snpm_glm_levels).
    if nargin < 2, order = []; end
    [levels, gi] = snpm_glm_levels(x, order);
    k = numel(levels);
    n = numel(gi);
    Dd = zeros(n, max(k-1, 0));
    for j = 2:k
        Dd(:, j-1) = (gi == j);
    end
end

function b = block_ids(x)
    [~, ~, b] = unique(x, 'stable');
end

function C = contrast_from_idx(ofint)
% Build a q x p contrast selecting the of-interest columns (identity rows).
    cols = find(ofint);
    q = numel(cols);
    p = numel(ofint);
    C = zeros(q, p);
    for i = 1:q
        C(i, cols(i)) = 1;
    end
end

function [cov, names] = get_covariates(meta, opts)
% Numeric covariates pass through; categorical/string/cellstr nuisance columns
% are dummy-coded (k-1 columns, first level = reference) so e.g. sex can be a
% covariate alongside age. All columns are mean-centred. names: one per
% column of cov ('covariate <col>' or 'covariate <col>_<level>').
    cov = [];
    names = {};
    if ~isfield(opts, 'covariate_cols'), return; end
    cols = snpm_glm_covariate_list(opts.covariate_cols);   % blank names dropped
    if isempty(cols), return; end
    parts = {};
    for i = 1:numel(cols)
        c = char(cols{i});
        v = meta.(c);
        if isnumeric(v) || islogical(v)
            parts{end+1} = double(v); %#ok<AGROW>
            names{end+1} = ['covariate ' c]; %#ok<AGROW>
        else
            % removecats: code from the levels PRESENT in these rows. Without
            % it a categorical column subset to two groups keeps a category
            % nobody has, and dummyvar returns an all-zero column for it.
            vc = removecats(categorical(v));
            D = double(dummyvar(vc));            % n x k
            cats = categories(vc);
            if size(D, 2) >= 2
                D = D(:, 2:end);                 % drop reference level
                names = [names, strcat({['covariate ' c '_']}, cats(2:end)')]; %#ok<AGROW>
            else
                D = [];
            end
            parts{end+1} = D; %#ok<AGROW>
        end
    end
    cov = [parts{:}];
    if ~isempty(cov)
        % A covariate that does not vary over these rows is the intercept
        % again and carries nothing; drop it (it would make X rank-deficient).
        % Tested on the raw column, before centring, so it is exact.
        varies = max(cov, [], 1, 'omitnan') > min(cov, [], 1, 'omitnan');
        cov = cov(:, varies);
        names = names(varies);
        cov = cov - mean(cov, 1, 'omitnan');
    end
end

function names = level_names(col, levels)
% '<col>_<level>' for every non-reference level (levels 2..k), as ref_dummy codes them
    names = cell(1, max(numel(levels) - 1, 0));
    for j = 2:numel(levels)
        names{j-1} = sprintf('%s_%s', col, tostr(levels(j)));
    end
end

function names = interaction_names(opts, gl, cl)
% names of interaction_cols: group outer, condition inner
    g = level_names(opts.group_col, gl);
    c = level_names(opts.condition_col, cl);
    names = cell(1, numel(g) * numel(c));
    k = 0;
    for a = 1:numel(g)
        for b = 1:numel(c)
            k = k + 1;
            names{k} = [g{a} ':' c{b}];
        end
    end
end

function Dx = interaction_cols(Dg, Dc)
% Pairwise products of group and condition dummies -> interaction columns.
    ng = size(Dg, 2); nc = size(Dc, 2);
    Dx = zeros(size(Dg, 1), ng * nc);
    col = 0;
    for a = 1:ng
        for b = 1:nc
            col = col + 1;
            Dx(:, col) = Dg(:, a) .* Dc(:, b);
        end
    end
end

function ph = pairwise_levels(levels)
% Every pair of levels, in level order: mean(a) - mean(b), a before b. The
% levels are in order of first appearance (unique(...,'stable')), and a / b
% are indices into that order.
    k = numel(levels);
    ph = struct('label', {}, 'a', {}, 'b', {});
    idx = 0;
    for a = 1:k
        for b = a+1:k
            idx = idx + 1;
            ph(idx).label = sprintf('%s vs %s', tostr(levels(a)), tostr(levels(b)));
            ph(idx).a = a;
            ph(idx).b = b;
        end
    end
end

function s = tostr(v)
    if iscell(v), v = v{1}; end
    if isnumeric(v), s = num2str(v); else, s = char(string(v)); end
end
