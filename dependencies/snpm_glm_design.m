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
%       for mixed2way; default 'interaction')
%
% OUTPUT D with fields:
%   X            : nObs x p full design
%   C            : q x p contrast (q==1 -> t, q>1 -> F)
%   contrast_type: 't' | 'F'
%   nuisance_idx : logical 1 x p, the columns held fixed under Freedman-Lane
%   eb           : nObs x 1 exchangeability-block id (subject for 'within')
%   perm_type    : 'free' (between) | 'within' (repeated measures)
%   effect_label : human-readable effect name
%   posthoc      : struct array (label, C) of pairwise t-contrasts (omnibus only)

    D = struct('X',[], 'C',[], 'contrast_type','t', 'nuisance_idx',[], ...
        'eb',[], 'perm_type','free', 'effect_label','', 'posthoc',struct('label',{},'C',{}));

    n = height(meta);
    cov = get_covariates(meta, opts);   % nObs x nCov, mean-centered (may be empty)

    switch lower(preset)
        case 'anova1'
            [Dg, gl] = ref_dummy(meta.(opts.group_col));
            X = [ones(n,1), Dg];
            ofint = false(1, size(X,2)); ofint(2:end) = true;
            D.X = X;
            D.C = contrast_from_idx(ofint);
            D.nuisance_idx = ~ofint;                 % intercept only
            D.perm_type = 'free';  D.eb = ones(n,1);
            D.effect_label = sprintf('group (%d levels)', numel(gl));
            D.posthoc = pairwise_contrasts(gl, 1);   % group dummies start at col 2

        case 'ancova'
            [Dg, gl] = ref_dummy(meta.(opts.group_col));
            X = [ones(n,1), Dg, cov];
            ng = size(Dg,2);
            ofint = false(1, size(X,2)); ofint(1+(1:ng)) = true;
            D.X = X;
            D.C = contrast_from_idx(ofint);
            D.nuisance_idx = ~ofint;                 % intercept + covariates
            D.perm_type = 'free';  D.eb = ones(n,1);
            D.effect_label = sprintf('group (%d levels), covariate-adjusted', numel(gl));
            D.posthoc = pairwise_contrasts(gl, 1);

        case 'regression'
            pred = meta.(opts.predictor_col);
            pred = pred - mean(pred, 'omitnan');
            X = [ones(n,1), pred, cov];
            ofint = false(1, size(X,2)); ofint(2) = true;
            D.X = X;
            D.C = contrast_from_idx(ofint);
            D.contrast_type = 't';
            D.nuisance_idx = ~ofint;                 % intercept + covariates
            D.perm_type = 'free';  D.eb = ones(n,1);
            D.effect_label = sprintf('slope of %s', opts.predictor_col);

        case 'rmanova'
            [Ds, ~]  = ref_dummy(meta.(opts.subject_col));
            [Dc, cl] = ref_dummy(meta.(opts.condition_col));
            X = [ones(n,1), Ds, Dc];
            nc = size(Dc,2);
            ofint = false(1, size(X,2)); ofint(end-nc+1:end) = true;
            D.X = X;
            D.C = contrast_from_idx(ofint);
            D.nuisance_idx = ~ofint;                 % intercept + subject means
            D.perm_type = 'within';
            D.eb = block_ids(meta.(opts.subject_col));
            D.effect_label = sprintf('condition (%d levels), within-subject', numel(cl));
            % Pairwise condition post-hoc: condition dummies sit after the
            % intercept + subject-mean columns, so they start at 1 + size(Ds,2).
            % Same reference coding as anova1, so mean_a - mean_b holds; the
            % within-subject permutation (perm_type/eb above) is inherited.
            D.posthoc = pairwise_contrasts(cl, 1 + size(Ds,2));

        case 'mixed2way'
            effect = 'interaction';
            if isfield(opts,'effect') && ~isempty(opts.effect), effect = lower(opts.effect); end
            % Subject intercepts absorb the subject and (nested) group main
            % effect, so the group main effect is not estimable here — it is a
            % between-subjects test on subject means (core_snpm_glm handles it).
            [Ds, ~]  = ref_dummy(meta.(opts.subject_col));
            [Dc, cl] = ref_dummy(meta.(opts.condition_col));
            ns = size(Ds,2); nc = size(Dc,2);
            switch effect
                case 'interaction'
                    % Full model; test the group x condition product columns.
                    [Dg, ~] = ref_dummy(meta.(opts.group_col));
                    Dx = interaction_cols(Dg, Dc); nx = size(Dx,2);
                    X = [ones(n,1), Ds, Dc, Dx];
                    ofint = false(1, size(X,2)); ofint(end-nx+1:end) = true;
                    D.effect_label = 'two-way mixed ANOVA (group x condition)';
                case 'condition'
                    % Reduced model WITHOUT the interaction columns, so the Dc
                    % coefficients are the condition effect POOLED over groups
                    % (the marginal main effect). With the interaction in the
                    % model and reference coding, Dc would instead be the
                    % condition effect in the reference group only.
                    X = [ones(n,1), Ds, Dc];
                    ofint = false(1, size(X,2)); ofint(1+ns+(1:nc)) = true;
                    D.effect_label = sprintf('condition main effect (%d levels)', numel(cl));
                    % Pairwise condition post-hoc: condition dummies start at
                    % 1+ns; contrast width already equals size(X,2).
                    D.posthoc = pairwise_contrasts(cl, 1 + ns);
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
function [Dd, levels] = ref_dummy(x)
% Reference-coded dummies (first level = reference): n x (k-1).
    [levels, ~, gi] = unique(x, 'stable');
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

function cov = get_covariates(meta, opts)
% Numeric covariates pass through; categorical/string/cellstr nuisance columns
% are dummy-coded (k-1 columns, first level = reference) so e.g. sex can be a
% covariate alongside age. All columns are mean-centred.
    cov = [];
    if ~isfield(opts, 'covariate_cols') || isempty(opts.covariate_cols), return; end
    cols = opts.covariate_cols;
    if ischar(cols) || isstring(cols), cols = cellstr(cols); end
    parts = {};
    for i = 1:numel(cols)
        v = meta.(char(cols{i}));
        if isnumeric(v) || islogical(v)
            parts{end+1} = double(v); %#ok<AGROW>
        else
            vc = categorical(v);
            D = double(dummyvar(vc));            % n x k
            if size(D, 2) >= 2, D = D(:, 2:end); else, D = []; end   % drop reference level
            parts{end+1} = D; %#ok<AGROW>
        end
    end
    cov = [parts{:}];
    if ~isempty(cov), cov = cov - mean(cov, 1, 'omitnan'); end
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

function ph = pairwise_contrasts(levels, dummy_start_col)
% t-contrasts for every group pair, in the reference-coded design where
% column 1 is the intercept and group dummies (levels 2..k) start at
% (dummy_start_col+1). mean_1 = intercept; mean_j = intercept + beta_j.
    k = numel(levels);
    p = dummy_start_col + (k - 1);    % intercept(s) + (k-1) dummies
    ph = struct('label', {}, 'C', {});
    idx = 0;
    for a = 1:k
        for b = a+1:k
            idx = idx + 1;
            c = zeros(1, p);
            % difference mean_a - mean_b in terms of dummy betas
            if a > 1, c(dummy_start_col + (a-1)) = c(dummy_start_col + (a-1)) + 1; end
            if b > 1, c(dummy_start_col + (b-1)) = c(dummy_start_col + (b-1)) - 1; end
            ph(idx).label = sprintf('%s vs %s', tostr(levels(a)), tostr(levels(b)));
            ph(idx).C = c;
        end
    end
end

function s = tostr(v)
    if iscell(v), v = v{1}; end
    if isnumeric(v), s = num2str(v); else, s = char(string(v)); end
end
