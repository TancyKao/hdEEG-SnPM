function [stat, wald, p_model] = snpm_lmm_fit(power, meta, spec, evaluable)
% Fit a linear mixed model at every channel and return the statistic of
% the effect of interest. This is the per-channel engine used by
% snpm_lmm_TFCE and snpm_lmm_cluster (analogous to the per-channel ttest
% inside the existing SnPM functions, but using fitlme).
%
% INPUTS
%   power : trials x nCh matrix of channel values (the per-channel predictor,
%           e.g. absolute band power). One row per trial/awakening.
%   meta  : table with one row per trial holding the model variables, e.g.
%           the dependent variable, Subject, group, time. The per-channel
%           predictor is injected as a variable named 'POWER' before fitting,
%           so spec.fixed must reference 'POWER'.
%   evaluable : (optional) 1 x nCh logical. Channels that are FALSE are not
%           fitted at all and come back stat=NaN, wald=0, p_model=NaN. The
%           mask is decided ONCE by the caller (snpm_lmm_TFCE /
%           snpm_lmm_cluster / core_snpm_lmm) from the channel's completeness
%           and is held fixed across the observed fit and every permuted fit,
%           so the observed and null maps are defined on one channel set.
%           Omit (or pass []) to fit every channel, which is the historical
%           behaviour and is identical on complete data.
%   spec  : struct describing the model:
%       .dv          (char) name of the dependent-variable column in meta
%       .fixed       (char) fixed-effects part of the formula, referencing
%                    'POWER' and any meta columns, e.g. 'POWER + group'
%       .random      (char) random-effects part, e.g. '(1|Subject) + (1|time)'
%       .effect      (char) the term whose effect is mapped, e.g. 'POWER'
%                    or 'group'
%       .effect_type (char) 'continuous' -> signed t statistic (1 df)
%                           'factor'     -> F statistic from anova (k-1 df)
%       .categorical (cellstr, optional) meta columns to force categorical.
%                    Default {'Subject','group'} when present.
%
% OUTPUTS (each 1 x nCh)
%   stat    : signed t (continuous) or F (factor) for the effect of interest.
%             NaN for channels that could not be fit, and for channels whose
%             statistic came back non-finite (fitlme does NOT throw when the
%             effect is inestimable on the rows it was given -- e.g. a group
%             factor that is constant among a channel's non-missing rows --
%             it returns FStat = NaN from anova(lme)).
%   wald    : t^2 (continuous) or F (factor). 0 whenever stat is not finite,
%             so an unfit or undefined channel cannot contribute to (or
%             poison, via NaN) cluster mass. This is the documented contract
%             and it is now enforced on BOTH the continuous and factor paths.
%   p_model : model p-value of the effect of interest. NaN for unfit channels.

    nCh = size(power, 2);

    if nargin < 4 || isempty(evaluable)
        evaluable = true(1, nCh);
    end
    evaluable = reshape(logical(evaluable), 1, []);
    if numel(evaluable) ~= nCh
        error('snpm_lmm_fit:evaluableSize', ...
            'evaluable mask has %d entries but power has %d channels.', ...
            numel(evaluable), nCh);
    end

    dv          = spec.dv;
    fixed       = spec.fixed;
    random      = spec.random;
    effect      = spec.effect;
    effect_type = lower(spec.effect_type);

    if isfield(spec, 'categorical') && ~isempty(spec.categorical)
        cat_vars = spec.categorical;
    else
        cat_vars = intersect({'Subject', 'group'}, meta.Properties.VariableNames);
    end

    % Force the requested columns to categorical once, up front.
    base_meta = meta;
    for v = 1:numel(cat_vars)
        cv = cat_vars{v};
        if ismember(cv, base_meta.Properties.VariableNames) && ~iscategorical(base_meta.(cv))
            base_meta.(cv) = categorical(base_meta.(cv));
        end
    end

    formula = sprintf('%s ~ %s + %s', dv, fixed, random);

    stat    = nan(1, nCh);
    wald    = zeros(1, nCh);
    p_model = nan(1, nCh);

    % minimum trials needed before a fit is attempted
    min_trials = 5;

    parfor ch = 1:nCh
        if ~evaluable(ch)
            continue   % masked out by the caller: stat NaN, wald 0, p NaN
        end
        col = power(:, ch);
        valid = ~isnan(col);
        if sum(valid) < min_trials
            continue
        end

        tbl = base_meta(valid, :);
        tbl.POWER = col(valid);

        try
            lme = fitlme(tbl, formula);
        catch
            continue   % singular / non-convergent fit -> leave as NaN/0
        end

        switch effect_type
            case 'continuous'
                % single-coefficient (1 df) signed t statistic
                names = lme.Coefficients.Name;
                row = find(strcmp(names, effect), 1);
                if isempty(row)
                    % fall back to first coefficient that contains the name
                    row = find(contains(names, effect), 1);
                end
                % A non-finite t means the coefficient is not estimable on the
                % rows supplied; leave the channel at stat=NaN / wald=0 rather
                % than propagating NaN into the cluster mass or the max null.
                if ~isempty(row) && isfinite(lme.Coefficients.tStat(row))
                    t = lme.Coefficients.tStat(row);
                    stat(ch)    = t;
                    wald(ch)    = t^2;
                    p_model(ch) = lme.Coefficients.pValue(row);
                end

            case 'factor'
                % omnibus (k-1 df) F statistic for the categorical term
                av = anova(lme);
                row = find(strcmp(av.Term, effect), 1);
                % anova(lme) RETURNS the term with FStat = NaN when the factor
                % is constant among the fitted rows (no error is thrown), so
                % the finiteness check is what keeps the wald contract honest.
                if ~isempty(row) && isfinite(av.FStat(row))
                    f = av.FStat(row);
                    stat(ch)    = f;
                    wald(ch)    = f;     % already a ratio; summed for cluster mass
                    p_model(ch) = av.pValue(row);
                end

            otherwise
                error('snpm_lmm_fit:badEffectType', ...
                    'spec.effect_type must be ''continuous'' or ''factor''');
        end
    end
end
