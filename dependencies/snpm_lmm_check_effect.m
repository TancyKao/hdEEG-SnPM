function effect = snpm_lmm_check_effect(meta, spec)
% Refuse a mixed-model effect that the permutation cannot test or the model
% does not contain, before any channel is fitted.
%
%   effect = snpm_lmm_check_effect([],   spec)   % before the file is read
%   effect = snpm_lmm_check_effect(meta, spec)   % after the missing-row drop
%
% spec fields: dv, fixed, effect, effect_type ('continuous'|'factor'), perm,
% group, and (with meta) categorical (the columns snpm_lmm_fit converts).
%
% THE NAMES. With meta, the coefficient and term names are those of the fixed
% part fitted to the metadata exactly as snpm_lmm_fit prepares it
% (snpm_lmm_categorize: same columns converted to categorical, rest as read)
% with the channel appended last as POWER. They come from ONE ordinary
% least-squares fit (fitlm) of the fixed part with a fixed synthetic POWER
% column (no random draws, no channel data): fitlm and fitlme build the fixed
% design through the same code, so the names, their order and the reference
% level of every factor are fitlme's -- including the first-appearance level
% order of a text column that is not converted (see snpm_lmm_categorize).
% fitlm is used rather than fitlme because it also names a rank-deficient
% design, which fitlme refuses.
%
% 1. group_label: the effect must involve the group, or relabelling subjects
%    leaves every permuted map equal to the observed one
%    (core_snpm_lmm:effectNotPermuted). A token of the effect (split on ':'
%    and '*') names the group when it IS the group column, or is a coefficient
%    of it ('<group>_<level>' as fitlme names it: group_B, group_true, inside
%    an interaction group_B:POWER). Before the file is read only the coarse
%    form (token starts with '<group>_') can be checked. A token '<group>_x'
%    that is neither a group coefficient nor another column (group_C for a
%    level not in the data) is left to check 2, whose message lists the real
%    coefficients; another column such as group_size is never the group.
% 2. meta given: the effect must name a coefficient (continuous) or a term
%    (factor) of the fixed part, else every channel's statistic would be NaN
%    (core_snpm_lmm:unknownEffect). Continuous matching is the one snpm_lmm_fit
%    applies (exact name, else first name containing the effect). Factor
%    matching ignores the order of the variables inside an interaction
%    ('POWER:group' finds 'group:POWER'); the returned effect is then the
%    term exactly as anova(lme) names it, so the fit's exact match finds it.
%
% Message texts are shared with the Python port (lmm.py); do not reword.

    effect = char(spec.effect);
    grp    = char(spec.group);
    is_grp_perm = strcmpi(spec.perm, 'group_label');
    toks = strtrim(regexp(effect, '[:*]', 'split'));

    if isempty(meta)
        if is_grp_perm && ~any(strcmp(toks, grp) | startsWith(toks, [grp '_']))
            not_permuted(grp, effect);
        end
        return
    end

    [coefs, terms, tbl] = fixed_names(meta, spec);
    if is_grp_perm
        % the group's coefficient names: '<group>_<category>' of the column as
        % the fit converts it (always categorical when group is in
        % spec.categorical, as core_snpm_lmm sets it); otherwise any coefficient
        % part '<group>_...' that is not the name of another column
        if iscategorical(tbl.(grp))
            grp_coef = strcat([grp '_'], categories(tbl.(grp)))';
        else
            parts = unique(split_parts(coefs));
            grp_coef = setdiff(parts(startsWith(parts, [grp '_'])), meta.Properties.VariableNames);
        end
        other_col = ismember(toks, meta.Properties.VariableNames) & ~strcmp(toks, grp);
        names_group = strcmp(toks, grp) | ismember(toks, grp_coef) | ...
            (startsWith(toks, [grp '_']) & ~other_col);   % group_C: left to unknownEffect
        if ~any(names_group)
            not_permuted(grp, effect);
        end
    end

    fixed = char(spec.fixed);
    if strcmpi(spec.effect_type, 'factor')
        want = sort(strtrim(strsplit(effect, ':')));
        hit = find(cellfun(@(t) isequal(sort(strsplit(t, ':')), want), terms), 1);
        if isempty(hit)
            error('core_snpm_lmm:unknownEffect', ...
                'Effect ''%s'' is not a term of the fixed part ''%s'' (terms: %s)', ...
                effect, fixed, strjoin(terms, ', '));
        end
        effect = terms{hit};
    elseif ~any(strcmp(coefs, effect)) && ~any(contains(coefs, effect))
        error('core_snpm_lmm:unknownEffect', ...
            'Effect ''%s'' is not a coefficient of the fixed part ''%s'' (coefficients: %s)', ...
            effect, fixed, strjoin(coefs, ', '));
    end
end

function not_permuted(grp, effect)
    error('core_snpm_lmm:effectNotPermuted', ...
        ['The group_label scheme relabels subjects between levels of ''%s'', but ' ...
         'the effect ''%s'' does not contain ''%s'', so every permuted map ' ...
         'would equal the observed one. Use within_subject for an effect that varies ' ...
         'within subjects, or test an effect of the group column.'], grp, effect, grp);
end

function p = split_parts(names)
    p = {};
    for k = 1:numel(names), p = [p, strsplit(names{k}, ':')]; end %#ok<AGROW>
end

function [coefs, terms, tbl] = fixed_names(meta, spec)
% Coefficient and term names of the fixed part, from one fitlm on the
% metadata prepared as snpm_lmm_fit prepares it, POWER appended last.
    tbl = snpm_lmm_categorize(meta, spec);
    n = height(tbl);
    tbl.POWER = sin((1:n)') + (1:n)' / max(n, 1);   % fixed, non-constant, no random draw
    dv = char(spec.dv);
    if strcmp(dv, 'POWER') || ~ismember(dv, tbl.Properties.VariableNames) || ~isnumeric(tbl.(dv))
        tbl.snpm_probe_response__ = cos(3 * (1:n)');
        dv = 'snpm_probe_response__';
    end
    ws = warning('off', 'all');
    restore = onCleanup(@() warning(ws));
    mdl = fitlm(tbl, sprintf('%s ~ %s', dv, char(spec.fixed)));
    coefs = mdl.CoefficientNames(:)';
    terms = mdl.Formula.TermNames(:)';
    terms = terms(~strcmp(terms, '(Intercept)'));
end
