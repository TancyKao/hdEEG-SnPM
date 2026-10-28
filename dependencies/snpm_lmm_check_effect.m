function effect = snpm_lmm_check_effect(meta, spec)
% Refuse a mixed-model effect that the permutation cannot test or the model
% does not contain, or a fixed part the metadata cannot estimate, before any
% channel is fitted. Returns the effect as fitlme names it.
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
% fitlm is used rather than fitlme because it also names (and returns the
% design of) a rank-deficient fixed part, which fitlme refuses.
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
%    (core_snpm_lmm:unknownEffect). A continuous effect e resolves to ONE
%    coefficient by the first rule that applies (parts = e split on ':'):
%      a. e is a coefficient name -> that coefficient;
%      b. the parts, as a sorted list, equal the sorted ':'-parts of exactly
%         one coefficient -> that coefficient ('POWER:group_B' -> 'group_B:POWER');
%      c. the parts, sorted, equal the sorted variables of exactly one term:
%         a term with one coefficient -> that coefficient ('group' -> 'group_B'
%         for a two-level group, 'POWER:group' -> 'group_B:POWER'); a term with
%         several coefficients -> core_snpm_lmm:ambiguousEffect (a 3-level
%         'group' is two coefficients, not one);
%      d. otherwise core_snpm_lmm:unknownEffect.
%    A coefficient's term: each ':'-part maps to its variable (the part itself
%    if it is a column, else the longest column c with the part starting with
%    c_), compared as a sorted list. Reference levels (group_A) and '*' are not
%    coefficients. There is no substring matching: 'time' never finds
%    'sleeptime', 'POWER' never finds 'group_B:POWER'. Factor matching ignores
%    the order of the variables inside an interaction ('POWER:group' finds
%    'group:POWER'). Either way the RETURNED effect is the coefficient or term
%    exactly as fitlme names it, and snpm_lmm_fit matches it exactly.
% 3. meta given: the fixed part must be estimable from the metadata
%    (core_snpm_lmm:fixedRankDeficient). X is the fixed design of the same
%    fitlm, columns in coefficient order, POWER at its placeholder value.
%    Deficient when rank(X, tol) < size(X, 2) with
%    tol = max(size(X)) * 2.220446049250313e-16 * s_max(X) (numpy's default;
%    at least as loose as fitlme's own rank test, so every design fitlme
%    refuses at a generic channel is refused here; no column scaling, because
%    fitlme does not scale). fitlm itself is no flag: it zeroes aliased
%    coefficients and picks them by pivoted QR. The aliased columns are named
%    by a left-to-right walk: column j is kept when it raises the rank of the
%    kept columns before it, else it is aliased, and its message lists the kept
%    coefficients with |w_k| * norm(x_k) / norm(x_j) > 1e-6 in the least-squares
%    fit of x_j on them ("is zero on every trial" when none). The POWER
%    placeholder is generic, so a deficiency seen here holds at every channel
%    (every map would be NaN); aliasing that exists only for particular channel
%    values (a channel equal to a covariate) is not seen and stays a per-channel
%    NaN. Checked on the observed labels under either scheme.
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

    [coefs, terms, tbl, X] = fixed_names(meta, spec);
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
    else
        effect = resolve_coefficient(effect, coefs, terms, tbl.Properties.VariableNames, fixed);
    end

    check_rank(X, coefs, fixed);
end

function out = resolve_coefficient(effect, coefs, terms, vars, fixed)
% One coefficient for a continuous effect: exact name, else the same parts in
% another order, else the variables of a one-coefficient term.
    if any(strcmp(coefs, effect))
        out = effect;
        return
    end
    want = sort(strtrim(strsplit(effect, ':')));
    hit = find(cellfun(@(c) isequal(sort(strsplit(c, ':')), want), coefs));
    if isscalar(hit)
        out = coefs{hit};
        return
    end
    t = find(cellfun(@(s) isequal(sort(strsplit(s, ':')), want), terms));
    if isscalar(t)
        cs = coefs(coef_terms(coefs, terms, vars) == t);
        if isscalar(cs)
            out = cs{1};
            return
        elseif numel(cs) > 1
            error('core_snpm_lmm:ambiguousEffect', ...
                ['Effect ''%s'' names the term ''%s'' of the fixed part ''%s'', which has %d coefficients ' ...
                 '(%s); a continuous effect maps one coefficient. Name one of them, or set ' ...
                 'lmm_effect_type = ''factor'' for the omnibus F.'], effect, terms{t}, fixed, numel(cs), strjoin(cs, ', '));
        end
    end
    error('core_snpm_lmm:unknownEffect', ...
        'Effect ''%s'' is not a coefficient of the fixed part ''%s'' (coefficients: %s)', ...
        effect, fixed, strjoin(coefs, ', '));
end

function ct = coef_terms(coefs, terms, vars)
% Index into terms of each coefficient (0 for the intercept / no term): each
% ':'-part maps to its variable -- the part itself if it is a column, else the
% longest column c with the part starting with c_ -- compared as a sorted list.
    ct = zeros(1, numel(coefs));
    for i = 1:numel(coefs)
        p = strsplit(coefs{i}, ':');
        v = p;
        for j = 1:numel(p)
            if ismember(p{j}, vars), continue, end
            c = vars(cellfun(@(v) startsWith(p{j}, [v '_']), vars));
            if ~isempty(c)
                [~, m] = max(cellfun(@numel, c));
                v{j} = c{m};
            end
        end
        t = find(cellfun(@(s) isequal(sort(strsplit(s, ':')), sort(v)), terms), 1);
        if ~isempty(t), ct(i) = t; end
    end
end

function check_rank(X, coefs, fixed)
% core_snpm_lmm:fixedRankDeficient when the fixed design is singular. Aliased
% columns are found left to right (not by pivoted QR, which is neither stable
% nor portable) and named with the kept columns they are a combination of.
    s = svd(X);
    if isempty(s), return, end
    tol = max(size(X)) * 2.220446049250313e-16 * s(1);
    if rank(X, tol) == size(X, 2), return, end
    kept = [];
    items = {};
    for j = 1:size(X, 2)
        if rank(X(:, [kept j]), tol) > numel(kept)
            kept(end+1) = j; %#ok<AGROW>
            continue
        end
        xj = X(:, j);
        of = {};
        if ~isempty(kept) && norm(xj) > 0
            w = X(:, kept) \ xj;
            rel = abs(w(:))' .* vecnorm(X(:, kept)) / norm(xj);
            of = coefs(kept(rel > 1e-6));
        end
        if isempty(of)
            items{end+1} = sprintf('%s is zero on every trial', coefs{j}); %#ok<AGROW>
        else
            items{end+1} = sprintf('%s is a linear combination of %s', coefs{j}, strjoin(of, ', ')); %#ok<AGROW>
        end
    end
    error('core_snpm_lmm:fixedRankDeficient', ...
        ['The fixed part ''%s'' cannot be estimated from these trials: %s. fitlme refuses a ' ...
         'rank-deficient design, so every channel would be NaN. Remove or recode the redundant ' ...
         'column, or drop the term.'], fixed, strjoin(items, '; '));
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

function [coefs, terms, tbl, X] = fixed_names(meta, spec)
% Coefficient and term names of the fixed part, and its design X (columns in
% coefficient order, rows fitlm used), from one fitlm on the metadata prepared
% as snpm_lmm_fit prepares it, POWER appended last.
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
    X = mdl.Design(mdl.ObservationInfo.Subset, :);
end
