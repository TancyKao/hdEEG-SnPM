function w = snpm_lmm_effect_within_share(meta, effect, subj_col)
% Largest within-subject share of variation among the metadata variables of a
% mixed-model effect, for the within_subject refusal
% core_snpm_lmm:effectConstantWithinSubject.
%
%   w = snpm_lmm_effect_within_share(meta, effect, subj_col)
%
% effect is split on ':' and '*'. Returns NaN ("nothing to decide", the run
% goes on) when any token is POWER (POWER's within-subject variation is per
% channel and handled by the channel exclusion) or when no token resolves to a
% metadata column. A token resolves to the column of that exact name, else to
% the column c for which the token starts with [c '_'] (a coefficient name such
% as group_B -> group; the longest such c when several match). Tokens that
% resolve to nothing are ignored.
%
% Per resolved column, w = snpm_lmm_within_ratio grouped by subject on the
% analysed rows: a numeric (or logical) column directly; any other column
% (categorical, string, cellstr) through one 0/1 indicator per level, taking
% the largest w over the indicators. The returned w is the largest over the
% resolved columns; the caller refuses when it is <= 1e-6.

    toks = strtrim(regexp(char(effect), '[:*]', 'split'));
    w = NaN;
    if any(strcmp(toks, 'POWER')), return; end
    vn = meta.Properties.VariableNames;
    cols = {};
    for t = toks
        if ismember(t{1}, vn)
            cols{end+1} = t{1}; %#ok<AGROW>
        else
            c = vn(startsWith(t{1}, strcat(vn, '_')));
            if ~isempty(c)
                [~, i] = max(cellfun(@numel, c));
                cols{end+1} = c{i}; %#ok<AGROW>
            end
        end
    end
    cols = unique(cols, 'stable');
    if isempty(cols), return; end
    subj = meta.(subj_col);
    ws = zeros(1, numel(cols));
    for k = 1:numel(cols)
        x = meta.(cols{k});
        if isnumeric(x) || islogical(x)
            ws(k) = snpm_lmm_within_ratio(double(x), subj);
        else
            [~, ~, li] = unique(cellstr(string(x)));
            D = double(li(:) == 1:max(li));
            ws(k) = max(snpm_lmm_within_ratio(D, subj));
        end
    end
    w = max(ws);
end
