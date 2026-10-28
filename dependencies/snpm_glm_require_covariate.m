function snpm_glm_require_covariate(preset, covariate_cols, D)
% Refuse an ANCOVA that has no covariate: it would be a one-way ANOVA labelled
% "covariate-adjusted" (core_snpm_glm:ancovaNoCovariate).
%
%   snpm_glm_require_covariate(preset, covariate_cols)     % before the data are read
%   snpm_glm_require_covariate(preset, covariate_cols, D)  % after snpm_glm_design
%
% Without D: refuses when covariate_cols names no column once blank or
% whitespace-only entries are dropped (snpm_glm_covariate_list: [], {}, '', "",
% {''}, {' '}, [""] all count as empty; {'age',''} is {'age'}). With D: also refuses when the design
% kept no covariate column, i.e. every listed covariate was constant (or a
% one-level category) on these rows and snpm_glm_design dropped it.
% Any preset other than 'ancova' is a no-op.
%
% Message text is shared with the Python port; do not reword.

    if ~strcmpi(preset, 'ancova'), return, end
    if nargin >= 3 && ~isempty(D)
        none = ~any(startsWith(D.col_names, 'covariate '));
    else
        none = isempty(snpm_glm_covariate_list(covariate_cols));
    end
    if none
        error('core_snpm_glm:ancovaNoCovariate', ...
            'ANCOVA needs at least one covariate; with none, use anova1.');
    end
end
