function snpm_glm_require_numeric_predictor(preset, meta, col)
% Refuse a regression whose predictor column is not numeric
% (core_snpm_glm:predictorNotNumeric), before the design is built; centring a
% text column would otherwise fail with an unrelated MATLAB error.
%
%   snpm_glm_require_numeric_predictor(preset, meta, col)
%
% Only numeric (non-logical) columns pass. Logical columns (e.g. TRUE/FALSE
% cells of an Excel sheet, which readtable reads as logical), text, string and
% categorical columns are refused, as the Python port refuses bool and text
% columns. Text is refused even when every value reads as a number: readtable
% parses numeric columns itself, so a column that arrives as text is text.
% Any preset other than 'regression' is a no-op. The quoted value is the
% column's first value that is neither missing nor an empty string ('' if
% there is none); a logical value is quoted in lowercase, 'true'/'false'
% (string(true) is "true"). Message text is shared with the Python port; do
% not reword.
    if ~strcmpi(preset, 'regression'), return, end
    v = meta.(col);
    if isnumeric(v), return, end   % logical is not numeric in MATLAB: refused
    s = string(v);
    % the first NON-MISSING value, as the Python port quotes it: an empty text
    % cell is missing (readtable reads an empty CSV cell as '', pandas as NaN),
    % as are <missing> strings and <undefined> categories
    s = s(~ismissing(s) & strlength(s) > 0);
    first = '';
    if ~isempty(s), first = char(s(1)); end
    error('core_snpm_glm:predictorNotNumeric', '%s', sprintf(['The predictor column ''%s'' is not ' ...
        'numeric (it holds text such as ''%s''); regression needs a numeric predictor. To compare ' ...
        'groups, use anova1 or ancova.'], col, first));
end
