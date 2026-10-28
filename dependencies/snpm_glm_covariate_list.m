function cols = snpm_glm_covariate_list(cols)
% Covariate column names as a 1 x k cellstr with blank or whitespace-only
% entries removed: {'age',''} -> {'age'}; [], {}, '', "", {''}, [""] -> {}.
% Shared by core_snpm_glm, snpm_glm_design, export_glm_report and
% snpm_glm_require_covariate, and matches the Python port.
    if isempty(cols), cols = {}; return, end
    cols = cellstr(string(cols));
    cols = strtrim(cols(:)');
    cols = cols(~cellfun(@isempty, cols));
end
