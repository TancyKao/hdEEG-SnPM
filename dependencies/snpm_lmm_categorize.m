function base_meta = snpm_lmm_categorize(meta, spec)
% The metadata table exactly as every mixed-model fit sees it: the columns in
% spec.categorical (default: Subject and group, when present) converted to
% categorical once; every other column left as read. Shared by snpm_lmm_fit
% (the fits) and snpm_lmm_check_effect (the coefficient names), so the names
% checked are the names fitted.
%
% Level order, and so the reference level fitlme drops, follows from this:
%   converted column     categories(categorical(x)): text sorted, numbers in
%                        numeric order, logical 'false','true'; names
%                        '<var>_<category>' (group_B, group_2, group_true)
%   cellstr / string     not converted: levels in order of FIRST APPEARANCE in
%                        the fitted rows, the first one is the reference
%   categorical as read  its own categories() order
%   logical, not converted  reference false, name '<var>_1'
%   numeric, not converted  a continuous predictor, name '<var>'

    if isfield(spec, 'categorical') && ~isempty(spec.categorical)
        cat_vars = cellstr(spec.categorical);
    else
        cat_vars = intersect({'Subject', 'group'}, meta.Properties.VariableNames);
    end
    base_meta = meta;
    for v = 1:numel(cat_vars)
        cv = cat_vars{v};
        if ismember(cv, base_meta.Properties.VariableNames) && ~iscategorical(base_meta.(cv))
            base_meta.(cv) = categorical(base_meta.(cv));
        end
    end
end
