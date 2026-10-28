function snpm_write_settings_sheet(outputSname, cluster_statistic, alpha)
% One-row 'settings' sheet in <outputSname>.xlsx recording how the clusters of
% a GLM, circular or mixed-model run were scored, so the workbook alone says
% which rule produced its cluster masses and p-values. Called by
% core_snpm_glm, core_snpm_circ and core_snpm_lmm right after
% func_genSnpmTable; the legacy t-test workbook does not get one.
%
%   cluster_statistic : 'sum' | 'mean' (snpm_cluster_statistic)
%   alpha             : the cluster-forming and significance level
    T = table({char(cluster_statistic)}, double(alpha), ...
        'VariableNames', {'cluster_statistic', 'alpha'});
    writetable(T, [outputSname '.xlsx'], 'Sheet', 'settings');
end
