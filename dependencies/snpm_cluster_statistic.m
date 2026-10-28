function [out, line] = snpm_cluster_statistic(varargin)
% The cluster statistic of the GLM, circular and mixed-model tiers: one place
% for the setting, its validation, its wording and the scoring rule.
%
% VALIDATE (called once by core_snpm_glm, core_snpm_circ, core_snpm_lmm and
% export_glm_report):
%   [cs, line] = snpm_cluster_statistic(params)   % struct: reads .cluster_statistic
%   [cs, line] = snpm_cluster_statistic(value)    % char/string, or [] for the default
%   [cs, line] = snpm_cluster_statistic(x, what)  % what = the summed quantity
%     cs   : 'sum' (default) | 'mean'
%     line : the one results-text line naming the rule. what defaults to
%            't² (or F)' (GLM, mixed model); the circular tier passes
%            'the statistic' (its Hotelling F, Watson U² or circular-linear F
%            is summed as is, not squared). Texts are shared with the Python
%            port; do not reword.
%   Any other value is refused with snpm:badClusterStatistic.
%
% SCORE (called by snpm_perm_correction/find_pclusters,
% snpm_lmm_cluster/local_find_pclusters and the parity replay in
% scripts/export_parity_fixtures.m):
%   score = snpm_cluster_statistic(wald, ch, cs)
%     wald : 1 x nCh Wald map (t^2 or F; U^2 or Hotelling F on the circular tier)
%     ch   : channel indices of one cluster
%     score = sum(wald(ch),'omitnan'), divided by numel(ch) under 'mean'.
%
% 'sum' is a summed-statistic cluster mass in the sense of Maris & Oostenveld
% (2007), taken here on the Wald scale (t^2 or F; they and FSL sum t, which
% does not affect validity): it grows with extent, so a broad moderate effect is not beaten by a
% single strong channel in the permutation maxima. 'mean' is the rule used
% before 2026-10-28 (Stephan et al. 2021's recipe); it reproduces older runs.
% The threshold, the Wald scale, the permutation draws and TFCE do not depend
% on this setting. The test stays exact under either rule: the same rule is
% applied to the observed map and to every permuted map.

    if nargin == 3
        wald = varargin{1};
        ch   = varargin{2};
        out  = sum(wald(ch), 'omitnan');
        if strcmp(varargin{3}, 'mean')
            out = out / numel(ch);
        end
        line = '';
        return
    end
    if nargin > 3
        error('snpm_cluster_statistic:nargin', ...
            'Call as snpm_cluster_statistic(params|value[, what]) or snpm_cluster_statistic(wald, ch, cs).');
    end
    what = 't² (or F)';
    if nargin == 2 && ~isempty(varargin{2}), what = char(varargin{2}); end

    v = varargin{1};
    if isstruct(v)
        if isfield(v, 'cluster_statistic')
            v = v.cluster_statistic;
        else
            v = [];
        end
    end
    if isempty(v)
        out = 'sum';
    elseif (ischar(v) && (isrow(v))) || (isstring(v) && isscalar(v))
        out = lower(strtrim(char(v)));
    else
        out = '<non-text value>';
    end
    if isempty(out), out = 'sum'; end   % "" or blanks: same as absent
    if ~any(strcmp(out, {'sum', 'mean'}))
        if ischar(v) || isstring(v), shown = char(v); else, shown = out; end
        error('snpm:badClusterStatistic', ...
            'cluster_statistic must be ''sum'' or ''mean''; got ''%s''.', shown);
    end
    if strcmp(out, 'sum')
        line = ['Cluster statistic: sum of ' what ' over the cluster'];
    else
        line = ['Cluster statistic: sum of ' what ' over the cluster, divided by the ' ...
                'number of channels (mean; reproduces runs before this change)'];
    end
end
