function w = snpm_lmm_within_ratio(V, subj)
% Within-subject share of the variation of each column of V.
%
%   w = snpm_lmm_within_ratio(V, subj)
%
%   V    : rows x k numeric (trials x channels, or one variable as a column)
%   subj : rows x 1 subject ID (any type unique() accepts)
%   w    : 1 x k,  w = SS_within / SS_total, with
%            SS_within = sum over rows of (v - mean of v within its subject)^2
%            SS_total  = sum over rows of (v - mean of v)^2
%          and w = 0 when SS_total <= 1e-20 * sum(v.^2) (a flat column, which is
%          constant within subject too). A column containing NaN gets NaN.
%
% Grouped by SUBJECT only, never by group or condition, so w does not depend
% on any permutation's labelling: core_snpm_lmm computes it ONCE, before the
% permutation loop, and the channels (or the run) it rules out are the same
% for the observed map and every permuted map. Threshold used there: 1e-6.

    V = double(V);
    [~, ~, si] = unique(subj, 'stable');
    si = si(:);
    q = max(si);
    cnt = accumarray(si, 1, [q, 1]);
    w = nan(1, size(V, 2));
    for j = 1:size(V, 2)
        v = V(:, j);
        if any(~isfinite(v)), continue; end
        m = accumarray(si, v, [q, 1]) ./ cnt;
        ssw = sum((v - m(si)).^2);
        sst = sum((v - mean(v)).^2);
        if sst <= 1e-20 * sum(v.^2)
            w(j) = 0;
        else
            w(j) = ssw / sst;
        end
    end
end
