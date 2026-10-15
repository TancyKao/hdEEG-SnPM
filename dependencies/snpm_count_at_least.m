function count = snpm_count_at_least(null, obs)
%SNPM_COUNT_AT_LEAST  Permutation count #{null >= obs}, ties within rounding counted.
%
%   count = SNPM_COUNT_AT_LEAST(null, obs)
%
%   null  : vector of null values (e.g. the max-statistic null), already
%           oriented for the tail (absolute values for a two-sided test).
%           NaN entries never count (a NaN comparison is false), as before.
%   obs   : scalar or vector of observed values, oriented the same way.
%   count : same size as obs; count(k) = #{null >= obs(k) - TIE_RTOL*|obs(k)|}.
%           The slack applies only where obs(k) is finite.
%
% Why (round-4 spec, approved 2026-10-14): the (b+1)/(N+1) estimator counts
% null values >= observed in exact arithmetic. In the GLM, LMM and circular
% tiers the observed statistic is computed on Y while an identity or
% equivalent labelling is computed on Zfit + R(P): equal mathematically,
% different in the last bits, with the rounding falling either way. A strict
% >= then drops about half the true ties and biases p DOWNWARD
% (anti-conservative). TIE_RTOL = 1e-12 leaves 25-50x headroom over the
% largest rounding ties seen (4e-14) and cross-language map differences
% (2.1e-14). A false tie between genuinely distinct values that close is
% negligible and errs conservative.
%
% Mirrors legacy._at_least in the Python port (TIE_RTOL there too). Used for
% the max-null FWE counts (legacy tMax/tMaxTFCE, snpm_perm_correction TFCE and
% cluster mass, LMM TFCE and cluster mass). Uncorrected p is parametric and is
% never counted against a null. Cluster-extent (integer size) counts in
% snpm_cluster_analysis stay strict.

    TIE_RTOL = 1e-12;
    null = reshape(null, [], 1);
    sz = size(obs);
    obs = reshape(obs, 1, []);
    slack = TIE_RTOL * abs(obs);
    slack(~isfinite(obs)) = 0;
    count = reshape(sum(null >= (obs - slack), 1), sz);
end
