function [t, degen] = snpm_ttest_novar(t, a, b, design)
%SNPM_TTEST_NOVAR  Zero-variance rule 2 for the legacy t-tests.
%
%   [t, degen] = SNPM_TTEST_NOVAR(t, a, b, design)
%
%   t      : 1 x nCh t map as returned by ttest / ttest2 for this labelling
%   a, b   : the two matrices that t was computed from, for THIS labelling
%            design 'paired'   -> a, b row-matched (nSubj x nCh), d = a - b
%            design 'unpaired' -> a = group-1 rows, b = group-2 rows
%   t      : (out) t with every degenerate channel set to 0
%   degen  : 1 x nCh logical, true where the error sum of squares of this
%            labelling counts as zero (snpm_zero_ss)
%
% Rule (Method Spec 2026-10-14, t-test tier rule 2): at any labelling where
% the error sum of squares meets the tolerance the t statistic is DEFINED as
% 0. Paired / one-sample: SS = sum((d - mean d).^2) against sum(d.^2).
% Unpaired: the pooled within-group SS against sum(x.^2) over both groups.
% Applied identically to the observed labelling and to every permutation, so
% it is one fixed statistic and the permutation test stays exact. Without it
% ttest returns +/-Inf (all differences equal) and ttest2 returns +/-Inf or a
% finite 1e16-scale t (0.1s vs 0.3s), and ClusterEnhancement then tries to
% build a threshold ladder dh:dh:max+dh that does not fit in memory.
%
% Computed on the usable (finite) values of each channel, the same values
% ttest/ttest2 use. A channel with too few usable values for a t at all
% (paired n < 2; unpaired a group empty or n1 + n2 < 3) is never marked here.
% Its builtin t is NOT necessarily NaN: paired/one-sample with n = 1 gives
% std = 0 and df = 0, so ttest returns +/-Inf (NaN only when the single
% difference is 0); ttest2 with an empty group returns NaN. Such channels
% never reach the engines: the missing-data degeneracy guard in
% core_snpm_analysis (m >= 2 usable pairs; m >= max(n1,n2)+2 for unpairedT)
% has excluded them before any statistic is computed.
%
% On data where no labelling is degenerate this changes nothing: t is only
% written where degen is true.

    switch design
        case 'paired'
            d  = a - b;
            ok = isfinite(d);
            n  = sum(ok, 1);
            d(~ok) = 0;
            m  = sum(d, 1) ./ max(n, 1);
            r  = (d - m) .* ok;
            ss = sum(r.^2, 1);
            v2 = sum(d.^2, 1);
            degen = snpm_zero_ss(ss, v2) & n >= 2;
        case 'unpaired'
            [ssa, v2a, na] = within_ss(a);
            [ssb, v2b, nb] = within_ss(b);
            degen = snpm_zero_ss(ssa + ssb, v2a + v2b) & na >= 1 & nb >= 1 & (na + nb) >= 3;
        otherwise
            error('snpm_ttest_novar:badDesign', 'design must be ''paired'' or ''unpaired''.');
    end
    degen = reshape(degen, 1, []);
    t(degen) = 0;
end

function [ss, v2, n] = within_ss(x)
    ok = isfinite(x);
    n  = sum(ok, 1);
    x(~ok) = 0;
    m  = sum(x, 1) ./ max(n, 1);
    r  = (x - m) .* ok;
    ss = sum(r.^2, 1);
    v2 = sum(x.^2, 1);
end
