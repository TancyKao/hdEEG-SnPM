function p = watson_u2_pvalue(u)
% Asymptotic upper-tail p-value for Watson's two-sample U^2 statistic.
%
%   p = watson_u2_pvalue(u)
%
% Elementwise; p has the same size as u. NaN in -> NaN out.
%
% Uses the asymptotic series (Watson 1961; Zar, Biostatistical Analysis,
% Ch. 27), truncated at 12 terms:
%
%   P(U^2 >= u) = 2 * sum_{k=1}^{12} (-1)^(k-1) * exp(-2 * k^2 * pi^2 * u)
%
% WHY A SERIES AND NOT THE PUBLISHED TABLE. The published critical values are
% a five-point step function (alpha = .10 .05 .025 .01 .001). Printing a
% five-valued p per channel into the Excel table and HTML report reads as a
% broken engine, and it cannot be used as a cluster-forming threshold at an
% arbitrary alpha. The series reproduces all five published critical values
% (U^2 = 0.152 / 0.187 / 0.221 / 0.268 / 0.385) to better than 5e-4, which is
% far inside the resolution any permutation test built on top of it can
% deliver, and it is continuous. Note that the toolbox's *inference* is the
% permutation p, not this one; this value is used for cluster forming and for
% descriptive reporting.
%
% SMALL-u GUARD. The series is asymptotic and only converges usefully once
% u is not tiny: as u -> 0 every term tends to 1 and a 12-term alternating
% truncation oscillates instead of tending to the correct limit of 1. Below
% U2_SERIES_FLOOR the exact limit (p = 1) is returned instead. The floor is
% set at 0.0065, which is where the 12-term truncation still agrees with a
% 5000-term evaluation to 1e-9; the true tail probability there is
% 0.99999996, so the substitution is a 4e-8 step, not a visible
% discontinuity. For any realistic sample size the observed U^2 is far above
% the floor.
%
% The result is clamped to [0, 1] — the alternating series can overshoot
% slightly on the low-u side.
%
% See also SNPM_CIRC_WATSONS_U2.

    U2_SERIES_FLOOR = 0.0065;

    p = nan(size(u));
    finite = isfinite(u);
    if ~any(finite(:))
        return
    end

    uu = u(finite);
    uu = uu(:)';                       % 1 x m

    k  = (1:12)';                      % 12 x 1
    sgn = (-1).^(k - 1);
    series = 2 * sum(sgn .* exp(-2 * (k.^2) * (pi^2) * uu), 1);

    series(uu <= U2_SERIES_FLOOR) = 1; % includes u <= 0 (degenerate)
    series = min(max(series, 0), 1);

    p(finite) = series;
end
