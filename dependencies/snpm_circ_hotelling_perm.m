function ctx = snpm_circ_hotelling_perm(gi, COV, usecov, valid)
%SNPM_CIRC_HOTELLING_PERM  Permutation-side Hotelling evaluator on a (cos,sin)
%   response, with a PER-CHANNEL nuisance covariate.
%
%   ctx = snpm_circ_hotelling_perm(gi, COV, usecov, valid)
%   [F, p, df2] = ctx.eval(Yc, Ys)
%
% WHY THIS EXISTS AND WHY IT IS NOT snpm_circ_hotelling.
% snpm_circ_hotelling takes ANGLES and embeds them itself. Under Freedman-Lane
% the permuted response is Zfit + R(pe,:) on the (cos,sin) embedding, which is
% NOT on the unit circle any more (that is the whole point: the nuisance fit is
% held fixed and only the residuals move). There is no angle that encodes it,
% so the permuted statistic has to be evaluated on the embedding directly.
% The algebra is identical to snpm_circ_hotelling's covariate-adjusted form and
% is pinned to it: on the OBSERVED embedding, ctx.eval(cos A, sin A) reproduces
% snpm_circ_hotelling(A, g, Cnuis) to machine precision (asserted by the
% caller's verification, see core_snpm_circ).
%
%   X = [1, cov_ch, g]   (or [1, g] where the covariate was dropped)
%   c = the g column;  Z = [yc, ys]
%   nu = N - rank(X);  B = (X'X)^-1 X' Z;  S = (Z - XB)'(Z - XB)/nu
%   T2 = (c'B) S^-1 (c'B)' / (c'(X'X)^-1 c);  F = T2 (nu-1)/(2 nu) ~ F(2, nu-1)
%
% X'X does not depend on the response, so (X'X)^-1 and c'(X'X)^-1 c are built
% ONCE here, outside the permutation loop, per channel. Only X'Z changes per
% permutation, and that is three vectorised reductions over subjects.
%
% TWO REGIMES, matching what the caller can guarantee:
%   usecov(ch) TRUE  - a per-channel covariate is in the design. The caller
%                      must have ensured the channel is COMPLETE (no missing
%                      subject), because Freedman-Lane fits the nuisance over
%                      whole columns; that is the same rule
%                      snpm_glm_fl_context enforces for the whole GLM tier.
%                      Solved by the vectorised 3-column path.
%   usecov(ch) FALSE - no covariate. The nuisance is the intercept alone, so
%                      Freedman-Lane reduces EXACTLY to permuting rows, missing
%                      cells ride along, and the classic two-sample Hotelling
%                      (pooled within-group covariance) is used with per-channel
%                      pairwise deletion.
%
% INPUTS
%   gi     : nSubj x 1 group indicator in {0,1} (no NaN)
%   COV    : nSubj x nCh per-channel covariate (e.g. log event count); may be []
%   usecov : 1 x nCh logical, TRUE where COV enters the design at that channel
%   valid  : 1 x nCh logical estimability mask, computed ONCE from the observed
%            labels by the caller and held FIXED across every permutation.
%            Channels outside it are returned NaN in every call, so the observed
%            map and every permuted map are defined on one channel set.
%
% OUTPUT ctx
%   .eval   function handle, [F, p, df2] = eval(Yc, Ys)
%   .valid  the mask, echoed back so the caller can assert it never changed
%   .nu     1 x nCh residual df of the full model (NaN where not evaluable)
%
% See also SNPM_CIRC_HOTELLING, SNPM_GLM_PERMUTE, SNPM_PERM_CORRECTION.

    gi = double(gi(:));
    nSubj = numel(gi);
    if isempty(COV), COV = zeros(nSubj, numel(usecov)); end
    nCh = size(COV, 2);
    usecov = reshape(logical(usecov), 1, nCh);
    valid  = reshape(logical(valid),  1, nCh);

    covCh = find(valid &  usecov);
    plnCh = find(valid & ~usecov);

    % ---- precompute (X'X)^-1 and c'(X'X)^-1 c for the covariate channels ----
    % X = [1, cov, g], complete rows, so X'X is fixed across permutations.
    nC = numel(covCh);
    XtXinv = zeros(3, 3, max(nC, 1));
    vgg    = zeros(1, nCh);
    okcov  = true(1, nC);
    s1 = nSubj; sg = sum(gi); sgg = sum(gi .^ 2);
    for k = 1:nC
        c = COV(:, covCh(k));
        M = [s1, sum(c), sg; sum(c), c' * c, c' * gi; sg, gi' * c, sgg];
        if rcond(M) < 1e-12
            okcov(k) = false;
            continue
        end
        Mi = M \ eye(3);
        XtXinv(:, :, k) = Mi;
        vgg(covCh(k))   = Mi(3, 3);
    end
    covCh  = covCh(okcov);
    XtXinv = XtXinv(:, :, okcov);
    COVc   = COV(:, covCh);
    nuCov  = nSubj - 3;

    nu = nan(1, nCh);
    nu(covCh) = nuCov;

    ctx = struct();
    ctx.valid = valid;
    ctx.nu    = nu;
    ctx.eval  = @evaluate;

    % ==================================================================
    function [F, p, df2] = evaluate(Yc, Ys)
        F   = nan(1, nCh);
        p   = nan(1, nCh);
        df2 = nan(1, nCh);

        % ---- covariate channels: vectorised 3-column solve --------------
        if ~isempty(covCh) && nuCov >= 2
            yc = Yc(:, covCh); ys = Ys(:, covCh);
            % X'y for each channel: [sum(y); cov'y; g'y]  (3 x 1 x m)
            rhsC = cat(1, reshape(sum(yc, 1),        1, 1, []), ...
                          reshape(sum(COVc .* yc, 1), 1, 1, []), ...
                          reshape(gi' * yc,           1, 1, []));
            rhsS = cat(1, reshape(sum(ys, 1),        1, 1, []), ...
                          reshape(sum(COVc .* ys, 1), 1, 1, []), ...
                          reshape(gi' * ys,           1, 1, []));
            Bc = reshape(pagemtimes(XtXinv, rhsC), 3, []);   % 3 x m
            Bs = reshape(pagemtimes(XtXinv, rhsS), 3, []);

            fitC = Bc(1, :) + COVc .* Bc(2, :) + gi * Bc(3, :);
            fitS = Bs(1, :) + COVc .* Bs(2, :) + gi * Bs(3, :);
            rc = yc - fitC; rs = ys - fitS;

            Scc = sum(rc .^ 2, 1) / nuCov;
            Sss = sum(rs .^ 2, 1) / nuCov;
            Scs = sum(rc .* rs, 1) / nuCov;
            dc  = Bc(3, :); ds = Bs(3, :);
            v   = vgg(covCh);

            detS = Scc .* Sss - Scs .^ 2;
            quad = (dc .^ 2 .* Sss - 2 * dc .* ds .* Scs + ds .^ 2 .* Scc) ./ detS;
            t2   = quad ./ v;
            f    = t2 * (nuCov - 1) / (2 * nuCov);
            f(~(detS > 1e-12)) = NaN;      % catches NaN too

            F(covCh)   = f;
            df2(covCh) = nuCov - 1;
        end

        % ---- plain channels: classic two-sample, pairwise deletion ------
        if ~isempty(plnCh)
            yc = Yc(:, plnCh); ys = Ys(:, plnCh);
            ok = isfinite(yc) & isfinite(ys);
            m1 = ok & (gi == 1); m2 = ok & (gi == 0);
            n1 = sum(m1, 1);     n2 = sum(m2, 1);
            N  = n1 + n2;        nuP = N - 2;

            yc0 = yc; yc0(~ok) = 0;
            ys0 = ys; ys0(~ok) = 0;
            s1c = sum(yc0 .* m1, 1) ./ max(n1, 1); s2c = sum(yc0 .* m2, 1) ./ max(n2, 1);
            s1s = sum(ys0 .* m1, 1) ./ max(n1, 1); s2s = sum(ys0 .* m2, 1) ./ max(n2, 1);

            rcm = yc0 - (m1 .* s1c + m2 .* s2c); rcm(~ok) = 0;
            rsm = ys0 - (m1 .* s1s + m2 .* s2s); rsm(~ok) = 0;
            Scc = sum(rcm .^ 2, 1) ./ nuP;
            Sss = sum(rsm .^ 2, 1) ./ nuP;
            Scs = sum(rcm .* rsm, 1) ./ nuP;

            dc = s1c - s2c; ds = s1s - s2s;
            detS = Scc .* Sss - Scs .^ 2;
            quad = (dc .^ 2 .* Sss - 2 * dc .* ds .* Scs + ds .^ 2 .* Scc) ./ detS;
            t2   = quad .* (n1 .* n2 ./ N);
            f    = t2 .* (nuP - 1) ./ (2 * nuP);
            bad  = ~(detS > 1e-12) | (n1 < 1) | (n2 < 1) | (nuP < 2);
            f(bad) = NaN;

            F(plnCh)   = f;
            df2(plnCh) = nuP - 1;
            df2(plnCh(bad)) = NaN;
        end

        good = isfinite(F) & isfinite(df2) & df2 > 0;
        p(good)   = fcdf(F(good), 2, df2(good), 'upper');
        F(~good)  = NaN;
        df2(~good) = NaN;
    end
end
