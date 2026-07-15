function Yperm = snpm_glm_permute(flctx)
% One Freedman-Lane permutation realization (Winkler et al. 2014, NeuroImage).
% To test an effect while controlling nuisance covariates, the nuisance fit
% is held fixed and the *residuals* of the data on the nuisance are permuted,
% then added back:  Yperm = Z*betaZ + permute(residuals).
%
% flctx (built once by the caller) holds:
%   .Zfit      : nObs x nCh nuisance fitted values  Z*(Z\Y)
%   .R         : nObs x nCh residuals               Y - Zfit
%   .eb        : nObs x 1 exchangeability-block id
%   .perm_type : 'free'  -> permute rows across all observations (between)
%                'within'-> permute rows within each block (repeated measures)
%
% Returns one permuted data matrix Yperm (nObs x nCh).

    R = flctx.R;
    n = size(R, 1);

    if strcmp(flctx.perm_type, 'free')
        pe = randperm(n);
    else
        pe = within_block_perm(flctx.eb);
    end

    Yperm = flctx.Zfit + R(pe, :);
end

function pe = within_block_perm(eb)
% Permutation index that shuffles rows only within each exchangeability block.
    pe = (1:numel(eb))';
    blocks = unique(eb, 'stable');
    for b = 1:numel(blocks)
        rows = find(eb == blocks(b));
        if numel(rows) > 1
            pe(rows) = rows(randperm(numel(rows)));
        end
    end
end
