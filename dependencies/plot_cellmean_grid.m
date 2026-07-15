function [pngs, labels, row_labels, col_labels] = plot_cellmean_grid(power, meta, ...
        group_col, condition_col, chanlocs, out_prefix)
% Group x condition cell-mean scalp topographies for a two-way mixed ANOVA.
% For each (group, condition) cell, averages power over that cell's rows to an
% nCh vector and renders a scalp topo on a SHARED colour scale (a common clim
% across all cells, so the cells are directly comparable by eye). One PNG per
% cell: "<out_prefix>_cell_<group>_<condition>.png".
%
% OUTPUTS
%   pngs        : cellstr of PNG basenames (filename+ext, for relative HTML refs)
%   labels      : parallel cellstr "Group / Condition"
%   row_labels  : group level labels (grid rows)
%   col_labels  : condition level labels (grid cols)

    gcol = meta.(group_col); ccol = meta.(condition_col);
    [glev, ~, gi] = unique(gcol, 'stable');
    [clev, ~, ci] = unique(ccol, 'stable');
    nG = numel(glev); nC = numel(clev);

    % First pass: all cell means -> shared colour scale (robust 2-98 pctile).
    means = cell(nG, nC);
    allvals = [];
    for g = 1:nG
        for c = 1:nC
            rows = (gi == g) & (ci == c);
            if any(rows), m = mean(power(rows, :), 1, 'omitnan'); else, m = nan(1, size(power, 2)); end
            means{g, c} = m;
            allvals = [allvals, m(isfinite(m))]; %#ok<AGROW>
        end
    end
    clim = [];
    if ~isempty(allvals)
        lo = prctile(allvals, 2); hi = prctile(allvals, 98);
        if isfinite(lo) && isfinite(hi) && hi > lo, clim = [lo hi]; end
    end

    row_labels = cell(1, nG); for g = 1:nG, row_labels{g} = lev2str(glev, g); end
    col_labels = cell(1, nC); for c = 1:nC, col_labels{c} = lev2str(clev, c); end

    pngs = {}; labels = {};
    for g = 1:nG
        for c = 1:nC
            lbl  = sprintf('%s / %s', row_labels{g}, col_labels{c});
            safe = matlab.lang.makeValidName(sprintf('%s_%s', row_labels{g}, col_labels{c}));
            png  = sprintf('%s_cell_%s.png', out_prefix, safe);
            try
                save_single_topo(means{g, c}, chanlocs, lbl, png, 'parula', clim);
                [~, fn, ext] = fileparts(png);
                pngs{end+1}   = [fn ext]; %#ok<AGROW>
                labels{end+1} = lbl;      %#ok<AGROW>
            catch ME
                warning(ME.identifier, 'cell-mean topo failed (%s): %s', lbl, ME.message);
            end
        end
    end
end

function s = lev2str(levels, i)
    if iscell(levels), v = levels{i}; else, v = levels(i); end
    if isnumeric(v) || islogical(v), s = num2str(v); else, s = char(string(v)); end
end
