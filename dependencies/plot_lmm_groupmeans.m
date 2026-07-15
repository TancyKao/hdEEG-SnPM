function levels = plot_lmm_groupmeans(y, group, ylab, out_png)
% Mean +/- SE bar chart of a response per group, for the LMM report's
% "categorical effect" descriptive. Y is summarised over the significant
% channels (computed by the caller); one bar per group level, error bars =
% standard error across the observations in that level.
%
% INPUTS
%   y      : nObs x 1 response (cluster-averaged power, or a DV column)
%   group  : nObs x 1 grouping (numeric/cellstr/categorical)
%   ylab   : y-axis label
%   out_png: full path of the PNG to write
% OUTPUT
%   levels : cellstr of the group level labels, in plot order
    y = y(:);
    [glev, ~, gi] = unique(group, 'stable'); nG = numel(glev);
    mu = nan(nG, 1); se = nan(nG, 1); levels = cell(1, nG);
    for g = 1:nG
        v = y(gi == g); v = v(isfinite(v));
        if ~isempty(v)
            mu(g) = mean(v);
            if numel(v) > 1, se(g) = std(v) / sqrt(numel(v)); else, se(g) = 0; end
        end
        levels{g} = lev2str(glev, g);
    end

    fig = figure('Position', [50 50 560 440], 'Color', 'w', 'PaperPositionMode', 'auto');
    cl = onCleanup(@() close(fig)); %#ok<NASGU>
    bar(1:nG, mu, 0.6, 'FaceColor', [0.12 0.43 0.55], 'EdgeColor', 'none'); hold on;
    errorbar(1:nG, mu, se, 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 10);
    hold off;
    set(gca, 'XTick', 1:nG, 'XTickLabel', levels);
    xlim([0.4, nG + 0.6]);
    ylabel(ylab, 'Interpreter', 'none');
    box on; grid on;
    print(fig, '-dpng', '-r300', out_png);
end

function s = lev2str(levels, i)
    if iscell(levels), v = levels{i}; else, v = levels(i); end
    if isnumeric(v) || islogical(v), s = num2str(v); else, s = char(string(v)); end
end
