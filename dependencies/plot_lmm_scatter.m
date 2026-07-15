function used = plot_lmm_scatter(x, y, group, xlab, ylab, out_png)
% Scatter of y vs x with a least-squares fit line, for the LMM report's
% "continuous effect" descriptive. Each point is a trial. If GROUP is
% non-empty, points are coloured by group and one fit line is drawn per group
% (the moderation / interaction view); otherwise a single overall fit line.
%
% INPUTS
%   x, y   : nTrials x 1 (predictor, response). NaNs dropped pairwise.
%   group  : nTrials x 1 grouping (numeric/cellstr/categorical) or [] for none
%   xlab, ylab : axis labels
%   out_png    : full path of the PNG to write
% OUTPUT
%   used   : number of points actually plotted
    x = x(:); y = y(:);
    ok = isfinite(x) & isfinite(y);
    if ~isempty(group), group = group(ok); end
    x = x(ok); y = y(ok);
    used = numel(x);

    fig = figure('Position', [50 50 640 460], 'Color', 'w', 'PaperPositionMode', 'auto');
    cl = onCleanup(@() close(fig)); %#ok<NASGU>
    hold on;
    if isempty(group)
        scatter(x, y, 26, [0.31 0.43 0.69], 'filled', 'MarkerFaceAlpha', 0.5);
        if used >= 2
            b = polyfit(x, y, 1); xs = linspace(min(x), max(x), 100);
            plot(xs, polyval(b, xs), '-', 'Color', [0.76 0.33 0.23], 'LineWidth', 2);
        end
    else
        [glev, ~, gi] = unique(group, 'stable'); nG = numel(glev);
        colors = lines(max(nG, 1));
        h = gobjects(nG, 1); labs = cell(nG, 1);
        for g = 1:nG
            m = (gi == g);
            scatter(x(m), y(m), 26, colors(g, :), 'filled', 'MarkerFaceAlpha', 0.5);
            if sum(m) >= 2
                b = polyfit(x(m), y(m), 1); xs = linspace(min(x(m)), max(x(m)), 100);
                h(g) = plot(xs, polyval(b, xs), '-', 'Color', colors(g, :), 'LineWidth', 2);
            else
                h(g) = plot(nan, nan, '-', 'Color', colors(g, :));
            end
            labs{g} = lev2str(glev, g);
        end
        legend(h, labs, 'Location', 'best', 'Interpreter', 'none');
    end
    hold off;
    xlabel(xlab, 'Interpreter', 'none'); ylabel(ylab, 'Interpreter', 'none');
    box on; grid on;
    print(fig, '-dpng', '-r300', out_png);
end

function s = lev2str(levels, i)
    if iscell(levels), v = levels{i}; else, v = levels(i); end
    if isnumeric(v) || islogical(v), s = num2str(v); else, s = char(string(v)); end
end
