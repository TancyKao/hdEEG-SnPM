function [out_png, used_channels] = plot_interaction_lines(power, meta, group_col, ...
        condition_col, subject_col, channels, out_png)
% Classic interaction line plot for a two-way mixed ANOVA. Averages power over
% CHANNELS (spatial mean per observation), then for each group x condition cell
% computes the mean and BETWEEN-SUBJECT standard error, and draws one line per
% group (x = condition levels in stable order, error bars = SE across subjects).
%
% Averaging over channels first, then aggregating to cell means with a
% subject-level SE, makes the error bars reflect between-subject variability
% (the relevant error for reading the interaction), not channel scatter.
%
% INPUTS
%   power        : nObs x nCh
%   meta         : table with the group / condition / subject columns
%   group_col, condition_col, subject_col : char column names
%   channels     : channel indices to average over (interaction cluster;
%                  empty -> all channels)
%   out_png      : full path of the PNG to write
% OUTPUTS
%   out_png        : the written path (echoed)
%   used_channels  : the channel indices actually averaged over

    nCh = size(power, 2);
    channels = channels(:).';
    channels = channels(channels >= 1 & channels <= nCh);
    if isempty(channels), channels = 1:nCh; end
    used_channels = channels;

    sig = mean(power(:, channels), 2, 'omitnan');       % nObs x 1

    gcol = meta.(group_col); ccol = meta.(condition_col); scol = meta.(subject_col);
    [glev, ~, gi] = unique(gcol, 'stable');
    [clev, ~, ci] = unique(ccol, 'stable');
    [slev, ~, si] = unique(scol, 'stable'); %#ok<ASGLU>
    nG = numel(glev); nC = numel(clev); nS = numel(unique(si));

    % each subject's (constant) group
    subj_g = zeros(nS, 1);
    for s = 1:nS, r = find(si == s, 1); subj_g(s) = gi(r); end

    mu = nan(nG, nC); se = nan(nG, nC);
    for g = 1:nG
        subs = find(subj_g == g);
        for c = 1:nC
            vals = nan(numel(subs), 1);
            for k = 1:numel(subs)
                rows = (si == subs(k)) & (ci == c);
                if any(rows), vals(k) = mean(sig(rows), 'omitnan'); end
            end
            vals = vals(~isnan(vals));
            if ~isempty(vals)
                mu(g, c) = mean(vals);
                if numel(vals) > 1, se(g, c) = std(vals) / sqrt(numel(vals)); else, se(g, c) = 0; end
            end
        end
    end

    condLabels = cell(1, nC); for c = 1:nC, condLabels{c} = lev2str(clev, c); end
    groupLabels = cell(1, nG); for g = 1:nG, groupLabels{g} = lev2str(glev, g); end

    fig = figure('Position', [50 50 560 420], 'Color', 'w', 'PaperPositionMode', 'auto');
    cl = onCleanup(@() close(fig)); %#ok<NASGU>
    colors = lines(max(nG, 1));
    x = 1:nC;
    h = gobjects(nG, 1);
    hold on;
    for g = 1:nG
        h(g) = errorbar(x, mu(g, :), se(g, :), '-o', 'Color', colors(g, :), ...
            'MarkerFaceColor', colors(g, :), 'LineWidth', 1.6, 'MarkerSize', 6, 'CapSize', 8);
    end
    hold off;
    set(gca, 'XTick', x, 'XTickLabel', condLabels);
    xlim([0.5, nC + 0.5]);
    xlabel(condition_col, 'Interpreter', 'none');
    ylabel('Mean power (selected channels)');
    legend(h, groupLabels, 'Location', 'best', 'Interpreter', 'none');
    title('Group x Condition interaction', 'Interpreter', 'none');
    box on; grid on;
    print(fig, '-dpng', '-r300', out_png);
end

function s = lev2str(levels, i)
    if iscell(levels), v = levels{i}; else, v = levels(i); end
    if isnumeric(v) || islogical(v), s = num2str(v); else, s = char(string(v)); end
end
