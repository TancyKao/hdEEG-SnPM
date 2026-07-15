function plot_topoInd(dat1, dat2, chanlocs, dat1name, dat2name, filepath, figtitle, varargin)
    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'cmap1', 'hot', @(x) ischar(x) || iscell(x) || isnumeric(x));
    addParameter(p, 'cmap2', 'hot', @(x) ischar(x) || iscell(x) || isnumeric(x));
    addParameter(p, 'cmap_diff', 'jet', @(x) ischar(x) || iscell(x) || isnumeric(x));
    addParameter(p, 'clim1', 'auto', @(x) strcmp(x,'auto') || (isnumeric(x) && length(x)==2));
    addParameter(p, 'clim2', 'auto', @(x) strcmp(x,'auto') || (isnumeric(x) && length(x)==2));
    addParameter(p, 'clim_diff', 'auto', @(x) strcmp(x,'auto') || (isnumeric(x) && length(x)==2));
    % diffdata: data for the 3rd (effect) map. Default [] -> dat1-dat2 (mean
    % difference). Provide it to plot a different statistic there, e.g. the
    % per-channel correlation r-map for a correlation analysis.
    addParameter(p, 'diffdata', [], @(x) isempty(x) || isnumeric(x));
    % difftitle: plain title for the 3rd map (filename still from figtitle).
    addParameter(p, 'difftitle', '', @(x) ischar(x) || isstring(x));
    % skip2: don't render the 2nd (dat2) map. Used for one-sample-vs-0, where
    % dat2 is an all-zero placeholder (constant data breaks topoplot's minmax).
    addParameter(p, 'skip2', false, @(x) islogical(x) || isnumeric(x));
    parse(p, varargin{:});
    
    % Common plot settings
    plot_settings = {'headrad', 0.6, 'style', 'map', 'electrodes', 'on', 'maplimits', 'minmax', 'whitebk', 'on'};
    
    % Create plots with individual colormaps
    create_topo_plot(dat1, chanlocs, dat1name, filepath, p.Results.cmap1, p.Results.clim1, plot_settings, true);
    if ~p.Results.skip2
        create_topo_plot(dat2, chanlocs, dat2name, filepath, p.Results.cmap2, p.Results.clim2, plot_settings, true);
    end
    if isempty(p.Results.diffdata), diffmap = dat1 - dat2; else, diffmap = p.Results.diffdata; end
    create_topo_plot(diffmap, chanlocs, figtitle, filepath, p.Results.cmap_diff, p.Results.clim_diff, plot_settings, true, p.Results.difftitle);
end

function create_topo_plot(data, chanlocs, title_name, filepath, cmap_input, clim_input, plot_settings, show_title, display_title)
    % display_title (optional): plain title shown on the figure, decoupled from
    % title_name (which sets the filename). Used so the effect map can read e.g.
    % "Correlation r-map" instead of the "A minus B" difference wording.
    if nargin < 9, display_title = ''; end
    fig = figure('Position', [50 50 450 450], 'Color', 'w', 'PaperPositionMode', 'auto');
    
    topoplot(data, chanlocs, plot_settings{:});
    
    % Handle different colormap inputs
    if ischar(cmap_input)
        colormap(cmap_input);
    elseif iscell(cmap_input)
        colormap(cmap_input{1});
    elseif isnumeric(cmap_input)
        colormap(cmap_input);
    end
    
    % Handle color limits
    if strcmp(clim_input, 'auto')
        caxis('auto');
    else
        caxis(clim_input);
    end
    
    axis off;
    set(gca, 'XLim', [-0.55 0.55], 'YLim', [-0.59 0.59]);
    delete(findobj(gca, 'Marker', '.'));
    
    if show_title && ~isempty(display_title)
        title(display_title, 'Interpreter', 'none');
    elseif show_title
        % Handle multi-line titles with VS separator
        if contains(title_name, ' VS ')
            title_parts = strsplit(title_name, ' VS ');
            if length(title_parts) == 2
                % Create two-line title
                title({title_parts{1}; ['minus ' title_parts{2}]}, 'Interpreter', 'none');
            else
                % Fallback to single line if unexpected format
                title(title_name, 'Interpreter', 'none');
            end
        else
            title(title_name, 'Interpreter', 'none');
        end
    end
    
    colorbar;
    filename = fullfile(filepath, [title_name, '_topo.png']);
    print(fig, '-dpng', '-r300', filename);
    close(fig);
end