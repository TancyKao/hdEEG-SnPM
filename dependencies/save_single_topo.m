function save_single_topo(data, chanlocs, title_str, out_png, cmap, clim)
% Single-map scalp topography saved to OUT_PNG. Mirrors create_topo_plot in
% dependencies/plot_topoInd.m (settings, sizing, dpi) so the per-group means,
% F-map, and group x condition cell means match the look of the legacy
% difference maps.
%
% CLIM (optional) sets explicit colour limits. Signed maps (e.g. post-hoc
% pairwise t) must pass symmetric limits [-m m] so 0 sits at the colormap
% centre; without it 'minmax'/'auto' scaling is asymmetric and a mixed +/-
% map can read as all-one-colour. Passing an explicit CLIM is also how the
% cell-mean grid keeps a shared scale across cells. Omit CLIM (or pass []) to
% keep the legacy data-driven 'minmax' scaling used by the 'hot' F-map /
% group means.
    if nargin < 6, clim = []; end
    if isempty(clim), maplim = 'minmax'; else, maplim = clim; end
    plot_settings = {'headrad', 0.6, 'style', 'map', 'electrodes', 'on', ...
        'maplimits', maplim, 'whitebk', 'on'};
    fig = figure('Position', [50 50 450 450], 'Color', 'w', 'PaperPositionMode', 'auto');
    cl = onCleanup(@() close(fig)); %#ok<NASGU>
    topoplot(data, chanlocs, plot_settings{:});
    if ischar(cmap) || iscell(cmap) || isnumeric(cmap)
        if iscell(cmap), colormap(cmap{1}); else, colormap(cmap); end
    end
    if isempty(clim), caxis('auto'); else, caxis(clim); end
    axis off;
    set(gca, 'XLim', [-0.55 0.55], 'YLim', [-0.59 0.59]);
    delete(findobj(gca, 'Marker', '.'));
    title(title_str, 'Interpreter', 'none');
    colorbar;
    print(fig, '-dpng', '-r300', out_png);
end
