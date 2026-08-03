function h = TopoplotSignificant_single(topodata, sigch_unP, sigch_correctedTFCE, clusters, chanlocs, insidegoodch, comparison, filepath, titletext, select_mode)
    % Plot topographical maps with significance markers for TFCE and Cluster corrections

    % select_mode ('label'|'positional') chooses the significance-marker index space:
    %   'label'      -> markers use already-filtered sequential channel indices
    %   'positional' -> markers are full-montage indices intersected with insidegoodch
    % Every recording system is label-based now; default keeps legacy callers working.
    if nargin < 10 || isempty(select_mode), select_mode = 'positional'; end
    is_178_channels = strcmp(select_mode, 'label');
    % Configuration for different plot types
    plot_configs = get_plot_configurations(comparison, 'data1', 'data2');
    
    % Create both plots
    create_significance_plot(topodata, sigch_unP, sigch_correctedTFCE, chanlocs, insidegoodch, ...
                           filepath, titletext, plot_configs, 'TFCE',is_178_channels);
    
    h = create_significance_plot(topodata, sigch_unP, clusters, chanlocs, insidegoodch, ...
                               filepath, titletext, plot_configs, 'Cluster',is_178_channels);
end

% function is_178_channels = detect_channel_system(topodata, chanlocs)
%     % Detect if we're using 178-channel system
%     is_178_channels = (length(chanlocs) == 178) && (sum(~isnan(topodata)) <= 178);
% end


function configs = get_plot_configurations(comparison, dat1, dat2)
    % Define plot configurations for different comparison types
    configs = struct();
    
    switch comparison
        case {'pairedT', 'unpairedT'}
            configs.colormap = 'jet';
            configs.clim = [-3, 3];
            configs.colorlabel = ['warm: ', dat1, ' > ', dat2];

        case 'onesampleT'
            configs.colormap = 'jet';
            configs.clim = [-3, 3];
            configs.colorlabel = ['warm: ', dat1, ' > 0'];

        case 'circ_phase_group'
            % Hotelling T^2 -> F on the (cos,sin) embedding: non-negative and
            % omnibus, so a sequential map with an auto scale, not a diverging one.
            configs.colormap = 'hot';
            configs.clim = 'auto';
            configs.colorlabel = 'Hotelling F (mean resultant differs)';

        case 'circ_phase_group_u2'
            configs.colormap = 'hot';
            configs.clim = 'auto';
            configs.colorlabel = 'Watson U^2 (distributions differ)';

        case 'circ_corrAngLinear'
            configs.colormap = 'hot';
            configs.clim = 'auto';
            configs.colorlabel = 'circular-linear F';

        case 'correlationS'
            configs.colormap = 'jet';
            configs.clim = [-1, 1];
            configs.colorlabel = 'warm: positive corr';

        case 'mixedmodel'
            % LMM statistic is signed t (continuous effect) or F (factor,
            % positive only); auto clim adapts to either scale.
            configs.colormap = 'jet';
            configs.clim = 'auto';
            configs.colorlabel = 'LMM effect (t / F)';

        case {'anova1', 'rmanova', 'mixed2way'}
            % omnibus F statistic (non-negative)
            configs.colormap = 'hot';
            configs.clim = 'auto';
            configs.colorlabel = 'F statistic';

        case {'ancova', 'regression'}
            % signed t (regression / 2-group) or F (>2-group ANCOVA); auto clim
            configs.colormap = 'jet';
            configs.clim = 'auto';
            configs.colorlabel = 'GLM effect (t / F)';

        case 'posthocT'
            % post-hoc pairwise t-contrast: signed statistic, diverging jet,
            % symmetric colour scale so 0 sits at the colormap centre.
            % t = mean_a - mean_b for label "A vs B" (warm = A > B).
            configs.colormap = 'jet';
            configs.clim = 'sym';
            configs.colorlabel = 't statistic (warm: first > second)';

        otherwise
            configs.colormap = 'jet';
            configs.clim = 'auto';
            configs.colorlabel = 'warm: positive corr';
    end
end

function h = create_significance_plot(topodata, sigch_unP, sigch_corrected, chanlocs, ...
    insidegoodch, filepath, titletext, configs, correction_type, is_178_channels)
    % Create individual significance plot

    % Resolve a symmetric colour scale for signed statistics (clim = 'sym'):
    % [-m m] centres 0 on the diverging colormap so a mixed +/- map isn't
    % dominated by one tail (the asymmetry plain 'auto'/'minmax' would give).
    if ischar(configs.clim) && strcmpi(configs.clim, 'sym')
        m = max(abs(topodata(~isnan(topodata))));
        if isempty(m) || m == 0 || ~isfinite(m), m = 1; end
        configs.clim = [-m m];
    end

    % Create figure
    fig = figure('Position', [20 50 500 500], 'Color', 'w', 'PaperPositionMode', 'auto');
    
    % Add legend text
    add_legend_text();
    
    % Create topographical plot
    topoplot(topodata, chanlocs, 'headrad', 0.60, 'style', 'map', 'electrodes', 'on', 'maplimits', 'minmax', 'whitebk', 'on');
    
    % Apply styling
    apply_plot_styling(configs, titletext, correction_type);
    
    % Set axis limits and get electrode positions
    [electrodes] = setup_electrodes();
    
    % Add significance markers
    add_significance_markers(electrodes, sigch_unP, sigch_corrected, chanlocs, insidegoodch, topodata, is_178_channels);

    % Save figure
    save_figure(fig, filepath, titletext, correction_type);
    
    h = fig;
    close(fig);
end

function add_legend_text()
    % Add legend for significance markers
    text(0.2, 0.58, 'B: uncorrectedSig', 'FontSize', 9, 'HorizontalAlignment', 'left');
    text(0.2, 0.55, 'W: correctedSig', 'FontSize', 9, 'HorizontalAlignment', 'left');
end

function apply_plot_styling(configs, titletext, correction_type)
    % Apply colormap, color limits, and title
    colormap(configs.colormap);
    
    if strcmp(configs.clim, 'auto')
        caxis('auto');
    else
        caxis(configs.clim);
    end
    
    colorbar;
    

    title_name = [titletext, ' ', correction_type];
    % Handle multi-line titles with VS separator
    if contains(title_name, ' VS ')
        title_parts = strsplit(title_name, ' VS ');
        if length(title_parts) == 2
            % Create two-line title
            title({title_parts{1}; ['VS ' title_parts{2}]; configs.colorlabel}, 'Interpreter', 'none');
        else
            % Fallback to single line if unexpected format
            title({title_name, configs.colorlabel}, 'Interpreter', 'none');
        end
    else
        title({title_name, configs.colorlabel}, 'Interpreter', 'none');
    end

    
end

function electrodes = setup_electrodes()
    % Set axis limits and extract electrode positions
    set(gca, 'XLim', [-0.55 0.55], 'YLim', [-0.59 0.59]);
    
    electrodes.x = get(findobj(gca, 'Marker', '.'), 'XData');
    electrodes.y = get(findobj(gca, 'Marker', '.'), 'YData');
    electrodes.z = get(findobj(gca, 'Marker', '.'), 'ZData');
    
    delete(findobj(gca, 'Marker', '.'));
end

function add_significance_markers(electrodes, sigch_unP, sigch_corrected, chanlocs, insidegoodch, topodata,is_178_channels)
    % Add significance markers based on channel configuration
    
    if is_178_channels
        add_markers_178_channels(electrodes, sigch_unP, sigch_corrected, insidegoodch);
    else
        add_markers_164_channels(electrodes, sigch_unP, sigch_corrected, insidegoodch);
    end
end

function add_markers_178_channels(electrodes, sigch_unP, sigch_corrected, insidegoodch)
    % Add significance markers for 178-channel configuration
    % Now handles properly filtered data from core_snpm_analysis.m
    % The data has been pre-filtered to match channel locations, so indices align correctly
    
    max_electrode_idx = length(electrodes.x);
    
    % Plot uncorrected significant channels
    if ~isempty(sigch_unP) && any(sigch_unP > 0)
        valid_indices = sigch_unP(sigch_unP > 0 & sigch_unP <= max_electrode_idx);
        if ~isempty(valid_indices)
            plot_markers(electrodes, valid_indices, [0 0 0], [0.5 0.5 0.5]);
        end
    end
    
    % Plot corrected significant channels
    if ~isempty(sigch_corrected) && any(sigch_corrected > 0)
        valid_indices = sigch_corrected(sigch_corrected > 0 & sigch_corrected <= max_electrode_idx);
        if ~isempty(valid_indices)
            plot_markers(electrodes, valid_indices, [1 1 1], [1 1 1]);
        end
    end
end

function add_markers_164_channels(electrodes, sigch_unP, sigch_corrected, insidegoodch)
    % Add significance markers for 164-channel configuration (original code)
    if length(electrodes.x) ~= length(insidegoodch)
        warning('Error with plotting sig channels: electrode count mismatch');
        return;
    end
    
    % Plot uncorrected significant channels
    if ~isempty(sigch_unP) && sum(sigch_unP) > 0
        [~, chi] = intersect(insidegoodch, sigch_unP);
        if ~isempty(chi)
            plot_markers(electrodes, chi, [0 0 0], [0.5 0.5 0.5]);
        end
    end
    
    % Plot corrected significant channels
    if ~isempty(sigch_corrected) && sum(sigch_corrected) > 0
        [~, chi3] = intersect(insidegoodch, sigch_corrected);
        if ~isempty(chi3)
            plot_markers(electrodes, chi3, [1 1 1], [1 1 1]);
        end
    end
end

function plot_markers(electrodes, indices, face_color, edge_color)
    % Plot significance markers at specified electrode positions
    if ~isempty(indices) && any(indices > 0)
        valid_indices = indices(indices > 0 & indices <= length(electrodes.x));
        if ~isempty(valid_indices)
            hold on;
            scatter(electrodes.x(valid_indices), electrodes.y(valid_indices), electrodes.z(valid_indices), ...
                   'filled', 'SizeData', 60, 'CData', face_color, 'MarkerEdgeColor', edge_color, 'LineWidth', 0.5);
        end
    end
end

function save_figure(fig, filepath, titletext, correction_type)
    % Save figure with appropriate filename
    titlename = [titletext, ' ', correction_type];
    filename = fullfile(filepath, [titlename, '.png']);
    print(fig, '-dpng', '-r300', filename);
end