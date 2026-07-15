function build_compu_montage(target_degree)
%BUILD_COMPU_MONTAGE  Create the Compumedics/Neuvo montage assets from the .sfp.
%
%   build_compu_montage           % target mean degree 6 (default)
%   build_compu_montage(degree)   % custom target mean neighbour degree
%
%   Reads electrode coordinates from Compumedics-257.sfp, keeps the analysis
%   scalp set (hdeeg_scalpchannels('compu257'), excludes cheek/neck + REF),
%   and writes two files used by snpm_montage_registry / snpm_setup_channels:
%     dependencies/compu257_chanlocs.mat        (var compu257_chanlocs)
%     dependencies/NeighborMatrix_compu257.mat  (var neighbors)
%
%   The neighbour matrix uses the same distance-threshold construction as
%   build_178_neighbors.m, so the cluster/TFCE machinery runs unchanged.
%
%   See also BUILD_178_NEIGHBORS, SNPM_MONTAGE_REGISTRY, SNPM_SETUP_CHANNELS.

    if nargin < 1 || isempty(target_degree), target_degree = 6; end
    here = fileparts(mfilename('fullpath'));

    % --- read all 257 sfp electrodes ---
    sfp = fullfile(here, 'Compumedics-257.sfp');
    loc = readlocs(sfp, 'filetype', 'sfp');
    sfpLabels = {loc.labels};

    % --- keep only the analysis scalp channels, in the canonical list order ---
    analyze = hdeeg_scalpchannels('compu257');
    [tf, idx] = ismember(analyze, sfpLabels);
    if ~all(tf)
        error('build_compu_montage:missingLabels', ...
            'compu257 labels not in the .sfp: %s', strjoin(analyze(~tf), ', '));
    end
    compu257_chanlocs = loc(idx);
    nCh = numel(idx);

    % --- distance-threshold adjacency from X/Y/Z (cf. build_178_neighbors) ---
    cl  = loc(idx);
    XYZ = [[cl.X]', [cl.Y]', [cl.Z]'];
    assert(size(XYZ,1) == nCh && ~any(isnan(XYZ(:))), 'Compumedics chanlocs missing X/Y/Z');

    Dmat = squareform(pdist(XYZ));
    Dmat(1:nCh+1:end) = Inf;                    % ignore self
    nn = min(Dmat, [], 2);
    base = median(nn);
    factors = 1.0:0.05:4.0;
    chosen = factors(end);
    for f = factors
        A = Dmat <= f * base;
        if mean(sum(A, 2)) >= target_degree, chosen = f; break; end
    end
    A = Dmat <= chosen * base;

    if ~is_connected(A)
        warning('build_compu_montage:disconnected', ...
            'Adjacency graph is not fully connected at threshold factor %.2f.', chosen);
    end

    deg = sum(A, 2);
    neighbors = nan(nCh, max(deg));
    for c = 1:nCh
        nb = find(A(c, :));
        neighbors(c, 1:numel(nb)) = nb;
    end

    fprintf('Compumedics montage: %d scalp channels; neighbours factor=%.2f, degree min=%d mean=%.1f max=%d, connected=%d\n', ...
        nCh, chosen, min(deg), mean(deg), max(deg), is_connected(A));

    save(fullfile(here, 'compu257_chanlocs.mat'), 'compu257_chanlocs');
    save(fullfile(here, 'NeighborMatrix_compu257.mat'), 'neighbors');
    fprintf('Saved compu257_chanlocs.mat and NeighborMatrix_compu257.mat\n');
end

function tf = is_connected(A)
    n = size(A,1);
    seen = false(n,1); stack = 1; seen(1) = true;
    while ~isempty(stack)
        v = stack(end); stack(end) = [];
        nb = find(A(v,:) & ~seen');
        seen(nb) = true; stack = [stack, nb]; %#ok<AGROW>
    end
    tf = all(seen);
end
