function neighbors = build_178_neighbors(outfile, target_degree)
% Build a channel-adjacency (neighbour) matrix for the 178-channel hd-EEG
% montage from the 3-D electrode coordinates in EEG178chanlocs.mat, so the
% SnPM cluster/TFCE machinery (make_neighbors_sparse + ClusterEnhancement)
% can run on 178-channel data. Only a 256-channel neighbour matrix shipped
% with the tool; this fills the gap (and lets the 178 path work correctly).
%
% Method: a distance-threshold adjacency (cf. Maris & Oostenveld 2007). The
% threshold is auto-tuned so the mean number of neighbours is ~target_degree,
% then the graph is checked for connectivity.
%
% OUTPUT format matches NeighborMatrix_256: an nCh x maxNeighbours matrix of
% neighbour channel indices, NaN-padded. Saved as variable `neighbors`.
%
% USAGE:  build_178_neighbors;                 % saves dependencies/NeighborMatrix_178.mat
%         build_178_neighbors(outfile, 6);     % custom path / target degree

    if nargin < 2 || isempty(target_degree), target_degree = 6; end
    if nargin < 1 || isempty(outfile)
        here = fileparts(mfilename('fullpath'));
        outfile = fullfile(here, 'NeighborMatrix_178.mat');
    end

    S = load('EEG178chanlocs');
    chanlocs = S.EEG178chanlocs;
    nCh = numel(chanlocs);
    XYZ = [[chanlocs.X]', [chanlocs.Y]', [chanlocs.Z]'];
    assert(size(XYZ,1) == nCh && ~any(isnan(XYZ(:))), 'EEG178chanlocs missing X/Y/Z');

    % pairwise Euclidean distances
    Dmat = squareform(pdist(XYZ));
    Dmat(1:nCh+1:end) = Inf;                 % ignore self

    % auto-tune the distance threshold to hit the target mean degree
    nn = min(Dmat, [], 2);                    % nearest-neighbour distance per channel
    base = median(nn);
    factors = 1.0:0.05:4.0;
    chosen = factors(end);
    for f = factors
        A = Dmat <= f * base;
        if mean(sum(A, 2)) >= target_degree
            chosen = f; break
        end
    end
    thr = chosen * base;
    A = Dmat <= thr;

    % connectivity check (BFS from node 1)
    if ~is_connected(A)
        warning('build_178_neighbors:disconnected', ...
            'Adjacency graph is not fully connected at threshold factor %.2f.', chosen);
    end

    deg = sum(A, 2);
    maxdeg = max(deg);
    neighbors = nan(nCh, maxdeg);
    for c = 1:nCh
        nb = find(A(c, :));
        neighbors(c, 1:numel(nb)) = nb;
    end

    fprintf('178 neighbours: threshold factor=%.2f, degree min=%d mean=%.1f max=%d, connected=%d\n', ...
        chosen, min(deg), mean(deg), max(deg), is_connected(A));

    save(outfile, 'neighbors');
    fprintf('Saved %s\n', outfile);
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
