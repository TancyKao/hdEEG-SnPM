function [chanlocs, neighbors] = build_routine10_montage(target_degree)
% Build the montage files for a 10-channel routine PSG EEG (Fz F3 F4 Cz C3 C4
% Pz Oz O1 O2): standard 10-20 positions from EEGLAB's Standard-10-20-Cap81
% and a neighbour matrix built the same way as build_178_neighbors (distance
% threshold tuned to a target mean degree, connectivity checked).
%
% OUTPUT (saved next to this file):
%   chanlocs_routine10.mat       variable `chanlocs`, in the label order above
%   NeighborMatrix_routine10.mat variable `neighbors`, nCh x maxNeighbours, NaN-padded
%
% USAGE:  build_routine10_montage;      % target mean degree 3
%         build_routine10_montage(4);

    if nargin < 1 || isempty(target_degree), target_degree = 3; end
    labels = {'Fz','F3','F4','Cz','C3','C4','Pz','Oz','O1','O2'};
    here = fileparts(mfilename('fullpath'));

    locs = readlocs(fullfile(here, '..', 'eeglab2022.1', 'sample_locs', 'Standard-10-20-Cap81.locs'));
    [ok, ix] = ismember(lower(labels), lower({locs.labels}));
    assert(all(ok), 'Standard-10-20-Cap81 lacks %s', strjoin(labels(~ok), ' '));
    chanlocs = locs(ix);
    for c = 1:numel(chanlocs), chanlocs(c).labels = labels{c}; chanlocs(c).type = 'EEG'; end

    nCh = numel(chanlocs);
    XYZ = [[chanlocs.X]', [chanlocs.Y]', [chanlocs.Z]'];
    Dmat = squareform(pdist(XYZ));
    Dmat(1:nCh+1:end) = Inf;
    base = median(min(Dmat, [], 2));
    chosen = 4.0;
    for f = 1.0:0.05:4.0
        if mean(sum(Dmat <= f * base, 2)) >= target_degree, chosen = f; break, end
    end
    A = Dmat <= chosen * base;
    deg = sum(A, 2);
    neighbors = nan(nCh, max(deg));
    for c = 1:nCh
        nb = find(A(c, :));
        neighbors(c, 1:numel(nb)) = nb;
        fprintf('  %-3s (%d): %s\n', labels{c}, numel(nb), strjoin(labels(nb), ' '));
    end
    fprintf('routine10 neighbours: threshold factor=%.2f, degree min=%d mean=%.1f max=%d, connected=%d\n', ...
        chosen, min(deg), mean(deg), max(deg), is_connected(A));

    save(fullfile(here, 'chanlocs_routine10.mat'), 'chanlocs');
    save(fullfile(here, 'NeighborMatrix_routine10.mat'), 'neighbors');
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
