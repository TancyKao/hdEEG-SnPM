function build_source2447_coords(outfile)
%BUILD_SOURCE2447_COORDS  Build the PLACEHOLDER source-voxel coordinate asset.
%
%   build_source2447_coords()          writes source2447_coords.mat next to
%                                       this file (dependencies/).
%   build_source2447_coords(outfile)   writes to OUTFILE instead.
%
%   Produces a chanlocs-like struct array `source2447_coords` with 2447 entries:
%       .labels   'src0001' .. 'src2447'  (zero-padded, graph-node order)
%       .X .Y .Z  a DETERMINISTIC placeholder 3-D layout derived from the source
%                 adjacency graph (graph-Laplacian eigenmap: the 2nd/3rd/4th
%                 eigenvectors, sign-normalised). These coordinates are NOT used
%                 for statistics (the neighbour graph already encodes adjacency)
%                 and topoplot rendering is bypassed for this system, so a
%                 placeholder layout is acceptable; it merely gives the
%                 significant-voxel table distinguishable, reproducible values.
%
%   PROVENANCE (also stored in the .mat as `provenance`):
%       *** PLACEHOLDER coordinates -- NOT real MNI space. ***
%   Replace with the user's GeoSource export: a 2447-row table with columns
%   label, X, Y, Z in the SAME node order as
%   NeighborMatrix_Sources_2447_Full.mat (row i = voxel i). Drop-in format:
%       T = readtable('geosource_mni.csv');   % columns: label, X, Y, Z
%       for i=1:2447, source2447_coords(i).labels=...; .X/.Y/.Z=T.X(i)/...; end
%       save('source2447_coords.mat','source2447_coords','provenance');
%   Keep node order == graph node order; snpm_assert_source enforces
%   src0001..src2447 in sequence.
%
%   See also SNPM_MONTAGE_REGISTRY, SNPM_ASSERT_SOURCE, HDEEG_SCALPCHANNELS.

if nargin < 1 || isempty(outfile)
    here = fileparts(mfilename('fullpath'));
    outfile = fullfile(here, 'source2447_coords.mat');
end

nb_file = 'NeighborMatrix_Sources_2447_Full.mat';
nb = load(nb_file);
neighbors = nb.neighbors;              % 2447 x 26, 1-indexed, NaN-padded
n = size(neighbors, 1);
assert(n == 2447, 'Expected 2447 source nodes, got %d.', n);

% --- build symmetric adjacency from the neighbour matrix ---
I = []; J = [];
for i = 1:n
    row = neighbors(i, ~isnan(neighbors(i, :)));
    I = [I; repmat(i, numel(row), 1)];   %#ok<AGROW>
    J = [J; row(:)];                     %#ok<AGROW>
end
A = sparse(I, J, 1, n, n);
A = max(A, A.');                          % symmetrise (26-connectivity is mutual)
A = spones(A);                            % binary adjacency
A = A - diag(diag(A));                    % no self-loops

% --- graph-Laplacian eigenmap (deterministic placeholder layout) ---
d = full(sum(A, 2));
L = spdiags(d, 0, n, n) - A;             % combinatorial Laplacian
% Full symmetric eig is deterministic (no random restarts); 2447x2447 is cheap.
[V, D] = eig(full(L));
[~, ord] = sort(diag(D), 'ascend');
V = V(:, ord);
% Skip the trivial constant eigenvector (smallest eigenvalue ~0); use the next
% three as a 3-D embedding. Sign-normalise each so the largest-magnitude entry
% is positive -> fully reproducible across runs/machines.
emb = V(:, 2:4);
for c = 1:size(emb, 2)
    [~, k] = max(abs(emb(:, c)));
    if emb(k, c) < 0, emb(:, c) = -emb(:, c); end
end
% Scale to a unit-ish sphere for readable magnitudes (purely cosmetic).
emb = emb ./ max(abs(emb(:)));

source2447_coords = struct('labels', cell(1, n), 'X', [], 'Y', [], 'Z', []);
for i = 1:n
    source2447_coords(i).labels = sprintf('src%04d', i);
    source2447_coords(i).X = emb(i, 1);
    source2447_coords(i).Y = emb(i, 2);
    source2447_coords(i).Z = emb(i, 3);
end

provenance = struct( ...
    'placeholder', true, ...
    'space', 'PLACEHOLDER graph-Laplacian eigenmap -- NOT real MNI coordinates', ...
    'source_graph', nb_file, ...
    'node_order', 'row i == voxel i == src####; matches the neighbour graph', ...
    'replace_with', 'GeoSource MNI xyz export (2447 rows: label,X,Y,Z in node order)', ...
    'built_by', mfilename, ...
    'built_on', datestr(now, 'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>

save(outfile, 'source2447_coords', 'provenance');
fprintf('Wrote %d source coordinates (PLACEHOLDER) to %s\n', n, outfile);
end
