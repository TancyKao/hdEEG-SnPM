function snpm_assert_source(chanlocs, neighbors, channel_labels)
%SNPM_ASSERT_SOURCE  Correctness guards for the source-space (2447-voxel) path.
%
%   snpm_assert_source(chanlocs, neighbors, channel_labels)
%
%   Fails fast (clear error) when the source-space invariants are violated,
%   mirroring the channel-count guard the GLM/LMM paths use. Call AFTER
%   snpm_setup_channels has matched labels and subset chanlocs/neighbours.
%
%   Checks:
%     1. Count invariant   n_sources == size(neighbors,1) == numel(chanlocs)
%                          == 2447.
%     2. Node-order        the matched data columns are exactly
%                          src0001 .. src2447 IN ORDER, and the (reordered)
%                          chanlocs labels agree -- so neighbour-graph row i,
%                          data column i and coordinate i all refer to the same
%                          voxel. Enforcing canonical order removes any silent
%                          misalignment between the statistic and the reported
%                          coordinates.
%
%   The non-negative-magnitude guard lives in the caller (it needs the RAW,
%   pre-transform data and the datatype), not here.
%
%   See also SNPM_SETUP_CHANNELS, BUILD_SOURCE2447_COORDS.

N = 2447;
channel_labels = cellstr(channel_labels(:));

n = numel(channel_labels);
if n ~= N || size(neighbors, 1) ~= N || numel(chanlocs) ~= N
    error('snpm:source:count', ...
        ['Source-space analysis requires exactly %d voxels, but got %d data ' ...
         'columns, %d neighbour rows, %d coordinate entries. The source matrix ' ...
         'must be subjects x %d (columns src0001..src%04d).'], ...
        N, n, size(neighbors, 1), numel(chanlocs), N, N);
end

expected = arrayfun(@(k) sprintf('src%04d', k), 1:N, 'UniformOutput', false);
if ~isequal(channel_labels(:).', expected)
    bad = find(~strcmp(channel_labels(:).', expected), 1);
    error('snpm:source:order', ...
        ['Source data columns must be src0001..src%04d in node order (graph ' ...
         'row i == voxel i). First mismatch at position %d: got "%s", ' ...
         'expected "%s". Reorder the columns to canonical node order.'], ...
        N, bad, channel_labels{bad}, expected{bad});
end

cl_labels = {chanlocs.labels};
if ~isequal(cl_labels(:).', expected)
    error('snpm:source:coordOrder', ...
        ['Coordinate table node order does not match the data/graph order. ' ...
         'The source coordinate asset must list src0001..src%04d in graph ' ...
         'node order (rebuild with build_source2447_coords).'], N);
end
end
