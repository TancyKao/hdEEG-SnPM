function S = snpm_setup_channels(montage_key, data_cell, channel_cols)
%SNPM_SETUP_CHANNELS  Label-based channel selection for a recording system.
%
%   S = snpm_setup_channels(montage_key, data_cell, channel_cols)
%
%   Selects the channels to analyze by matching a data file's column names
%   (CHANNEL_COLS) against the recording system's chanlocs labels, then
%   subsets/reorders every data matrix, the chanlocs, and the neighbour matrix
%   into one shared channel order. This is the single channel-setup path for
%   all pipelines (core_snpm_analysis / core_snpm_glm / core_snpm_lmm), so a
%   recording system is added by editing snpm_montage_registry.m only.
%
%   Inputs
%     montage_key  recording-system key or legacy alias (see snpm_montage_registry)
%     data_cell    cell array of data matrices (subjects/trials x channels) to
%                  filter identically: {power} for GLM/LMM, {data_x,data_y} for
%                  the t-test pipeline. All must have the same column count/order.
%     channel_cols cellstr of the data's channel column names (CSV header order)
%
%   Output struct S
%     .data            cell array, same shape as DATA_CELL, columns filtered+reordered
%     .chanlocs        chanlocs struct, reordered to match the kept channels
%     .neighbors       neighbour matrix remapped to the kept order
%                      (nCh x maxNb, NaN-padded, size(neighbors,1)==nCh)
%     .insidegoodch    1:nCh (data is already filtered, so indices are sequential)
%     .channel_mapping struct describing the selection (for reports/debug)
%     .select_mode     'label' (every system is label-based)
%
%   See also SNPM_MONTAGE_REGISTRY, MAKE_NEIGHBORS_SPARSE.

m = snpm_montage_registry(montage_key);

% --- load montage assets (chanlocs + neighbour matrix) ---
try
    cl = load(m.chanloc_file);
    nb = load(m.neighbor_file);
catch ME
    error('snpm:montage:assetLoad', ...
        'Cannot load montage assets for "%s" (%s / %s): %s', ...
        m.key, m.chanloc_file, m.neighbor_file, ME.message);
end
chanlocs  = cl.(m.chanloc_var);
neighbors = nb.(m.neighbor_var);

% --- match data columns against the montage channel labels ---
[is_matched, mapping] = ismember(channel_cols, {chanlocs.labels});
vcols = find(is_matched);
if isempty(vcols)
    error('snpm:montage:noMatch', ...
        'No data columns match the %s montage labels.', m.display);
end
sel = mapping(is_matched);            % chanloc indices, in data-column order

% --- subset/reorder chanlocs + neighbours to the kept channels ---
chanlocs  = chanlocs(sel);
neighbors = remap_neighbors(neighbors, sel);

% --- filter every data matrix identically ---
data_out = cell(size(data_cell));
for k = 1:numel(data_cell)
    data_out{k} = data_cell{k}(:, vcols);
end

nch = numel(vcols);
if size(neighbors, 1) ~= nch
    error('snpm:montage:sizeMismatch', ...
        'Selected %d channels but neighbour matrix has %d rows for %s.', ...
        nch, size(neighbors, 1), m.display);
end

S = struct();
S.data         = data_out;
S.chanlocs     = chanlocs;
S.neighbors    = neighbors;
S.insidegoodch = 1:nch;
S.select_mode  = 'label';
S.channel_mapping = struct( ...
    'original_data_columns',   vcols, ...
    'original_chanloc_indices', sel, ...
    'filtered_data_size',      nch, ...
    'channel_labels',          {channel_cols(vcols)});
end

% ------------------------------------------------------------------------
function nb_re = remap_neighbors(nb, sel)
% Subset a neighbour matrix to the kept channels `sel` (original indices, in
% the new column order) and renumber the neighbour entries to the new indices,
% dropping neighbours that were not kept. Identity when sel == 1:nCh.
% (Lifted from core_snpm_glm.m's remap_neighbors, the proven reference.)
    remap = nan(size(nb, 1), 1);
    remap(sel) = 1:numel(sel);
    nb_sel = nb(sel, :);
    nb_re = nan(size(nb_sel));
    for i = 1:numel(sel)
        row = nb_sel(i, ~isnan(nb_sel(i, :)));
        row = remap(row);
        row = row(~isnan(row))';
        nb_re(i, 1:numel(row)) = row;
    end
    nb_re = nb_re(:, any(~isnan(nb_re), 1));   % trim all-NaN columns
    if isempty(nb_re), nb_re = nan(numel(sel), 1); end
end
