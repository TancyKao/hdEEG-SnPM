function out = snpm_montage_registry(key)
%SNPM_MONTAGE_REGISTRY  Single source of truth for hd-EEG recording systems.
%
%   reg = snpm_montage_registry()       returns the struct array of all systems
%                                       (used to build the GUI dropdown).
%   m   = snpm_montage_registry(key)    returns the one system matching KEY,
%                                       resolving by key or by legacy alias
%                                       (case-insensitive); errors otherwise.
%
%   Each system struct has fields:
%     .key            internal key ('egi' | 'compu' | 'source2447')
%     .display        plain-language name for the GUI
%     .labels         cellstr of channel labels to ANALYZE (the scalp set)
%     .chanloc_file   .mat file holding the chanlocs struct (loaded by name)
%     .chanloc_var    variable name inside that .mat
%     .neighbor_file  .mat file holding the neighbour matrix
%     .neighbor_var   variable name inside that .mat (always 'neighbors')
%     .legacy_aliases cellstr of older keys that resolve to this system
%     .is_source      true for a source-space (cortical voxel) system whose
%                     spatial units are NOT scalp electrodes: scalp topoplots
%                     are bypassed and a significant-voxel list is emitted
%                     instead. false for scalp electrode montages.
%
%   Channel selection is LABEL-BASED for every system: snpm_setup_channels.m
%   matches a data file's column names against the montage chanlocs labels.
%   The legacy GUI strings '164 channels' / '178 channels' both resolve to EGI.
%
%   See also SNPM_SETUP_CHANNELS, HDEEG_SCALPCHANNELS.

systems = local_systems();

if nargin < 1 || isempty(key)
    out = systems;
    return;
end

key = char(key);

% 1) exact key match
for i = 1:numel(systems)
    if strcmpi(key, systems(i).key)
        out = systems(i);
        return;
    end
end
% 2) legacy alias match
for i = 1:numel(systems)
    if any(strcmpi(key, systems(i).legacy_aliases))
        out = systems(i);
        return;
    end
end

error('snpm:montage:unknownKey', ...
    'Unknown recording system "%s". Known systems: %s.', ...
    key, strjoin({systems.key}, ', '));
end

% ------------------------------------------------------------------------
function S = local_systems()
% EGI 256 (HydroCel): the egi257 scalp set IS the existing EEG178chanlocs
% montage (verified label-for-label), so EGI reuses the shipped 178 assets.
S(1) = struct( ...
    'key',            'egi', ...
    'display',        'EGI 256 (HydroCel)', ...
    'labels',         {hdeeg_scalpchannels('egi257')}, ...
    'chanloc_file',   'EEG178chanlocs', ...
    'chanloc_var',    'EEG178chanlocs', ...
    'neighbor_file',  'NeighborMatrix_178', ...
    'neighbor_var',   'neighbors', ...
    'is_source',      false, ...
    'legacy_aliases', {{'164 channels', '178 channels', 'egi257', 'egi256'}});

% Compumedics 257 (Neuvo): chanlocs + neighbours are built from
% Compumedics-257.sfp by build_compu_montage.m (the compu257 scalp set).
S(2) = struct( ...
    'key',            'compu', ...
    'display',        'Compumedics 257 (Neuvo)', ...
    'labels',         {hdeeg_scalpchannels('compu257')}, ...
    'chanloc_file',   'compu257_chanlocs', ...
    'chanloc_var',    'compu257_chanlocs', ...
    'neighbor_file',  'NeighborMatrix_compu257', ...
    'neighbor_var',   'neighbors', ...
    'is_source',      false, ...
    'legacy_aliases', {{'compu257', 'neuvo', 'neuvo256', 'compumedics'}});

% Source space (cortical voxels): the 2447-node sLORETA/template-head source
% grid. Spatial units are cortical voxels, NOT scalp electrodes -- the inverse
% solution is computed UPSTREAM (GeoSource); the toolbox only ingests the
% subjects x 2447 reconstructed band-power matrix and runs the same
% permutation t / TFCE / cluster-mass engines on the source adjacency graph.
% chanloc_file holds src0001..src2447 labels + a PLACEHOLDER graph-derived
% coordinate layout (source2447_coords.mat, built by build_source2447_coords.m)
% -- swap in the real GeoSource MNI xyz table (same 2447 node order) when
% available. is_source=true routes rendering to a significant-voxel list
% instead of a scalp topoplot.
S(3) = struct( ...
    'key',            'source2447', ...
    'display',        'Source 2447 (cortical voxels)', ...
    'labels',         {hdeeg_scalpchannels('source2447')}, ...
    'chanloc_file',   'source2447_coords', ...
    'chanloc_var',    'source2447_coords', ...
    'neighbor_file',  'NeighborMatrix_Sources_2447_Full', ...
    'neighbor_var',   'neighbors', ...
    'is_source',      true, ...
    'legacy_aliases', {{'source', 'sources', 'source2447', 'loreta2447', 'voxel2447'}});
end
