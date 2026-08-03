function outcsv = write_source_voxel_table(results_struct, chanlocs, outname)
%WRITE_SOURCE_VOXEL_TABLE  Statistically-sufficient deliverable for source space.
%
%   outcsv = write_source_voxel_table(results_struct, chanlocs, outname)
%
%   For the source-space (2447-voxel) system there is no scalp layout to plot,
%   so instead of a topoplot we emit the full per-voxel statistic map with
%   coordinates and cluster membership: one row per voxel, flagged for
%   uncorrected / TFCE / cluster significance. Writes "<outname>_sigvoxels.csv"
%   AND a 'sigVoxels' sheet in "<outname>.xlsx" (created by func_genSnpmTable).
%
%   Columns: idx, label, X, Y, Z, stat, p, pTFCE, sig_uncorr, sig_TFCE,
%            in_sig_cluster, cluster_id, cluster_p.
%
%   Coordinates come from CHANLOCS (.X/.Y/.Z). Until the real GeoSource MNI
%   export is dropped in (see build_source2447_coords) these are PLACEHOLDER
%   graph-layout values -- the stats, flags and cluster membership are exact.
%
%   See also CORE_SNPM_ANALYSIS, CORE_SNPM_LMM, BUILD_SOURCE2447_COORDS.

    S = results_struct;
    if isfield(S, 'T') && isfield(S.T, 'real_T')
        stat = S.T.real_T(:);
    elseif isfield(S, 'real_T')
        stat = S.real_T(:);
    else
        stat = nan(numel(chanlocs), 1);
    end
    if isfield(S, 'p') && isfield(S.p, 'real')
        p_real = S.p.real(:);
    else
        p_real = nan(numel(chanlocs), 1);
    end
    if isfield(S, 'p') && isfield(S.p, 'correctedTFCE')
        p_tfce = S.p.correctedTFCE(:);
    else
        p_tfce = nan(numel(chanlocs), 1);
    end

    n = numel(chanlocs);
    idx    = (1:n).';
    label  = string({chanlocs.labels}).';
    X = getcoord(chanlocs, 'X');
    Y = getcoord(chanlocs, 'Y');
    Z = getcoord(chanlocs, 'Z');

    % Pad/trim statistic vectors to n voxels defensively.
    stat   = fitn(stat, n);
    p_real = fitn(p_real, n);
    p_tfce = fitn(p_tfce, n);

    sig_uncorr = false(n, 1);
    if isfield(S, 'uncorrsigch') && ~isempty(S.uncorrsigch)
        sig_uncorr(intersectidx(S.uncorrsigch, n)) = true;
    end
    sig_tfce = false(n, 1);
    if isfield(S, 'correctTFCEsigch') && ~isempty(S.correctTFCEsigch)
        sig_tfce(intersectidx(S.correctTFCEsigch, n)) = true;
    end

    % Cluster membership: assign id + p to member voxels; flag significant ones.
    cluster_id = zeros(n, 1);
    cluster_p  = nan(n, 1);
    in_sig_cluster = false(n, 1);
    if isfield(S, 'Clusters') && ~isempty(S.Clusters)
        for c = 1:numel(S.Clusters)
            ch = S.Clusters(c).channels;
            ch = intersectidx(ch, n);
            if isempty(ch), continue; end
            cluster_id(ch) = c;
            cp = S.Clusters(c).p;
            cluster_p(ch)  = cp;
            if cp <= 0.05, in_sig_cluster(ch) = true; end
        end
    end

    Tbl = table(idx, label, X, Y, Z, stat, p_real, p_tfce, ...
        double(sig_uncorr), double(sig_tfce), double(in_sig_cluster), ...
        cluster_id, cluster_p, ...
        'VariableNames', {'idx','label','X','Y','Z','stat','p','pTFCE', ...
        'sig_uncorr','sig_TFCE','in_sig_cluster','cluster_id','cluster_p'});

    outcsv = [outname, '_sigvoxels.csv'];
    writetable(Tbl, outcsv);

    % Also drop it into the xlsx alongside the stat sheets, if present.
    try
        writetable(Tbl, [outname, '.xlsx'], 'Sheet', 'sigVoxels');
    catch
        % xlsx may not exist yet / be locked; the CSV is the primary output.
    end

    nsig = sum(sig_tfce | in_sig_cluster | sig_uncorr);
    fprintf('Source voxel table: %d/%d voxels flagged significant -> %s\n', ...
        nsig, n, outcsv);
end

% ------------------------------------------------------------------------
function v = getcoord(chanlocs, f)
    n = numel(chanlocs);
    v = nan(n, 1);
    if isfield(chanlocs, f)
        for i = 1:n
            val = chanlocs(i).(f);
            if ~isempty(val) && isnumeric(val), v(i) = val(1); end
        end
    end
end

function v = fitn(v, n)
    v = v(:);
    if numel(v) >= n, v = v(1:n); else, v(end+1:n, 1) = NaN; end
end

function ix = intersectidx(ix, n)
    ix = ix(:);
    ix = ix(isfinite(ix) & ix >= 1 & ix <= n);
    ix = unique(round(ix));
end
