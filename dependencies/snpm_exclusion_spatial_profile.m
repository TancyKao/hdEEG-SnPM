function profile = snpm_exclusion_spatial_profile(evaluable, chanlocs)
%SNPM_EXCLUSION_SPATIAL_PROFILE  WHERE the excluded channels fell, not just how many.
%
%   profile = SNPM_EXCLUSION_SPATIAL_PROFILE(evaluable, chanlocs)
%
%   evaluable   1 x nCh logical, TRUE for the channels that survive the
%               complete-column rule (see core_snpm_analysis).
%   chanlocs    EEGLAB chanlocs struct array for those same nCh channels, in the
%               SAME order (post snpm_setup_channels, so index i is channel i).
%
%   WHY THIS EXISTS. The exclusion COUNT is a power statement: it says how much
%   of the montage is gone. The exclusion PATTERN is a validity statement: it
%   says WHICH hypotheses can still be tested. Losing 12 channels scattered over
%   the head costs a little sensitivity everywhere; losing 12 contiguous frontal
%   channels leaves a frontal hypothesis with almost no coverage while the
%   headline "12 of 178 excluded" reads identically in both cases. A reader who
%   only sees the count cannot tell those apart, so the count alone is not an
%   honest report of what the analysis can still answer.
%
%   WHAT IT REPORTS. Two orthogonal splits of the montage, taken at the MEDIAN
%   of the montage's own coordinates so each half holds about half the channels
%   by construction (this makes the reference fraction a property of the
%   montage, not of an assumed head model):
%     anterior / posterior   split on the chanlocs X axis (EEGLAB: +X = nose)
%     left / right           split on the chanlocs Y axis (EEGLAB: +Y = left)
%   For each axis it counts how many EXCLUDED channels fell on each side and
%   scores that count against the null "the exclusions are scattered over the
%   montage independently of position", i.e. Binomial(n_excluded, p) with p the
%   montage fraction on that side. The p-value is the exact two-sided binomial
%   probability (sum of all outcomes no more likely than the observed one).
%
%   LOPSIDED -- the threshold, and why. profile.lopsided is TRUE when EITHER
%   axis has an exact two-sided binomial p < 0.05, uncorrected across the two
%   axes. Rationale:
%     * It is a DIAGNOSTIC WARNING, not an inferential test, and nothing
%       downstream is conditioned on it. A false alarm costs the reader one
%       glance at the exclusion list; a miss lets a regionally-gutted montage be
%       reported as if it still covered the head. The asymmetric cost is why the
%       two axes are not Bonferroni-corrected.
%     * The exact test is self-regulating at small n, which is the property a
%       fixed proportion cutoff ("more than 2/3 on one side") lacks: 2 of 3
%       excluded channels being anterior is meaningless and is not flagged
%       (p = 1), while 5 of 5 anterior (p = 0.0625) is still not flagged and
%       20 of 24 (p < 0.002) is. No separate minimum-count guard is needed.
%     * The split is deliberately coarse. Detecting true spatial CONTIGUITY
%       (e.g. a neighbour-graph clustering test on the excluded set) would be a
%       sharper instrument, but it needs its own permutation null and would be a
%       second inferential procedure inside a data-quality report. The two-axis
%       split catches the case that actually threatens a sleep-EEG hypothesis --
%       a whole region gone -- with arithmetic a reader can check by hand.
%
%   OUTPUT struct (the call site uses .available/.summary/.lopsided):
%     .available   logical. FALSE when the profile cannot be computed -- no
%                  usable coordinates, a chanlocs/evaluable length mismatch, or
%                  placeholder (non-anatomical) source coordinates. The caller
%                  must then print nothing rather than an unfounded claim; this
%                  function never errors on missing coordinates.
%     .summary     one-line printable string (also used as the warning text).
%     .lopsided    logical, see threshold above.
%     .reason      why .available is FALSE ('' when it is TRUE).
%     .n_excluded  .n_channels
%     .axes        1x2 struct array (anterior/posterior, left/right) with
%                  .label .side_pos .side_neg .n_pos .n_neg .montage_frac_pos .p
%     .alpha       the flagging threshold actually used (0.05).
%
%   See also SNPM_EXCLUDED_CHANNEL_INFO, CORE_SNPM_ANALYSIS.

ALPHA = 0.05;

evaluable = reshape(logical(evaluable), 1, []);
nCh       = numel(evaluable);
excluded  = ~evaluable;
n_excl    = sum(excluded);

profile = struct('available', false, 'summary', '', 'lopsided', false, ...
    'reason', '', 'n_excluded', n_excl, 'n_channels', nCh, 'alpha', ALPHA);
profile.axes = struct([]);   % filled per axis below; assigned separately
                             % because struct() would consume an empty struct

% ---- guards. Every one of these returns a profile the caller can print
% nothing from; none of them is an error, because a missing coordinate is a
% reason to stay quiet, not a reason to kill an otherwise valid analysis.
if isempty(chanlocs) || numel(chanlocs) ~= nCh
    profile.reason = sprintf(['coordinate list has %d entries for %d channels; ' ...
        'cannot place the exclusions'], numel(chanlocs), nCh);
    return
end

if is_placeholder_layout(chanlocs)
    % build_source2447_coords emits a deterministic graph-Laplacian eigenmap as
    % X/Y/Z until real MNI coordinates are supplied. Those axes are not
    % anterior-posterior or left-right in any anatomical sense, so an
    % "exclusions are frontal" claim read off them would be fabricated.
    profile.reason = ['source coordinates are the PLACEHOLDER graph embedding, ' ...
        'not anatomical axes; spatial profile suppressed'];
    return
end

ap = coord_vector(chanlocs, 'X');       % EEGLAB: +X toward the nose
lr = coord_vector(chanlocs, 'Y');       % EEGLAB: +Y toward the left ear
if ~(all(isfinite(ap)) && all(isfinite(lr)))
    % Polar fallback: many chanlocs assets carry theta/radius only.
    [ap, lr] = polar_fallback(chanlocs);
end
if ~(all(isfinite(ap)) && all(isfinite(lr)))
    profile.reason = 'chanlocs carry no usable X/Y (or theta/radius) coordinates';
    return
end

profile.available = true;

if n_excl == 0
    profile.summary = sprintf(['Exclusion pattern: none -- all %d channels are ' ...
        'complete across the analysed subjects.'], nCh);
    return
end

a1 = axis_split(ap, excluded, 'anterior/posterior', 'anterior', 'posterior');
a2 = axis_split(lr, excluded, 'left/right',         'left',     'right');
profile.axes = [a1, a2];

profile.lopsided = (a1.p < ALPHA) || (a2.p < ALPHA);

if profile.lopsided
    if a1.p <= a2.p, worst = a1; else, worst = a2; end
    if worst.n_pos >= worst.n_neg, side = worst.side_pos; else, side = worst.side_neg; end
    verdict = sprintf(['LOPSIDED: the excluded channels pile up %s (exact binomial ' ...
        'p = %s < %.2f), so this is a regional loss of coverage, not an even ' ...
        'thinning -- check that the channels a %s hypothesis needs are still in ' ...
        'the analysed set before reading the map there.'], ...
        side, fmt_p(worst.p), ALPHA, side);
else
    verdict = sprintf(['Consistent with scatter across the montage (neither axis ' ...
        'reaches p < %.2f).'], ALPHA);
end

profile.summary = sprintf(['Exclusion pattern: %d of %d channels excluded -- %s %d/%d ' ...
    '(montage %.0f%%/%.0f%%, p = %s), %s %d/%d (montage %.0f%%/%.0f%%, p = %s), ' ...
    'median-split on the montage coordinates. %s'], ...
    n_excl, nCh, ...
    a1.label, a1.n_pos, a1.n_neg, 100*a1.montage_frac_pos, 100*(1-a1.montage_frac_pos), fmt_p(a1.p), ...
    a2.label, a2.n_pos, a2.n_neg, 100*a2.montage_frac_pos, 100*(1-a2.montage_frac_pos), fmt_p(a2.p), ...
    verdict);
end

% ------------------------------------------------------------------ helpers

function a = axis_split(coord, excluded, label, side_pos, side_neg)
% Median split of the montage, then an exact two-sided binomial test of the
% excluded channels' side counts against the montage's own side fractions.
    med  = median(coord);
    pos  = coord > med;
    frac = mean(pos);
    k    = sum(pos & excluded);
    n    = sum(excluded);
    a = struct('label', label, 'side_pos', side_pos, 'side_neg', side_neg, ...
        'n_pos', k, 'n_neg', n - k, 'montage_frac_pos', frac, ...
        'p', binom_two_sided(k, n, frac));
end

function p = binom_two_sided(k, n, prob)
% Exact two-sided binomial p: the total probability of every outcome that is no
% more likely than the observed one. Computed from gammaln so it needs no
% toolbox function and is stable for the montage sizes here (n <= a few
% thousand). A degenerate reference fraction (all channels on one side) makes
% the test uninformative -> p = 1, never a flag.
    if n <= 0, p = 1; return; end
    if ~(prob > 0 && prob < 1), p = 1; return; end
    x   = (0:n)';
    lpm = gammaln(n+1) - gammaln(x+1) - gammaln(n-x+1) + x*log(prob) + (n-x)*log(1-prob);
    pmf = exp(lpm);
    obs = pmf(k+1);
    p   = sum(pmf(pmf <= obs * (1 + 1e-9)));
    p   = min(1, max(0, p));
end

function s = fmt_p(p)
    if p < 0.001, s = '<0.001'; else, s = sprintf('%.3f', p); end
end

function v = coord_vector(chanlocs, field)
% Per-channel scalar coordinate, NaN where absent/empty/non-finite. Written as a
% loop rather than [chanlocs.field] because an empty coordinate on ONE channel
% silently shortens the concatenation and would misalign every channel after it.
    v = nan(1, numel(chanlocs));
    if ~isfield(chanlocs, field), return; end
    for i = 1:numel(chanlocs)
        c = chanlocs(i).(field);
        if isnumeric(c) && isscalar(c) && isfinite(c)
            v(i) = double(c);
        end
    end
end

function [ap, lr] = polar_fallback(chanlocs)
% EEGLAB polar layout: theta in degrees (0 = nose, positive counter-clockwise
% toward the left), radius in [0, 0.5]. Projecting to Cartesian recovers the
% same two axes the X/Y branch uses; only the ordering along each axis matters
% for a median split, so the absolute scale is irrelevant.
    th = coord_vector(chanlocs, 'theta');
    rd = coord_vector(chanlocs, 'radius');
    ap = rd .* cosd(th);
    lr = rd .* sind(th);
end

function tf = is_placeholder_layout(chanlocs)
% TRUE for the 2447-voxel source asset, whose X/Y/Z are a reproducible graph
% embedding rather than anatomy (see build_source2447_coords: provenance
% .placeholder = true). Detected from the canonical src#### labels, which
% snpm_assert_source already enforces for that path.
    tf = false;
    if ~isfield(chanlocs, 'labels'), return; end
    labs = {chanlocs.labels};
    ok = cellfun(@(s) ischar(s) || isstring(s), labs);
    if ~all(ok) || isempty(labs), return; end
    tf = all(~cellfun(@isempty, regexp(cellstr(string(labs)), '^src\d+$', 'once')));
end
