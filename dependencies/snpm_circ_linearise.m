function [Y, info, Ttab] = snpm_circ_linearise(A, measure, opts)
% Tier-1 circular-to-linear transform: angles in, a plain subjects x channels
% matrix out, ready for any existing GLM preset (anova1 / ancova / regression
% / rmanova / mixed2way) or legacy t-test.
%
%   [Y, info]       = snpm_circ_linearise(A, measure)
%   [Y, info]       = snpm_circ_linearise(A, measure, opts)
%   [Y, info, Ttab] = snpm_circ_linearise(...)      % table with the stamp column
%
% Tier 1 buys the whole existing design vocabulary (covariates, >2 groups,
% repeated measures, interactions) at the price of an approximation: it treats
% a bounded circular quantity as an unbounded linear one. That approximation
% is only defensible when the data are tightly concentrated, which is exactly
% what the guards below enforce. When the guards fire the channel is dropped
% rather than silently linearised.
%
% MEASURES
%   'signed'   (alias 'circ_signed_deviation')
%       y_i = circ_dist(alpha_i, theta_ref) on (-pi, pi], with theta_ref the
%       POOLED grand mean direction over ALL subjects in BOTH groups combined,
%       computed per channel. Using the pooled mean rather than a per-group
%       mean is what keeps the transform permutation-exact: the pooled mean is
%       a function of the unordered multiset of angles at that channel, so it
%       is invariant under relabelling and the transform commutes with the
%       permutation. A per-group reference would be recomputed inside each
%       permutation and would bias the null.
%       Interpretation: signed phase lead/lag relative to the group as a whole.
%
%   'unsigned' (alias 'circ_unsigned_distance')
%       y_i = abs(circ_dist(alpha_i, theta0)) in [0, pi], with theta0 an
%       A-PRIORI reference (default 0 = the up-state peak), NOT estimated from
%       the data.
%       Interpretation: how far this subject's preferred phase sits from the
%       up-state peak, regardless of direction.
%       *** This measure REQUIRES an event-count covariate. *** It is a
%       folded, dispersion-like quantity: a subject whose angle is estimated
%       from few events has a noisier angle, and noise inflates |deviation|
%       upward (it cannot go below 0). Any group difference in event count
%       therefore produces a spurious group difference in the measure. The
%       provenance stamp records this so a downstream guard can hard-error
%       with identifier core_snpm:circUnsignedNeedsCountCovariate when the
%       measure is pointed at a preset with no covariate slot.
%
% PER-CHANNEL GUARDS (a channel that fails ANY of these is dropped whole)
%   1. WRAPAROUND: any subject with abs(circ_dist(alpha_i, theta_ref)) >
%      opts.max_dev_deg (default 150 deg). Beyond that the subject is nearly
%      antipodal to the pooled mean and the sign of its linearised value is an
%      artefact of which side of the cut it landed on.
%   2. CONCENTRATION: pooled resultant length R < opts.min_R (default 0.3).
%      Below that there is no meaningful mean direction to deviate from.
%   3. MISSINGNESS: ANY subject missing at that channel.
%   The wraparound and concentration guards are both evaluated against the
%   POOLED mean direction for BOTH measures — if the sample is that dispersed,
%   neither linearisation is trustworthy, whichever reference the measure uses.
%
% *** WHOLE-COLUMN NaN, NEVER PER-CELL — THIS IS NOT A STYLE CHOICE. ***
% A channel that fails any guard, including a single missing subject, is
% returned as an ALL-NaN column. Per-cell NaN would be actively harmful here:
% the Freedman-Lane permutation context (snpm_glm_permute) computes Z\Y over
% whole columns, so one NaN cell turns that column all-NaN in every PERMUTED
% realisation while the OBSERVED statistic (snpm_glm_stat handles partial
% missingness per column) stays finite. The channel then carries a real
% observed value but contributes nothing to the max-statistic null, i.e. it is
% scored against a null built from the other channels only, and its corrected
% p is not a corrected p at all. Dropping the whole column makes the channel
% NaN on both sides, which is the honest and correct outcome. This behaviour
% was validated against the Freedman-Lane path.
%
% INPUTS
%   A       : nSubj x nCh angles. RADIANS by default; set opts.units='deg'
%             for degrees (converted on entry, output is always radians unless
%             opts.output_units='deg').
%   measure : 'signed' | 'unsigned' (or the long aliases above).
%   opts    : struct, all optional
%     .theta0          a-priori reference for the unsigned measure, radians.
%                      DEFAULT 0 (up-state peak).
%     .max_dev_deg     wraparound guard, degrees. DEFAULT 150.
%     .min_R           concentration guard. DEFAULT 0.3.
%     .units           'rad' (default) | 'deg' — units of A.
%     .output_units    'rad' (default) | 'deg' — units of Y.
%     .convention      free text describing the recording convention, stamped
%                      verbatim. DEFAULT 'radians_ccw_from_zero'.
%     .count_file      name of the per-subject event-count file this measure
%                      needs as a covariate. Stamped verbatim; leave empty
%                      only if you know the downstream guard will not need it.
%     .channel_labels  1 x nCh cellstr for the table output. DEFAULT Ch001...
%     .subject_ids     nSubj x 1 ids for the table output. Optional.
%
% OUTPUTS
%   Y     : nSubj x nCh linear matrix, all-NaN in dropped channels
%   info  : .measure, .theta_ref (1 x nCh), .R (1 x nCh pooled resultant),
%           .valid (1 x nCh logical), .reason (1 x nCh cellstr),
%           .n_dropped, .stamp (struct), .stamp_string (char)
%   Ttab  : table with (optional Subject) + one column per channel + the
%           reserved column 'circ_provenance' repeated on every row
%
% ------------------------------------------------------------------
% PROVENANCE STAMP FORMAT (a downstream guard reads this; grammar is fixed)
%
% The stamp is a single-line char, CSV-round-trip safe: no commas, no quotes,
% no newlines. It is key=value pairs joined by the pipe character '|', keys
% and values restricted to [A-Za-z0-9_.:+-] (anything else is replaced by '_').
% Example:
%
%   circ_stamp_version=1|measure=circ_unsigned_distance|units=rad|
%   convention=radians_ccw_from_zero|rotation=theta0_0|theta_ref=apriori|
%   requires_count_covariate=1|count_file=spindle_counts.csv
%
% (shown wrapped; it is one line in the file). Decode with:
%
%   kv = strsplit(s, '|');
%   for i = 1:numel(kv)
%       t = strsplit(kv{i}, '=');
%       stamp.(t{1}) = t{2};
%   end
%
% The field the guard keys on is requires_count_covariate ('1' or '0'); when
% it is '1' and the chosen preset has no covariate slot carrying the file in
% count_file, error with identifier core_snpm:circUnsignedNeedsCountCovariate.
% The stamp also travels as info.stamp (a struct with the same fields), so an
% in-memory caller need not parse anything.
% ------------------------------------------------------------------
%
% See also SNPM_CIRC_HOTELLING, SNPM_GLM_PERMUTE, CIRC_DIST, CIRC_MEAN, CIRC_R.

    if nargin < 3, opts = struct(); end

    measure = canonical_measure(measure);

    theta0       = getopt(opts, 'theta0', 0);
    max_dev_deg  = getopt(opts, 'max_dev_deg', 150);
    min_R        = getopt(opts, 'min_R', 0.3);
    units        = lower(getopt(opts, 'units', 'rad'));
    out_units    = lower(getopt(opts, 'output_units', 'rad'));
    convention   = getopt(opts, 'convention', 'radians_ccw_from_zero');
    count_file   = getopt(opts, 'count_file', '');

    if ~ismember(units, {'rad', 'deg'})
        error('snpm_circ_linearise:units', 'opts.units must be ''rad'' or ''deg''.');
    end
    if ~ismember(out_units, {'rad', 'deg'})
        error('snpm_circ_linearise:units', ...
            'opts.output_units must be ''rad'' or ''deg''.');
    end

    if strcmp(units, 'deg')
        A = A * pi / 180;
    end
    [nSubj, nCh] = size(A);

    % ---- pooled reference direction and concentration, per channel ----
    % circ_mean / circ_r are NOT NaN-aware, so they are computed on the
    % complex resultant directly; channels with any missing subject are
    % dropped anyway (guard 3).
    present   = ~isnan(A);
    Az        = exp(1i * A);
    Az(~present) = 0;
    npres     = sum(present, 1);
    Rsum      = sum(Az, 1);
    theta_ref = angle(Rsum);                        % 1 x nCh, (-pi, pi]
    Rlen      = abs(Rsum) ./ max(npres, 1);         % 1 x nCh mean resultant
    Rlen(npres == 0) = NaN;
    theta_ref(npres == 0) = NaN;

    % ---- deviation from the pooled mean, used by the guards ----
    dev_pooled = angle(exp(1i * (A - theta_ref)));  % nSubj x nCh, (-pi, pi]
    max_abs_dev = nan(1, nCh);
    if nSubj > 0
        d = abs(dev_pooled);
        d(~present) = NaN;
        max_abs_dev = max(d, [], 1, 'omitnan');
        max_abs_dev(npres == 0) = NaN;
    end

    % ---- guards -------------------------------------------------------
    anyMissing = any(~present, 1);
    wrapFail   = max_abs_dev > (max_dev_deg * pi / 180);
    concFail   = ~(Rlen >= min_R);

    valid  = ~(anyMissing | wrapFail | concFail);
    reason = repmat({''}, 1, nCh);
    for ch = 1:nCh
        r = {};
        if anyMissing(ch)
            r{end+1} = sprintf('missing_subjects=%d', nSubj - npres(ch)); %#ok<AGROW>
        end
        if wrapFail(ch)
            r{end+1} = sprintf('wraparound_max_dev=%.1fdeg', ...
                max_abs_dev(ch) * 180 / pi); %#ok<AGROW>
        end
        if concFail(ch)
            r{end+1} = sprintf('low_concentration_R=%.3f', Rlen(ch)); %#ok<AGROW>
        end
        reason{ch} = strjoin(r, ';');
    end

    % ---- the measure --------------------------------------------------
    switch measure
        case 'circ_signed_deviation'
            Y = dev_pooled;
            rotation = 'pooled_grand_mean_per_channel';
            theta_ref_kind = 'pooled_grand_mean';
            needs_count = 0;
        case 'circ_unsigned_distance'
            Y = abs(angle(exp(1i * (A - theta0))));
            rotation = sprintf('theta0_%g', theta0);
            theta_ref_kind = 'apriori';
            needs_count = 1;
        otherwise
            error('snpm_circ_linearise:measure', 'Unhandled measure ''%s''.', measure);
    end

    % WHOLE-COLUMN NaN for every channel that failed any guard.
    Y(:, ~valid) = NaN;
    Y(~present)  = NaN;      % belt and braces; those columns are already out

    if strcmp(out_units, 'deg')
        Y = Y * 180 / pi;
    end

    % ---- provenance stamp ---------------------------------------------
    stamp = struct( ...
        'circ_stamp_version', '1', ...
        'measure',            measure, ...
        'units',              out_units, ...
        'convention',         sanitize(convention), ...
        'rotation',           sanitize(rotation), ...
        'theta_ref',          theta_ref_kind, ...
        'requires_count_covariate', num2str(needs_count), ...
        'count_file',         sanitize(count_file));

    stamp_string = encode_stamp(stamp);

    info = struct('measure', measure, 'theta_ref', theta_ref, 'R', Rlen, ...
        'max_abs_dev_deg', max_abs_dev * 180 / pi, 'valid', valid, ...
        'reason', {reason}, 'n_dropped', sum(~valid), ...
        'n_present', npres, 'stamp', stamp, 'stamp_string', stamp_string);

    if nargout >= 3
        labels = getopt(opts, 'channel_labels', {});
        if isempty(labels)
            labels = arrayfun(@(k) sprintf('Ch%03d', k), 1:nCh, 'UniformOutput', false);
        end
        if numel(labels) ~= nCh
            error('snpm_circ_linearise:labelCount', ...
                'opts.channel_labels has %d entries, expected %d.', numel(labels), nCh);
        end
        Ttab = array2table(Y, 'VariableNames', matlab.lang.makeValidName(labels));
        subj = getopt(opts, 'subject_ids', []);
        if ~isempty(subj)
            if numel(subj) ~= nSubj
                error('snpm_circ_linearise:subjectCount', ...
                    'opts.subject_ids has %d entries, expected %d.', numel(subj), nSubj);
            end
            if isnumeric(subj)
                subj = arrayfun(@(v) sprintf('sub%03d', v), subj(:), ...
                    'UniformOutput', false);
            end
            Ttab = addvars(Ttab, string(subj(:)), 'Before', 1, ...
                'NewVariableNames', 'Subject');
        end
        Ttab.circ_provenance = repmat(string(stamp_string), nSubj, 1);
    end
end

% ======================================================================
function m = canonical_measure(measure)
    if ~(ischar(measure) || isstring(measure))
        error('snpm_circ_linearise:measure', 'measure must be text.');
    end
    switch lower(char(measure))
        case {'signed', 'signed_deviation', 'circ_signed_deviation'}
            m = 'circ_signed_deviation';
        case {'unsigned', 'unsigned_distance', 'circ_unsigned_distance'}
            m = 'circ_unsigned_distance';
        otherwise
            error('snpm_circ_linearise:measure', ...
                ['measure must be ''signed'' (circ_signed_deviation) or ' ...
                 '''unsigned'' (circ_unsigned_distance); got ''%s''.'], char(measure));
    end
end

% ======================================================================
function s = encode_stamp(stamp)
% key=value pairs joined by '|'. Field order is fixed so the string is
% reproducible; a decoder must not rely on order.
    f = fieldnames(stamp);
    parts = cell(1, numel(f));
    for i = 1:numel(f)
        parts{i} = sprintf('%s=%s', f{i}, stamp.(f{i}));
    end
    s = strjoin(parts, '|');
end

% ======================================================================
function v = sanitize(v)
% Restrict to the CSV-safe alphabet the stamp grammar allows. Empty stays
% empty but is encoded as 'none' so a key=value pair is never truncated.
    if isstring(v), v = char(v); end
    if isempty(v), v = 'none'; return; end
    if ~ischar(v), v = num2str(v); end
    v = regexprep(v, '[^A-Za-z0-9_.:+\-]', '_');
end

% ======================================================================
function v = getopt(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
        v = s.(f);
    else
        v = d;
    end
end
