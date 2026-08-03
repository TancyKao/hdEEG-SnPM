function [A, logC, usecov, collinear, rep] = snpm_circ_count_covariate( ...
        A, gi, Cnt, params, chanlocs)
%SNPM_CIRC_COUNT_COVARIATE  Load, validate and condition the per-channel
%   log(event count) covariate for a two-group circular phase analysis.
%
%   [A, logC, usecov, collinear, rep] = snpm_circ_count_covariate( ...
%           A, gi, Cnt, params, chanlocs)
%
% File loading, subject/channel alignment and montage subsetting belong to the
% orchestrator; this takes the already-aligned nSubj x nCh count matrix Cnt (in
% the same subject order as A and gi) and does the statistics.
%
% WHY THE COUNT IS NOT OPTIONAL. A subject's preferred phase is estimated from
% their detected events, so its precision -- and therefore the resultant length
% the Hotelling test reads -- depends on how many events there were. Without the
% count there is no way to separate a genuine phase difference from a
% detection-rate difference.
%
% WHY IT IS PER CHANNEL. Event count varies strongly across the montage and the
% confound acts channel-wise, so a per-subject scalar under-adjusts exactly at
% the posterior channels where the problem is worst.
% params.circ_count_covariate_mode = 'per_channel' (default) | 'per_subject'.
%
% WHAT THIS RETURNS AND WHY
%   A         : the angle matrix with every cell NaN'd that has an angle but no
%               usable count. Silently exempting those cells would defeat the
%               guard, so the angle goes with the count. Above 5% of angle cells
%               that is a hard error: the signature means the wrong file, not
%               sparse data.
%   logC      : nSubj x nCh log(count), NaN where unusable
%   usecov    : 1 x nCh, FALSE where sd(log count) < 0.05 (constant counts mean
%               there is no confound to adjust, and the error df rises by one
%               there) or where the covariate is not estimable
%   collinear : 1 x nCh, log(count) disjoint by group or |point-biserial| above
%               params.circ_collinearity_max (default 0.8). Those channels are
%               not estimable and the caller must NaN them. Above 25% of
%               channels this errors for the whole analysis, naming
%               Helfrich-style event-count subsampling as the only valid route.
%   rep       : diagnostic record, including the per-channel point-biserial
%               correlation, which is PRINTED BEFORE any circular result.
%
% CHANLOCS is used only to name offending channels in the validation errors.
%
% See also CORE_SNPM_CIRC, SNPM_CIRC_FL_CONTEXT, SNPM_CIRC_HOTELLING.

    nCh = size(A, 2);
    if ~isequal(size(Cnt), size(A))
        error('core_snpm:circCountsShape', ...
            'Count matrix is %dx%d but the angle matrix is %dx%d.', ...
            size(Cnt,1), size(Cnt,2), size(A,1), size(A,2));
    end
    bad = isfinite(Cnt) & (Cnt < 0 | abs(Cnt - round(Cnt)) > 1e-9);
    if any(bad(:))
        [r, c] = find(bad, 5);
        items = arrayfun(@(i) sprintf('row %d, %s = %g', r(i), ...
            chanlocs(c(i)).labels, Cnt(r(i), c(i))), 1:numel(r), 'UniformOutput', false);
        error('core_snpm:circCountsInvalid', ...
            ['Event counts must be non-negative integers. Offending entries ' ...
             '(first %d): %s.'], numel(r), strjoin(items, '; '));
    end

    % A cell with an angle but no usable count cannot be adjusted, and silently
    % exempting it would defeat the guard -> the ANGLE is dropped too.
    haveA   = isfinite(A);
    noCount = haveA & ~(isfinite(Cnt) & Cnt > 0);
    nA = sum(haveA(:)); nMiss = sum(noCount(:));
    frac = 0; if nA > 0, frac = nMiss / nA; end
    if frac > 0.05
        error('core_snpm:circCountsMissing', ...
            ['%d of %d angle cells (%.1f%%) have no usable event count. Above 5%% ' ...
             'that signature means the wrong file was supplied, not sparse data: ' ...
             'check that the count files correspond to the same detection run as ' ...
             'the angle files.'], nMiss, nA, 100*frac);
    end
    if nMiss > 0
        warning('core_snpm:circCountsMissing', ...
            ['%d of %d angle cells (%.2f%%) have no usable event count; those ' ...
             'angles are dropped because the count adjustment cannot be applied ' ...
             'to them.'], nMiss, nA, 100*frac);
    end
    A(noCount) = NaN;

    mode = lower(getopt(params, 'circ_count_covariate_mode', 'per_channel'));
    Cuse = Cnt; Cuse(~(isfinite(Cuse) & Cuse > 0)) = NaN;
    if strcmp(mode, 'per_subject')
        logC = repmat(log(mean(Cuse, 2, 'omitnan')), 1, nCh);
    elseif strcmp(mode, 'per_channel')
        logC = log(Cuse);
    else
        error('core_snpm:circCountCovariateMode', ...
            'params.circ_count_covariate_mode must be ''per_channel'' or ''per_subject''.');
    end

    % ---- mandatory diagnostic, printed BEFORE any circular result -----
    pb = point_biserial(logC, gi);
    fprintf('\n--- Event-count / group confound diagnostic (read this first) ---\n');
    fprintf('  point-biserial r(log count, group): median %.2f, max |r| %.2f\n', ...
        median(pb, 'omitnan'), max(abs(pb)));
    fprintf(['  |r| above ~0.5 means most of the group contrast is being absorbed ' ...
        'by the\n  covariate and the analysis is near-uninformative whatever p emerges.\n']);
    fprintf('  channels with |r| > 0.5: %d/%d\n\n', sum(abs(pb) > 0.5), nCh);

    % ---- degenerate / non-estimable covariate --------------------------
    sd = std(logC, 0, 1, 'omitnan');
    usecov = ~(sd < 0.05) & any(isfinite(logC), 1);   % constant counts: no confound
    nconst = sum(~usecov);

    rmax = getopt(params, 'circ_collinearity_max', 0.8);
    disjoint = false(1, nCh);
    for ch = 1:nCh
        a = logC(gi == 1, ch); b = logC(gi == 0, ch);
        a = a(isfinite(a)); b = b(isfinite(b));
        if isempty(a) || isempty(b), disjoint(ch) = true; continue; end
        disjoint(ch) = (min(a) > max(b)) || (min(b) > max(a));
    end
    collinear = usecov & (disjoint | abs(pb) > rmax);
    if mean(collinear) > 0.25
        error('core_snpm:circCountCollinear', ...
            ['log(event count) is collinear with group at %d of %d channels ' ...
             '(%.0f%%, threshold 25%%): |point-biserial| > %.2f or the two groups '  ...
             'occupy disjoint count ranges. The count adjustment is simply not ' ...
             'estimable there, and no covariate can separate a phase effect from ' ...
             'a detection-rate effect. The only valid route is Helfrich-style ' ...
             'event-count subsampling: draw the same number of events per subject ' ...
             'in both groups before computing the preferred phase, then re-run.'], ...
            sum(collinear), nCh, 100*mean(collinear), rmax);
    end
    usecov(collinear) = false;

    rep = struct('mode', mode, 'n_dropped_constant', nconst, ...
        'n_collinear', sum(collinear), 'pb_max', max(abs(pb)), ...
        'pb_median', median(pb, 'omitnan'), 'pb', pb, ...
        'n_cells_no_count', nMiss, 'frac_cells_no_count', frac, ...
        'collinearity_max', rmax);
    fprintf(['Count covariate (%s): used at %d channels; dropped at %d (constant ' ...
        'log count, sd<0.05); %d non-estimable (collinear with group).\n'], ...
        mode, sum(usecov), nconst, sum(collinear));
end

function r = point_biserial(X, gi)
% Per-column correlation between a continuous column and a 0/1 indicator.
    n = size(X, 2); r = nan(1, n);
    g = double(gi(:));
    for k = 1:n
        v = X(:, k); ok = isfinite(v);
        if sum(ok) < 3 || numel(unique(g(ok))) < 2, continue; end
        c = corrcoef(v(ok), g(ok));
        r(k) = c(1, 2);
    end
end

% ======================================================================
function v = getopt(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
