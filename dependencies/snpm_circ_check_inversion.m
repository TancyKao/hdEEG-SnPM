function snpm_circ_check_inversion(A, gate, chanlocs, imp)
%SNPM_CIRC_CHECK_INVERSION  Hard-stop on a 180-degree phase-convention error.
%
%   snpm_circ_check_inversion(A, gate, chanlocs, imp)
%
% Runs on the POOLED sample across both groups, which makes it label-invariant
% and therefore incapable of contaminating the inference, and is restricted to
% FRONTAL channels that pass the estimability gate (the anterior third by the
% chanlocs X coordinate, EEGLAB convention: +X toward the nose).
%
% A pooled grand mean preferred phase inside [110, 250] degrees -- within 70
% degrees of exactly antiphase -- is a hard error with identifier
% core_snpm:circPhaseConventionInverted. Coupling to the DOWN state at that
% magnitude is not a finding; it is the known TurtleWave pre-v4.0 inversion.
% Once the user declares turtlewave_pre_v4 the check re-runs after the 180
% degree rotation and must then pass. A deviation of more than 70 degrees that
% is still outside the inverted window is a WARNING, not an error.
%
% Real-world calibration: an export sitting at +158.7 degrees errors here and
% repairs to -21.3 degrees once declared.
%
% INPUTS  A angles (radians, already rotated), gate 1 x nCh logical estimability
%         mask, chanlocs, imp the struct returned by SNPM_CIRC_IMPORT.
%
% See also SNPM_CIRC_IMPORT, CORE_SNPM_CIRC.
% Pooled grand mean direction over frontal, estimable channels. Label-invariant
% (it pools both groups), so running it cannot contaminate the inference.
    ch = find(gate);
    if isempty(ch), return; end
    if isfield(chanlocs, 'X')
        x = arrayfun(@(c) firstnum(c.X), chanlocs(ch));
        if any(isfinite(x))
            thr = quantile(x(isfinite(x)), 2/3);
            sel = ch(x >= thr);
            if numel(sel) >= 3, ch = sel; end
        end
    end
    a = A(:, ch); a = a(isfinite(a));
    if isempty(a), return; end
    mu  = mod(rad2deg(angle(sum(exp(1i * a)))), 360);
    dev = mod(mu + 180, 360) - 180;           % wrap to (-180, 180]
    fprintf(['Convention check: pooled grand mean direction over %d frontal ' ...
        'estimable channels = %.1f deg (%+.1f deg from the up-state peak)\n'], ...
        numel(ch), mu, dev);
    if abs(dev) >= 110
        error('core_snpm:circPhaseConventionInverted', ...
            ['The pooled grand mean preferred phase over frontal channels is ' ...
             '%.1f deg (%+.1f deg from the slow-oscillation up-state peak), i.e. ' ...
             'inside the inverted window [110, 250] deg -- within 70 deg of ' ...
             'exactly antiphase. Coupling to the DOWN state at this magnitude is ' ...
             'not a finding, it is the known TurtleWave pre-v4.0 phase inversion. ' ...
             'The declared convention was ''%s'' (rotation %+.1f deg applied). ' ...
             'Repaired it would read %.1f deg. Two remedies: (1) re-export from ' ...
             'TurtleWave v4.0 or later, which fixes the sign at source; or ' ...
             '(2) if you are certain the export predates v4.0, declare it -- ' ...
             'params.circ_convention = ''turtlewave_pre_v4'' -- and this check ' ...
             'will re-run after the 180 deg rotation and must then pass.'], ...
            mu, dev, imp.convention, imp.rotation_deg, mod(mu + 180, 360));
    elseif abs(dev) > 70
        warning('core_snpm:circPhaseConventionSuspect', ...
            ['Pooled grand mean preferred phase is %+.1f deg from the up-state ' ...
             'peak (convention ''%s''). That is a large deviation; confirm the ' ...
             'zero convention of your export before reading the result.'], ...
            dev, imp.convention);
    end
end

% ======================================================================
function v = firstnum(x)
    if isempty(x), v = NaN; else, v = double(x(1)); end
end
