function [A, info] = snpm_circ_import(Araw, opts)
%SNPM_CIRC_IMPORT  Units, zero-convention and angular-resolution gate for angles.
%
%   [A, info] = snpm_circ_import(Araw, opts)
%
% Takes a raw subjects x channels angle matrix straight off a detector export
% and returns radians on [0, 2*pi) with zero at the SLOW-OSCILLATION UP-STATE
% PEAK, plus a record of exactly what was done to it. Every circular analysis
% goes through here so the three live traps below are handled in one place.
%
% TRAP 1 - UNITS. There is NO default. Reading degrees as radians silently
% scrambles every angle and nothing downstream can detect it, so opts.units is
% REQUIRED and out-of-range values are a hard error, never a silent rescale.
%
% TRAP 2 - WHERE ZERO SITS. Detectors disagree about which point of the slow
% oscillation their angle calls zero. opts.convention is REQUIRED:
%   'literature_uppeak0'  0 = up-state peak (the standard)     rotation   0 deg
%   'yasa_uppeak0'        YASA, 0 = up-state peak              rotation   0 deg
%   'luna_zerocross0'     Luna COUPL_ANGLE, up-peak at 270 deg rotation  90 deg
%   'turtlewave_pre_v4'   TurtleWave < v4.0, exactly inverted  rotation 180 deg
%   'custom'              rotation = opts.zero_offset_deg
% The applied rotation is returned in info.rotation_deg and must be recorded in
% the results struct and the report header.
%
% TRAP 3 - ANGULAR RESOLUTION. Some exports are whole-degree
% (Luna's COUPL_ANGLE, TurtleWave's degree column). Resolution is measured on
% the finite angles in degrees as max|deg - round(deg)|:
%   < 1e-6            whole-degree export
%   < 0.05 deg        suspiciously quantised
% What that costs depends ENTIRELY on the statistic, so the gate is
% path-dependent and the caller chooses via opts.resolution_gate:
%   'error' - Watson U^2 only. U^2 is ECDF-based, so ties are structural:
%             whole-degree input moves the uncorrected p by up to 0.084 and the
%             topographic rank correlation to 0.958.
%   'warn'  - Hotelling and the Tier-1 transform (DEFAULT). Both are smooth
%             functions of the angle with no ranks anywhere; a <=0.5 deg
%             perturbation against a 15-50 deg between-subject circular SD is a
%             fraction of a percent. Erroring here would obstruct a real and
%             harmless case.
%
% NORMALISATION is applied uniformly with mod(., 2*pi) on entry, with no
% per-measure exception. The Tier-1 SIGNED measure references the pooled grand
% mean, which co-rotates with the data, so it is rotation-invariant (the same
% property that makes it permutation-exact); only the UNSIGNED measure depends
% on the declared convention being right, which is exactly what the inversion
% check below protects.
%
% INVERSION CHECK (run separately by the caller via the 'checkinversion' mode,
% because it needs the pooled sample across BOTH groups and the estimability
% mask). See SNPM_CIRC_CHECK_INVERSION below - it is exposed as a second output
% mode of this file's sibling helper in core_snpm_circ.
%
% INPUTS
%   Araw : subjects x channels raw angles (NaN = missing)
%   opts : struct
%     .units            REQUIRED 'rad'|'radians'|'deg'|'degrees'
%     .convention       REQUIRED, one of the five keys above
%     .zero_offset_deg  required only for 'custom'
%     .resolution_gate  'warn' (default) | 'error'
%     .label            free text used in messages (e.g. 'Group A angles')
%
% OUTPUTS
%   A    : subjects x channels radians on [0, 2*pi), zero = up-state peak
%   info : .units .convention .rotation_deg .resolution_deg .quantised
%          .whole_degree .n_finite .label
%
% See also SNPM_CIRC_HOTELLING, SNPM_CIRC_WATSONS_U2, SNPM_CIRC_LINEARISE.

    if nargin < 2, opts = struct(); end
    label = getopt(opts, 'label', 'angles');
    gate  = lower(getopt(opts, 'resolution_gate', 'warn'));
    if ~ismember(gate, {'warn', 'error'})
        error('core_snpm:circResolutionGate', ...
            'opts.resolution_gate must be ''warn'' or ''error''.');
    end

    % ---- units: required, never guessed -------------------------------
    units = getopt(opts, 'units', '');
    if isempty(units) || ~(ischar(units) || isstring(units))
        error('core_snpm:circUnitsRequired', ...
            ['params.circ_units is REQUIRED for a circular analysis and has no ' ...
             'default: it must be ''rad'' or ''deg''. Reading degrees as radians ' ...
             'scrambles every angle and no downstream check can recover it.']);
    end
    units = lower(char(units));
    switch units
        case {'rad', 'radians'}, units = 'rad';
        case {'deg', 'degrees'}, units = 'deg';
        otherwise
            error('core_snpm:circUnitsRequired', ...
                'params.circ_units must be ''rad'' or ''deg''; got ''%s''.', units);
    end

    fin  = isfinite(Araw);
    vals = Araw(fin);
    if isempty(vals)
        error('core_snpm:circNoAngles', '%s contains no finite angles.', label);
    end
    mx = max(abs(vals));

    if strcmp(units, 'rad')
        if mx > 2*pi + 1e-9
            error('core_snpm:circUnitsOutOfRange', ...
                ['%s was declared as RADIANS but contains |angle| = %.3f > 2*pi. ' ...
                 'This is almost certainly a degree export. Fix the declaration ' ...
                 '(params.circ_units = ''deg''); the units are never converted ' ...
                 'silently.'], label, mx);
        end
        Adeg = Araw * 180 / pi;
    else
        if mx > 360 + 1e-9
            error('core_snpm:circUnitsOutOfRange', ...
                ['%s was declared as DEGREES but contains |angle| = %.3f > 360. ' ...
                 'Check the export; the units are never converted silently.'], ...
                label, mx);
        end
        Adeg = Araw;
    end

    % ---- angular-resolution gate --------------------------------------
    dvals = Adeg(fin);
    res   = max(abs(dvals - round(dvals)));
    whole = res < 1e-6;
    quant = res < 0.05;
    if quant
        if whole
            what = sprintf('whole-degree (max|deg - round(deg)| = %.3g)', res);
        else
            what = sprintf('quantised to ~%.3g deg', res);
        end
        if strcmp(gate, 'error')
            error('core_snpm:circResolutionTooCoarse', ...
                ['%s is %s, which Watson''s U^2 cannot use. U^2 is built on the ' ...
                 'empirical distribution functions, so tied angles are structural, ' ...
                 'not cosmetic: whole-degree input moves the uncorrected p by up to ' ...
                 '0.084 and drops the topographic rank correlation against the ' ...
                 'full-precision map to 0.958. Three remedies, any one is enough: ' ...
                 '(1) export preferred_phase_rad -- TurtleWave stores it as a float, ' ...
                 'so this costs nothing; (2) re-derive the circular mean at full ' ...
                 'precision instead of reading Luna''s whole-degree COUPL_ANGLE; ' ...
                 '(3) use the Hotelling path (comparison = ''circ_phase_group''), ' ...
                 'which is the default and is unaffected because it is a smooth ' ...
                 'function of the angle with no ranks anywhere.'], label, what);
        end
        warning('core_snpm:circResolutionCoarse', ...
            ['%s is %s. The Hotelling / Tier-1 paths are smooth functions of the ' ...
             'angle, so a <=0.5 deg perturbation against a 15-50 deg between-subject ' ...
             'circular SD costs a fraction of a percent -- proceeding. Watson''s U^2 ' ...
             'would REJECT this input.'], label, what);
    end

    % ---- convention -> rotation ---------------------------------------
    conv = getopt(opts, 'convention', '');
    if isempty(conv) || ~(ischar(conv) || isstring(conv))
        error('core_snpm:circConventionRequired', ...
            ['params.circ_convention is REQUIRED and has no default. Choose one of: ' ...
             'literature_uppeak0 | luna_zerocross0 | yasa_uppeak0 | ' ...
             'turtlewave_pre_v4 | custom (with params.circ_zero_offset_deg).']);
    end
    conv = lower(char(conv));
    switch conv
        case {'literature_uppeak0', 'yasa_uppeak0'}
            rot = 0;
        case 'luna_zerocross0'
            rot = 90;      % Luna COUPL_ANGLE puts the up-state peak at 270 deg
        case 'turtlewave_pre_v4'
            rot = 180;     % TurtleWave < v4.0 is exactly inverted
        case 'custom'
            rot = getopt(opts, 'zero_offset_deg', []);
            if isempty(rot) || ~isnumeric(rot) || ~isscalar(rot) || ~isfinite(rot)
                error('core_snpm:circConventionRequired', ...
                    ['convention ''custom'' needs a finite scalar ' ...
                     'params.circ_zero_offset_deg (degrees to add so that 0 lands ' ...
                     'on the slow-oscillation up-state peak).']);
            end
        otherwise
            error('core_snpm:circConventionRequired', ...
                ['Unknown params.circ_convention ''%s''. Choose one of: ' ...
                 'literature_uppeak0 | luna_zerocross0 | yasa_uppeak0 | ' ...
                 'turtlewave_pre_v4 | custom.'], conv);
    end

    A = mod((Adeg + rot) * pi / 180, 2*pi);
    A(~fin) = NaN;

    info = struct('units', units, 'convention', conv, 'rotation_deg', mod(rot, 360), ...
        'resolution_deg', res, 'quantised', quant, 'whole_degree', whole, ...
        'n_finite', numel(vals), 'label', label);
end

% ======================================================================
function v = getopt(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
        v = s.(f);
    else
        v = d;
    end
end
