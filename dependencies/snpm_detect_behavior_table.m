function B = snpm_detect_behavior_table(varnames)
%SNPM_DETECT_BEHAVIOR_TABLE  Header-only classifier for a correlation Data 2.
%
%   B = SNPM_DETECT_BEHAVIOR_TABLE(VARNAMES) inspects a table's column NAMES
%   only (no data is read) and reports which of three shapes it has:
%
%     'montage'   a subjects x channels table whose column names are the
%                 recording system's channel labels -- the ordinary Data 2 for
%                 a channel-wise correlation. Matched to Data 1 BY LABEL.
%     'behavior'  Subject plus exactly ONE non-channel measure column (a
%                 reaction time, a memory score): one value per subject, which
%                 core_snpm_analysis broadcasts across Data 1's channels.
%     'unknown'   several data columns, but too few of them are recognised
%                 channel labels to call it a montage. This is NOT an error
%                 here: it is what a valid montage export looks like when its
%                 columns are renamed or re-cased ('e1', 'E1_density'), and it
%                 is also what a behaviour file with a stray extra column looks
%                 like. Only the CALLER can tell those apart, because only the
%                 caller knows how many channels Data 1 has -- see the
%                 resolution rule below.
%
%   Inputs
%     VARNAMES  cellstr / string array / char of the table's variable names
%               (table.Properties.VariableNames).
%
%   Output struct B
%     .kind             'montage' | 'behavior' | 'unknown'
%     .is_behavior      true iff kind is 'behavior' (kept for older callers)
%     .col_name         the behavioural column name ('' unless kind='behavior')
%     .channel_cols     columns recognised as channel labels
%     .nonchannel_cols  data columns that are not recognised channel labels
%     .data_cols        every non-ID column, in file order (channel or not)
%     .subject_cols     the ID column(s) removed from consideration
%     .channel_fraction numel(channel_cols) / numel(data_cols), 0 when empty
%
%   CLASSIFICATION
%     ID column   the FIRST column whose lowercased name contains 'subject',
%                 and only that one. This matches what core_snpm_analysis does
%                 (find(contains(...,'subject'), 1)) and it is why a file with
%                 'Subject' and 'SubjectiveSleepiness' is read as one ID column
%                 plus one behavioural measure -- treating every 'subject*'
%                 column as an ID silently deleted sleepiness/alertness
%                 measures, which are named that way throughout sleep research.
%     channel     ^E\d+$  OR  Cz (case-insensitive)  OR  a member of ANY
%                 recording system's label set in snpm_montage_registry() (so
%                 Compumedics names like Fp1 and source-space src#### count).
%     kind        no data columns          -> error snpm:behavior:noDataColumns
%                 1 data column, unrecognised -> 'behavior'
%                 1 data column, recognised   -> 'montage'
%                 >1 data column, MAJORITY recognised (fraction > 0.5)
%                                             -> 'montage'
%                 >1 data column, otherwise   -> 'unknown'
%
%   WHY A MAJORITY AND NOT "ANY". Deciding 'montage' on the presence of a
%   single recognised label mis-classified whole real exports: a lowercased EGI
%   file (e1..e256, cz) or a suffixed one (E1_density..., Cz) has exactly ONE
%   recognised column, and calling that a 1-channel montage made the engine
%   reject a 178-column file for "having 1 channel". A montage table is
%   overwhelmingly channels; anything else is 'unknown'.
%
%   HOW A CALLER RESOLVES 'unknown' (core_snpm_analysis does exactly this):
%     numel(B.data_cols) == Data 1's data-column count
%         -> treat as a montage table matched to Data 1 BY POSITION, warn that
%            the labels could not be checked, and carry on (this is the
%            behaviour every such file had before the classifier existed).
%     otherwise
%         -> a real input error: the file is neither the same montage as Data 1
%            nor a single behavioural measure. Report both accepted formats and
%            name a couple of the columns actually found.
%
%   RETIRED: snpm:behavior:multipleColumns. This function used to raise it for
%   any table with no recognised channel column and more than one data column.
%   That conflated the two cases above and turned a valid renamed-montage run
%   into a hard error, so the ambiguity is now RETURNED ('unknown') instead of
%   thrown. A caller still catching that identifier will simply never see it.
%
%   Dependency-light on purpose (registry + regexp only): the GUI calls it on
%   every file pick.
%
%   See also SNPM_MONTAGE_REGISTRY, SNPM_SETUP_CHANNELS, CORE_SNPM_ANALYSIS.

    if nargin < 1 || isempty(varnames)
        varnames = {};
    end
    if ischar(varnames)
        varnames = {varnames};
    elseif isstring(varnames)
        varnames = cellstr(varnames);
    end
    if ~iscell(varnames)
        error('snpm:behavior:badInput', ...
            'snpm_detect_behavior_table expects a cellstr of variable names, got %s.', ...
            class(varnames));
    end
    varnames = varnames(:)';
    varnames = cellfun(@char, varnames, 'UniformOutput', false);

    % --- ID column: the FIRST 'subject' match only (core's rule) ---
    is_subject = false(1, numel(varnames));
    subj_idx = find(cellfun(@(s) contains(lower(s), 'subject'), varnames), 1);
    if ~isempty(subj_idx)
        is_subject(subj_idx) = true;
    end

    % --- channel columns: same test the GUI uses ---
    reg = snpm_montage_registry();
    all_labels = [reg.labels];
    is_chan = cellfun(@(s) ~isempty(regexp(s, '^E\d+$', 'once')) || strcmpi(s, 'Cz'), varnames) ...
              | ismember(varnames, all_labels);
    is_chan = is_chan & ~is_subject;

    channel_cols    = varnames(is_chan);
    nonchannel_cols = varnames(~is_chan & ~is_subject);
    data_cols       = varnames(~is_subject);

    B = struct();
    B.kind             = 'unknown';
    B.is_behavior      = false;
    B.col_name         = '';
    B.channel_cols     = channel_cols;
    B.nonchannel_cols  = nonchannel_cols;
    B.data_cols        = data_cols;
    B.subject_cols     = varnames(is_subject);
    B.channel_fraction = 0;

    if isempty(data_cols)
        error('snpm:behavior:noDataColumns', ...
            ['The table has no data column at all (columns: {%s}). Supply either a ' ...
             'subjects x channels table whose column names match the recording system, ' ...
             'or a Subject column plus ONE behavioural measure column.'], ...
            strjoin(varnames, ', '));
    end

    B.channel_fraction = numel(channel_cols) / numel(data_cols);

    if isscalar(data_cols)
        if isempty(channel_cols)
            B.kind        = 'behavior';
            B.is_behavior = true;
            B.col_name    = data_cols{1};
        else
            B.kind = 'montage';   % a genuine one-channel table
        end
        return;
    end

    if B.channel_fraction > 0.5
        B.kind = 'montage';
    else
        B.kind = 'unknown';
    end
end
