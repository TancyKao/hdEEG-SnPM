function [G] = global_stat_test(data_x, data_y, alpha,comparison,tail,ncov)
%GLOBAL_STAT_TEST  Whole-head omnibus on the channel-averaged signal.
%
%   G = GLOBAL_STAT_TEST(data_x, data_y, alpha, comparison, tail, ncov)
%
%   G.stat / G.pval           the channel-averaged statistic and its p
%   G.n_channels              HOW MANY channels that average was taken over
%   G.n_channels_total        how many were offered
%   G.channels_kept           1 x nCh logical, the retained set itself
%
% NCOV (optional, default 0): number of covariates residualized out of the data
% before this global test. For the partial-correlation cases the descriptive
% global p-value must use df = n - 2 - ncov, not the n - 2 that corr() reports.
% ncov = 0 leaves the ordinary (zero-order) correlation p unchanged.
%
% COMMON-CHANNEL RULE (2026-08-03). The average is taken over the channels that
% are FINITE IN EVERY ANALYSED ROW OF BOTH data_x AND data_y, and over nothing
% else. It used to be mean(data,2,'omitnan'), i.e. every unit averaged over
% whatever channels it happened to have.
%
% WHY. Scalp power has a strong spatial gradient -- roughly 4x front-to-back in
% NREM delta -- so which channels a unit is averaged over sets a per-unit offset
% of the same order as the effects this toolbox looks for. Two units averaged
% over different channel sets are not measuring the same quantity, and the
% difference between their means is then partly a difference between their
% channel sets. Measured on an 8x8 lattice with a 40 -> 10 uV^2 gradient, NO true
% effect anywhere, 4000 replicates, nominal 0.05 (+-0.0034 = 3 SE):
%
%   unpaired, 4 of 10 group-1 subjects lose an occipital band   0.1118
%   unpaired, 2 of 10                                           0.0740
%   paired, SAME region missing in both conditions              0.0505  (cancels)
%   paired, DIFFERENT region missing in condition A vs B        1.0000  <-- always
%
% The paired 1.0000 is the one that matters: two nights with different bad
% channels is routine in sleep EEG, and the condition-A-minus-condition-B
% difference then contains a fixed frontal-minus-occipital offset that no
% amount of data can average away. Under the rule below every one of those
% scenarios returns to nominal (measured 0.0470-0.0532), because the offset is
% now identical in both arms and cancels in the contrast.
%
% The rule is deliberately the same one core_snpm_lmm.m:318 already applies to
% its global test (average over the analysed channel set), so the three tiers
% agree.
%
% NO-OP FOR THE CORRELATION PATH. core_snpm_analysis applies the complete-column
% mask BEFORE this call, so an excluded channel is already all-NaN in both
% matrices and a retained one is finite in every row. The intersection computed
% here is therefore exactly the retained set, and averaging over it is bitwise
% identical to the previous mean(...,'omitnan') -- adding an exact zero in the
% omitnan accumulation is exact. Pinned by test_global_common_channels (C5).
%
% FINITE, NOT MERELY NON-NaN: datatype 'logscale' turns a zero-power cell into
% log10(0) = -Inf, which passes ~isnan and would drag the whole average to -Inf.

if nargin < 6 || isempty(ncov), ncov = 0; end

% Ensure comparison and tail are strings without extra whitespace
comparison = strtrim(comparison);
tail = strtrim(tail);
compstring = [comparison tail];

% Debug output
%fprintf('DEBUG global_stat_test: comparison="%s", tail="%s", compstring="%s"\n', comparison, tail, compstring);

if strcmp(comparison,'pairedT') ==1 | strcmp(comparison,'correlationP') ==1 | strcmp(comparison,'correlationS') ==1
    if size(data_x,1) ~= size(data_y,1)
        error('Group sizes must match for paired/correlation tests');
    end
end

[avg_dataX, avg_dataY, keep] = common_channel_average(data_x, data_y, compstring);
G.n_channels       = sum(keep);
G.n_channels_total = numel(keep);
G.channels_kept    = keep;



switch compstring
    
    case 'pairedTleft'
        [~,p,~,STATS] = ttest(avg_dataX,avg_dataY,alpha,'left');
        G.stat = STATS.tstat;
        G.pval = p;
    case 'pairedTright'
        [~,p,~,STATS] = ttest(avg_dataX,avg_dataY,alpha,'right');
        G.stat = STATS.tstat;
        G.pval = p;
    case 'pairedTboth'
        [~,p,~,STATS] = ttest(avg_dataX,avg_dataY);
        G.stat = STATS.tstat;
        G.pval = p;
    case 'onesampleTboth'
        [~,p,~,STATS] = ttest(avg_dataX,avg_dataY,alpha);
        G.stat = STATS.tstat;
        G.pval = p;
    case 'unpairedTleft'
        [~,p,~,STATS] = ttest2(avg_dataX,avg_dataY,alpha,'left');
        G.stat = STATS.tstat;
        G.pval = p;
    case 'unpairedTright'
        [~,p,~,STATS] = ttest2(avg_dataX,avg_dataY,alpha,'right');
        G.stat = STATS.tstat;
        G.pval = p;
    case 'unpairedTboth'
        [~,p,~,STATS] = ttest2(avg_dataX,avg_dataY,alpha);
        G.stat = STATS.tstat;
        G.pval = p;
    case 'correlationPboth'
        [r_corr, p] = corr(avg_dataX,avg_dataY);
        G.stat = r_corr;
        G.pval = adjust_corr_p(r_corr, numel(avg_dataX), ncov, p);
    case 'correlationSboth'
        [r_corr, p] = corr(avg_dataX,avg_dataY,'Type','Spearman');
        G.stat = r_corr;
        G.pval = adjust_corr_p(r_corr, numel(avg_dataX), ncov, p);
    otherwise
        error('Error - improper comparison: %s (comparison: %s, tail: %s)', compstring, comparison, tail);
end

end

function [ax, ay, keep] = common_channel_average(data_x, data_y, compstring)
% Average each row over the channels that are usable in EVERY analysed row of
% BOTH matrices. See the COMMON-CHANNEL RULE note in the header.
%
% data_x and data_y need not have the same number of ROWS (unpairedT), but they
% must have the same number of COLUMNS -- they are the same montage.
%
% When no channel survives there is nothing to average, so the average is a
% column of NaN and the builtins below return a NaN statistic and a NaN p
% (verified: ttest, ttest2 and corr all propagate rather than error). Returning
% NaN is the point: silently averaging over per-row channel sets is exactly the
% behaviour this rule exists to remove, so an omnibus that cannot be computed
% must say so rather than report a number built from incomparable quantities.

    if size(data_x, 2) ~= size(data_y, 2)
        error('global_stat_test:channelCountMismatch', ...
            ['data_x has %d channels and data_y has %d. The global omnibus averages ' ...
             'both arms over one common channel set, so they must be the same montage.'], ...
            size(data_x, 2), size(data_y, 2));
    end

    keep = all(isfinite(data_x), 1) & all(isfinite(data_y), 1);

    if ~any(keep)
        warning('global_stat_test:noCommonChannel', ...
            ['Global omnibus (%s): NO channel is usable in every analysed row of both ' ...
             'arms, so there is no common set to average over and the whole-head statistic ' ...
             'is undefined (returned as NaN). Averaging each row over its own available ' ...
             'channels instead would make the arms measure different things -- with a ' ...
             'front-to-back power gradient that alone can reject on null data. The ' ...
             'channel-wise maps are unaffected; inspect the per-channel missingness.'], ...
            compstring);
        ax = NaN(size(data_x, 1), 1);
        ay = NaN(size(data_y, 1), 1);
        return
    end

    ax = mean(data_x(:, keep), 2);
    ay = mean(data_y(:, keep), 2);
end

function pv = adjust_corr_p(r, n, ncov, p_default)
% Two-sided correlation p with df = n - 2 - ncov (partial correlation loses one
% df per covariate residualized out). With ncov = 0 return corr()'s own p so the
% ordinary zero-order correlation path is unchanged. Spearman uses the same
% t-approximation corr() applies for its asymptotic p.
    if ncov <= 0
        pv = p_default;
        return;
    end
    df = n - 2 - ncov;
    if df <= 0
        pv = NaN;
        return;
    end
    t  = r .* sqrt(df ./ max(1 - r.^2, eps));
    pv = 2 * tcdf(-abs(t), df);
end



   