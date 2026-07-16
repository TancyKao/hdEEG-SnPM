function [G] = global_stat_test(data_x, data_y, alpha,comparison,tail,ncov)
% NCOV (optional, default 0): number of covariates residualized out of the data
% before this global test. For the partial-correlation cases the descriptive
% global p-value must use df = n - 2 - ncov, not the n - 2 that corr() reports.
% ncov = 0 leaves the ordinary (zero-order) correlation p unchanged.

if nargin < 6 || isempty(ncov), ncov = 0; end

avg_dataX = mean(data_x,2,'omitnan');
avg_dataY = mean(data_y,2,'omitnan');

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



   