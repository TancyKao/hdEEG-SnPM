function [G] = global_stat_test(data_x, data_y, alpha,comparison,tail)

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
        [~,p,~,STATS] = ttest(avg_dataX,avg_dataY,alpha,'right');
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
        [~,p,~,STATS] = ttest2(avg_dataX,avg_dataY,alpha,'right');
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
        G.pval = p;
    case 'correlationSboth'
        [r_corr, p] = corr(avg_dataX,avg_dataY,'Type','Spearman');
        G.stat = r_corr;
        G.pval = p;
    otherwise
        error('Error - improper comparison: %s (comparison: %s, tail: %s)', compstring, comparison, tail);
end






   