function [Clusters, real_stat] = snpm_cluster_analysis(data_x,data_y,threshold,neighbors,alpha,comparison,tail,permutation_overide)
%% Outputs
% 1: Clusters   struct array of clusters with .channels .threshold .p .permutations
% 2: real_stat  the OBSERVED statistic map (1 x nChannels) this engine formed
%               its clusters from. Exposed so a test can assert that the cluster
%               engine and the TFCE engine compute the SAME observed statistic
%               for a given comparison. That assertion is not cosmetic: the
%               correlationS arm spent its whole life scoring a Spearman
%               observed map against a Pearson permutation null, and the visible
%               signature was exactly a disagreement between these two engines.
%               Optional -- single-output callers are unaffected.
%% Inputs
% 1: data_x
% 2: data_y
% 3: threshold 
% 4: neighbors (channel x neighbor matrix)
% 5: alpha
% 6: comparison ('unpairedT','pairedT'
% 7: tail ('both','right','left') in this case a left tail looks for
% negative threshold clusters
% 8: permutation_overide

% set data for channels that are not interesting to NaN (i.e. data for outside of 185 channels is set to NaN prior to inputing into matrix
%data_x(:,setdiff(1:256,inside185ch))=NaN;
%data_y(:,setdiff(1:256,inside185ch))=NaN;
% if you don't do this, clusters outside will be considered

%neighbors file is a channels by neighbors matrix (NaN if no neighbors, #
%columns is maximum number of neighbors
sparse_channel_adjacency_matrix = make_neighbors_sparse(neighbors,size(neighbors,1));

real_stat = [];   % second output; set from REALSTATS once the switch has run

if nargin < 5
    alpha = 0.05;
    disp('using default alpha value of 0.05');
end

if nargin < 6
    comparison = 'unpairedT';
    disp('using default unpaired T test');
end

if nargin < 7
    tail = 'both';
    disp('using default both tails');
    
end

if strcmp(comparison,'pairedT') ==1
    if size(data_x,1) ~= size(data_y,1)
        disp('group sizes must match')
        return;
    end
end

% use positive and negative thresholds if checking both
if strcmp(tail,'both') && length(threshold) == 1
    thresholds = [-threshold threshold];
elseif length(threshold) > 1
    disp ('continuing but >1 threshold is not technically valid')
else 
    thresholds = threshold;
end

compstring=[comparison tail];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this switch runs the permutation/combinations depending on the comparison
switch compstring

    case 'pairedTleft'

        nSubj=size(data_x,1);
        requested_permutations = 2^nSubj;

        if nargin == 8
            requested_permutations = permutation_overide;
        end

        % switch values so it does a left sided ttest
        temp = data_x;
        data_x = data_y;
        data_y = temp;
        clear temp;

        % Nichols & Holmes sign-flip labellings (exact if small, else random)
        signflip = snpm_signflip_patterns(nSubj, requested_permutations);
        possible_permutations = size(signflip,1);

        max_cluster_sizes = zeros(possible_permutations,length(thresholds));

        for permIndex=1:possible_permutations
            %print out the status
            if mod(permIndex, 1000) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end

            %% Calculate T value for this grouping (actual combos
            %creates groups from actual_combos which switches order of pairing or
            %keeps the same while maintaining pairing

            actual_combos=signflip(permIndex,:);

            data_x_temp=data_x(logical(actual_combos),:);
            data_x_temp=cat(1,data_x_temp,data_y(logical(~actual_combos),:));
            
            data_y_temp=data_y(logical(actual_combos),:);
            data_y_temp=cat(1,data_y_temp,data_x(logical(~actual_combos),:));
            
            %Calculate T value for this grouping (ttest is paired ttest)
            [~,~,~,STATS] = ttest(data_x_temp,data_y_temp,alpha,'right');
            % find clusters above threshold
            for ti = 1:length(thresholds)
                temp_clusters = snpm_find_clusters_graphalgs(STATS.tstat,thresholds(ti),sparse_channel_adjacency_matrix);
                if ~isempty(temp_clusters)
                    max_cluster_sizes(permIndex,ti) = max(cellfun('length',temp_clusters));
                end
            end
        end
        
        % calculate real ttest
        [~,~,~,REALSTATS] = ttest(data_x,data_y,alpha,'right');
    case 'pairedTright'

        nSubj=size(data_x,1);
        requested_permutations = 2^nSubj;

        if nargin == 8
            requested_permutations = permutation_overide;
        end

        % Nichols & Holmes sign-flip labellings (exact if small, else random)
        signflip = snpm_signflip_patterns(nSubj, requested_permutations);
        possible_permutations = size(signflip,1);

        max_cluster_sizes = zeros(possible_permutations,length(thresholds));

        for permIndex=1:possible_permutations
            %print out the status
            if mod(permIndex, 1000) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end


            actual_combos=signflip(permIndex,:);

            data_x_temp=data_x(logical(actual_combos),:);
            data_x_temp=cat(1,data_x_temp,data_y(logical(~actual_combos),:));
            
            data_y_temp=data_y(logical(actual_combos),:);
            data_y_temp=cat(1,data_y_temp,data_x(logical(~actual_combos),:));
            
            %Calculate T value for this grouping (ttest is paired ttest)
            [~,~,~,STATS] = ttest(data_x_temp,data_y_temp,alpha,'right');
            
            % find clusters above threshold
            for ti = 1:length(thresholds)
                temp_clusters = snpm_find_clusters_graphalgs(STATS.tstat,thresholds(ti),sparse_channel_adjacency_matrix);
                if ~isempty(temp_clusters)
                    max_cluster_sizes(permIndex,ti) = max(cellfun('length',temp_clusters));
                end
            end
        end
        
        %Calculate real ttest
        [~,~,~,REALSTATS] = ttest(data_x,data_y,alpha,'right');
    case 'pairedTboth'

        nSubj=size(data_x,1);
        requested_permutations = 2^nSubj;

        if nargin == 8
            requested_permutations = permutation_overide;
        end

        % Nichols & Holmes sign-flip labellings (exact if small, else random)
        signflip = snpm_signflip_patterns(nSubj, requested_permutations);
        possible_permutations = size(signflip,1);

        max_cluster_sizes = zeros(possible_permutations,length(thresholds));

        for permIndex=1:possible_permutations
            %print out the status
            if mod(permIndex, 1000) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end


            actual_combos=signflip(permIndex,:);

            data_x_temp=data_x(logical(actual_combos),:);
            data_x_temp=cat(1,data_x_temp,data_y(logical(~actual_combos),:));
            
            data_y_temp=data_y(logical(actual_combos),:);
            data_y_temp=cat(1,data_y_temp,data_x(logical(~actual_combos),:));
            
            %Calculate T value for this grouping (ttest is paired ttest)
            [~,~,~,REALSTATS] = ttest(data_x_temp,data_y_temp,alpha,'both');
            
            % find clusters above threshold
            for ti = 1:length(thresholds)
                temp_clusters = snpm_find_clusters_graphalgs(REALSTATS.tstat,thresholds(ti),sparse_channel_adjacency_matrix);
                if ~isempty(temp_clusters)
                    max_cluster_sizes(permIndex,ti) = max(cellfun('length',temp_clusters));
                end
            end
        end
        
        % find real ttest  values
        [~,~,~,REALSTATS] = ttest(data_x,data_y);    
        
    case 'onesampleTboth'

        nSubj=size(data_x,1);
        requested_permutations = 2^nSubj;

        if nargin == 8
            requested_permutations = permutation_overide;
        end

        % Nichols & Holmes sign-flip labellings (exact if small, else random)
        signflip = snpm_signflip_patterns(nSubj, requested_permutations);
        possible_permutations = size(signflip,1);

        max_cluster_sizes = zeros(possible_permutations,length(thresholds));

        for permIndex=1:possible_permutations
            %print out the status
            if mod(permIndex, 1000) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end


            actual_combos=signflip(permIndex,:);

            data_x_temp=data_x(logical(actual_combos),:);
            data_x_temp=cat(1,data_x_temp,data_y(logical(~actual_combos),:));
            
            data_y_temp=data_y(logical(actual_combos),:);
            data_y_temp=cat(1,data_y_temp,data_x(logical(~actual_combos),:));
            
            %Calculate T value for this grouping (ttest is paired ttest)
            [~,~,~,REALSTATS] = ttest(data_x_temp,data_y_temp,alpha,'both');
            
            % find clusters above threshold
            for ti = 1:length(thresholds)
                temp_clusters = snpm_find_clusters_graphalgs(REALSTATS.tstat,thresholds(ti),sparse_channel_adjacency_matrix);
                if ~isempty(temp_clusters)
                    max_cluster_sizes(permIndex,ti) = max(cellfun('length',temp_clusters));
                end
            end
        end
        
        % find real ttest  values
        [~,~,~,REALSTATS] = ttest(data_x,data_y); 
    case 'unpairedTleft'

        nSubj=size(data_x,1)+size(data_y,1);
        nGrp=size(data_y,1);
        G = nchoosek(nSubj,nGrp);

        data =  cat(1,data_y,data_x); clear data_*;

        % Honour an explicit override; otherwise enumerate exactly when the group
        % is small, else fall back to a Monte-Carlo ceiling. snpm_relabel_assignments
        % samples WITH REPLACEMENT (no uniqueness loop) so requested >= G never hangs.
        if nargin == 8
            requested_permutations = permutation_overide;
        elseif G <= 300000
            requested_permutations = G;          % exact enumeration
        else
            requested_permutations = 50000;      % default Monte-Carlo ceiling
        end

        actual_combos = snpm_relabel_assignments(nSubj,nGrp,requested_permutations);
        possible_permutations = size(actual_combos,1);

        max_cluster_sizes = zeros(possible_permutations,length(thresholds));

        for permIndex=1:possible_permutations
            %print out the status
            if mod(permIndex, 1000) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end

            group1 = actual_combos(permIndex,1:nGrp);
            group2 = actual_combos(permIndex,nGrp+1:end);

            %Calculate T value for this grouping (ttest2 is unpaired ttest)
            [~,~,~,REALSTATS] = ttest2(data(group1,:),data(group2,:),alpha,'right');

            % find clusters above threshold
            for ti = 1:length(thresholds)
                temp_clusters = snpm_find_clusters_graphalgs(REALSTATS.tstat,thresholds(ti),sparse_channel_adjacency_matrix);
                if ~isempty(temp_clusters)
                    max_cluster_sizes(permIndex,ti) = max(cellfun('length',temp_clusters));
                end
            end
        end
        % find real tvalue
        [~,~,~,REALSTATS] = ttest2(data(1:nGrp,:),data(nGrp+1:end,:),alpha,'right');
    case 'unpairedTright'

        nSubj=size(data_x,1)+size(data_y,1);
        nGrp=size(data_x,1);
        G = nchoosek(nSubj,nGrp);

        data =  cat(1,data_x,data_y); clear data_*;

        % See 'unpairedTleft'. snpm_relabel_assignments samples WITH REPLACEMENT
        % so requested >= G can no longer hang.
        if nargin == 8
            requested_permutations = permutation_overide;
        elseif G <= 300000
            requested_permutations = G;          % exact enumeration
        else
            requested_permutations = 50000;      % default Monte-Carlo ceiling
        end

        actual_combos = snpm_relabel_assignments(nSubj,nGrp,requested_permutations);
        possible_permutations = size(actual_combos,1);

        max_cluster_sizes = zeros(possible_permutations,length(thresholds));

        for permIndex=1:possible_permutations
            %print out the status
            if mod(permIndex, 1000) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end

            group1 = actual_combos(permIndex,1:nGrp);
            group2 = actual_combos(permIndex,nGrp+1:end);

            %Calculate T value for this grouping (ttest2 is unpaired ttest)
            [~,~,~,REALSTATS] = ttest2(data(group1,:),data(group2,:),alpha,'right');

            % find clusters above threshold
            for ti = 1:length(thresholds)
                temp_clusters = snpm_find_clusters_graphalgs(REALSTATS.tstat,thresholds(ti),sparse_channel_adjacency_matrix);
                if ~isempty(temp_clusters)
                    max_cluster_sizes(permIndex,ti) = max(cellfun('length',temp_clusters));
                end
            end
        end
        % find real tvalue
        [~,~,~,REALSTATS] = ttest2(data(1:nGrp,:),data(nGrp+1:end,:),alpha,'right');
    case 'unpairedTboth'

        nSubj=size(data_x,1)+size(data_y,1);
        nGrp=size(data_x,1);
        G = nchoosek(nSubj,nGrp);
        data =  cat(1,data_x,data_y); clear data_*;

        % See 'unpairedTleft'. snpm_relabel_assignments samples WITH REPLACEMENT
        % so requested >= G can no longer hang.
        if nargin == 8
            requested_permutations = permutation_overide;
        elseif G <= 300000
            requested_permutations = G;          % exact enumeration
        else
            requested_permutations = 50000;      % default Monte-Carlo ceiling
        end

        actual_combos = snpm_relabel_assignments(nSubj,nGrp,requested_permutations);
        possible_permutations = size(actual_combos,1);

        max_cluster_sizes = zeros(possible_permutations,length(thresholds));

        for permIndex=1:possible_permutations
            %print out the status
            if mod(permIndex, 1000) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end

            group1 = actual_combos(permIndex,1:nGrp);
            group2 = actual_combos(permIndex,nGrp+1:end);

            %Calculate T value for this grouping (ttest2 is unpaired ttest)
            [~,~,~,REALSTATS] = ttest2(data(group1,:),data(group2,:),alpha);

            % find clusters above threshold
            for ti = 1:length(thresholds)
                temp_clusters = snpm_find_clusters_graphalgs(REALSTATS.tstat,thresholds(ti),sparse_channel_adjacency_matrix);
                if ~isempty(temp_clusters)
                    max_cluster_sizes(permIndex,ti) = max(cellfun('length',temp_clusters));
                end
            end
        end

        % real t
        [~,~,~,REALSTATS] = ttest2(data(1:nGrp,:),data(nGrp+1:end,:));
    case 'correlationPboth'

        nSubj=size(data_x,1);
        % Honour an explicit override; otherwise use the sane default ceiling.
        % snpm_row_permutations enumerates exactly when requested >= nSubj!
        % (nSubj<=10) and otherwise samples WITH REPLACEMENT (no uniqueness loop),
        % so requested >= nSubj! cannot hang.
        if nargin == 8
            requested_permutations = permutation_overide;
        else
            requested_permutations = min(factorial(nSubj), 50000);
        end

        actual_combos = snpm_row_permutations(nSubj, requested_permutations);
        possible_permutations = size(actual_combos,1);

        max_cluster_sizes = zeros(possible_permutations,length(thresholds));

        % COMPLETE-COLUMN RULE (see snpm_corr_columns). A channel with any
        % missing cell is NaN in the observed map AND in every permuted map, so
        % the analysed row set cannot depend on the permutation. The old code
        % recomputed a per-channel `idxfinite` INSIDE this loop, which made the
        % number of complete pairs a function of the permutation.
        %
        % isfinite, not ~isnan: log10(0) = -Inf under datatype 'logscale' passes
        % ~isnan but still makes snpm_corr_columns return NaN, so an ~isnan mask
        % here would disagree with the statistic AND with the mask
        % core_snpm_analysis reports to the user. Same test in the Spearman arm
        % below, in snpm_single_threshold_with_TFCE, and in core_snpm_analysis.
        evaluable_col = all(isfinite(data_x) & isfinite(data_y), 1);

        for permIndex=1:possible_permutations
            %print out the status
            if mod(permIndex, 1000) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end

            randompermutations = actual_combos(permIndex,:);

            r_corr = snpm_corr_columns(data_x(randompermutations,:), data_y);
            r_corr(~evaluable_col) = NaN;   % defensive: unreachable
            REALSTATS.tstat=r_corr;  %this is an obvious misnomer until more appropriate terminology is established

            for ti = 1:length(thresholds)
                temp_clusters = snpm_find_clusters_graphalgs(REALSTATS.tstat,thresholds(ti),sparse_channel_adjacency_matrix);
                if ~isempty(temp_clusters)
                    max_cluster_sizes(permIndex,ti) = max(cellfun('length',temp_clusters));
                end
            end
        end

        % OBSERVED MAP VIA THE IDENTICAL CODE PATH AS THE PERMUTED MAPS.
        % Same snpm_corr_columns call as inside the loop, so the identity
        % relabelling (row 1 of actual_combos) reproduces it bitwise. Do not
        % replace this with a per-channel corr() loop: a one-ulp difference
        % between the observed statistic and its own image under the identity
        % relabelling flips a >= against the max-statistic null, changing b by
        % one and moving p by 1/(N+1) -- 0.005 at 200 permutations.
        REALSTATS.tstat = snpm_corr_columns(data_x, data_y);
        REALSTATS.tstat(~evaluable_col) = NaN;

    case 'correlationSboth'

        nSubj=size(data_x,1);
        % See 'correlationPboth'. snpm_row_permutations enumerates exactly when
        % requested >= nSubj! (nSubj<=10) and otherwise samples WITH REPLACEMENT
        % (no uniqueness loop), so requested >= nSubj! cannot hang.
        if nargin == 8
            requested_permutations = permutation_overide;
        else
            requested_permutations = min(factorial(nSubj), 50000);
        end

        actual_combos = snpm_row_permutations(nSubj, requested_permutations);
        possible_permutations = size(actual_combos,1);

        max_cluster_sizes = zeros(possible_permutations,length(thresholds));

        % COMPLETE-COLUMN RULE (see snpm_corr_columns) -- as in correlationPboth.
        % isfinite, not ~isnan, and read off the RAW matrices rather than the
        % ranks: tiedrank maps -Inf to a finite rank, so a mask taken after
        % ranking would silently readmit a channel the transform had destroyed.
        evaluable_col = all(isfinite(data_x) & isfinite(data_y), 1);

        % SPEARMAN NULL, NOT A PEARSON ONE.
        % This loop used to correlate the RAW data (corrcoef / snpm_corr_columns
        % on data_x, data_y) while the observed map below was Spearman, so the
        % observed value and its null were realisations of two different
        % functionals and no exchangeability argument connected them. Both
        % permutation nulls happen to have variance 1/(n-1) whatever the
        % marginals, so the error is NOT a scale error and does not show up as
        % an obviously wrong spread -- it is entirely in the tail SHAPE. At
        % n = 30 / 4000 permutations / 400 replicates on a single channel the
        % Pearson-null quantile ratios q95_P/q95_S and q99_P/q99_S were
        % 0.998/0.992 (gaussian) and 0.999/0.998 (uniform), i.e. invisible for
        % well-behaved data, but 0.899/1.537 with a single extreme point: the
        % Pearson permutation distribution goes BIMODAL there, because most
        % relabellings park the outlier somewhere unremarkable while a few align
        % it with an extreme partner. Single-channel rejection at nominal 0.05
        % was 0.0827 -- anti-conservative by 65% -- and the selection effect runs
        % the wrong way, because users pick correlationS precisely when they have
        % outliers or non-normal marginals.
        %
        % Rank ONCE, outside the loop. This is exact, not an approximation:
        % rank(x(perm)) == rank(x)(perm), so permuting ranks is the same thing as
        % ranking permuted data. tiedrank gives midranks, which is what
        % corr(...,'Type','Spearman') uses. The TFCE engine
        % (snpm_single_threshold_with_TFCE) already did it this way; this arm is
        % now identical to it.
        rank_x = tiedrank(data_x);
        rank_y = tiedrank(data_y);

        for permIndex=1:possible_permutations
            %print out the status
            if mod(permIndex, 1000) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end

            randompermutations = actual_combos(permIndex,:);

            r_corr = snpm_corr_columns(rank_x(randompermutations,:), rank_y);
            % NOT redundant in the Spearman arm: tiedrank maps -Inf to a finite
            % rank (unlike NaN, which it propagates), so a column made non-finite
            % by log10(0) would otherwise produce a real-looking r.
            r_corr(~evaluable_col) = NaN;
            REALSTATS.tstat=r_corr;  %this is an obvious misnomer until more appropriate terminology is established

            for ti = 1:length(thresholds)
                temp_clusters = snpm_find_clusters_graphalgs(REALSTATS.tstat,thresholds(ti),sparse_channel_adjacency_matrix);
                if ~isempty(temp_clusters)
                    max_cluster_sizes(permIndex,ti) = max(cellfun('length',temp_clusters));
                end
            end
        end

        % OBSERVED MAP VIA THE IDENTICAL CODE PATH AS THE PERMUTED MAPS.
        % Same rank_x/rank_y, same snpm_corr_columns call -- not a separate
        % "mathematically equivalent" routine. Do not "simplify" this into a
        % corr(...,'Type','Spearman') loop. With a max-statistic null and the
        % (b+1)/(N+1) estimator a one-ulp difference between the observed
        % statistic and its own image under the identity relabelling can flip a
        % >= comparison, change b by one, and move p by 1/(N+1) -- 0.005 at 200
        % permutations, enough to cross 0.05.
        REALSTATS.tstat = snpm_corr_columns(rank_x, rank_y);
        REALSTATS.tstat(~evaluable_col) = NaN;

    otherwise
        % FATAL. This used to print 'ERROR - improproper comparison' and then
        % CARRY ON: REALSTATS is never set for an unsupported comparison, so
        % execution continued into the cluster-forming code below with whatever
        % REALSTATS happened to be in scope (nothing, or a stale value from a
        % previous case) -- a typo'd comparison could silently produce a
        % Clusters struct. A comparison this engine has no permutation scheme
        % for is not a warning condition.
        %
        % Nothing legitimate reaches here from core_snpm_analysis: it calls
        % global_stat_test with the same comparison/tail FIRST, and that function
        % errors on exactly the same set of unsupported compstrings, so any
        % combination that would land here has already aborted the run. This
        % branch protects direct callers (export_report, scripts, tests).
        error('snpm_cluster_analysis:unsupportedComparison', ...
            ['No permutation scheme for comparison ''%s'' with tail ''%s'' (compstring ''%s''). ' ...
             'Supported here: pairedT and unpairedT (both|left|right), onesampleT (both), ' ...
             'correlationP and correlationS (both). The circular analyses ' ...
             '(circ_phase_group, circ_phase_group_u2, circ_corrAngLinear), the GLM presets ' ...
             '(anova1, ancova, regression, rmanova, mixed2way) and mixedmodel are NOT ' ...
             'compstring cases -- they run in their own pipelines and never reach this engine. ' ...
             'circ_wheeler_watson_Test and circ_WatsonsU2Test are retired.'], ...
            comparison, tail, compstring);
end
clear ti temp_clusters group*

% Expose the observed map (second output). Every case above ends by leaving the
% observed statistic in REALSTATS.tstat, produced by the same call the loop used.
real_stat = REALSTATS.tstat;

%T.real_T = REALSTATS.tstat;

% find clusters above thresholds
cl=1;
Clusters=struct([]);
for ti = 1:length(thresholds)
    temp_clusters = snpm_find_clusters_graphalgs(REALSTATS.tstat,thresholds(ti),sparse_channel_adjacency_matrix);
    if ~isempty(temp_clusters)
        for tcl = 1:length(temp_clusters)
            Clusters(cl,1).channels = temp_clusters{tcl};
            Clusters(cl,1).threshold=thresholds(ti);
            cl = cl +1;
        end
    else
        Clusters(cl,1).channels = [];
        Clusters(cl,1).threshold=thresholds(ti);
        Clusters(cl,1).p=1;
        cl = cl +1;
    end
    
end

% determine p-value based on
% Minimum-bias FWE p (Nichols & Holmes 2001; Phipson & Smyth 2010):
%   p = (#{max-null cluster size >= observed cluster size} + 1)/(N + 1)
% applied uniformly to every permutation tier so p is never zero. An empty
% (no-channel) cluster has size 0, so every null is >= it and p = 1.
if strcmp(tail,'both') && length(threshold) == 1
    max_cluster_sizes =  sort(max(max_cluster_sizes,[],2),1,'descend');
    for cli = 1 : length(Clusters)
        potential_clusters=Clusters(cli).channels;
        Clusters(cli).p = (sum(max_cluster_sizes >= length(potential_clusters)) + 1)/(possible_permutations + 1);
        Clusters(cli).permutations=possible_permutations;
    end

else
    max_cluster_sizes =  sort(max_cluster_sizes,1,'descend');
    for cli = 1 : length(Clusters)
        potential_clusters=Clusters(cli).channels;
        ti = Clusters(cli).threshold==thresholds;
        Clusters(cli).p = (sum(max_cluster_sizes(:,ti) >= length(potential_clusters)) + 1)/(possible_permutations + 1);
        Clusters(cli).permutations=possible_permutations;
    end

end

