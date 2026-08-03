function [T,p,E,H]= snpm_single_threshold_with_TFCE(data_x,data_y,neighbors,E,H,alpha,comparison,tail,permutation_overide)

%data_x, data_y = subjects x inputs

%%Optional inputs and defaults
%%alpha=0.05
%%comparison='unpairedT' ('pairedT', potentially add more)

%%Outputs
%%T struct chkT,realT,critT, tMax
%Threshold K means cluser thresholds Hi_Threshold,Med_Threshold,Low_Threshold

if nargin<4  % default E H
    E=0.5;
    H=2;
end

if nargin<6   % default alpha
    alpha = 0.05;
    display('using default alpha value of 0.05');
end

if nargin<7 % default comparison
    comparison = 'unpairedT';
    display('using default unpaired T test');
end

if nargin<8 % default tail
    tail = 'both';
    display('using default both tails');
end

sparse_channel_adjacency_matrix = make_neighbors_sparse(neighbors,size(neighbors,1));

number_of_inputs=size(data_x,2);
compstring=[comparison tail];
% this switch sets up the comparisons
switch compstring

    case 'pairedTleft'

        nSubj=size(data_x,1);
        requested_permutations = 2^nSubj;
        if nargin == 9
            requested_permutations = permutation_overide;
        end
        % switch values so it does a left sided test
        temp = data_x;
        data_x = data_y;
        data_y = temp;
        clear temp;

        % Nichols & Holmes sign-flip labellings (exact if small, else random)
        signflip = snpm_signflip_patterns(nSubj, requested_permutations);
        possible_permutations = size(signflip,1);

        TFCEdata=zeros(possible_permutations,number_of_inputs);
        tVals=zeros(possible_permutations,number_of_inputs);

        for permIndex = 1:possible_permutations
            %print out the status
            if mod(permIndex, 100) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end
            %Calculate T value for this grouping (actual combos
            %creates groups from actual_combos which switches order of pairing or
            %keeps the same while maintaining pairing


            actual_combos=signflip(permIndex,:);

            data_x_temp=data_x(logical(actual_combos),:);
            data_x_temp=cat(1,data_x_temp,data_y(logical(~actual_combos),:));
            
            data_y_temp=data_y(logical(actual_combos),:);
            data_y_temp=cat(1,data_y_temp,data_x(logical(~actual_combos),:));
            
            %Calculate T value for this grouping (ttest is paired ttest)
            [~,~,~,STATS] = ttest(data_x_temp,data_y_temp,alpha,'right');
            tVals(permIndex,:) = STATS.tstat;
            TFCEdata(permIndex,:) = ClusterEnhancement(STATS.tstat,sparse_channel_adjacency_matrix,E,H);
        end
        clear  p STATS group1 group
        
        %calculate real t
        [~,p.real,~,REALSTATS] = ttest(data_x,data_y,alpha,'right');
        T.tMax = sort(max(tVals,[],2),'descend');
        T.chk_T = mean(tVals);
        T.tMaxTFCE = sort(max(TFCEdata,[],2),'descend');
        
        T.real_T = REALSTATS.tstat;
        T.real_TFCE = ClusterEnhancement(REALSTATS.tstat,sparse_channel_adjacency_matrix,E,H);
            
    case 'pairedTright'

        nSubj=size(data_x,1);
        requested_permutations = 2^nSubj;
        if nargin == 9
            requested_permutations = permutation_overide;
        end

        % Nichols & Holmes sign-flip labellings (exact if small, else random)
        signflip = snpm_signflip_patterns(nSubj, requested_permutations);
        possible_permutations = size(signflip,1);

        TFCEdata=zeros(possible_permutations,number_of_inputs);
        tVals=zeros(possible_permutations,number_of_inputs);

        for permIndex = 1:possible_permutations
            %print out the status
            if mod(permIndex, 100) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end

            actual_combos=signflip(permIndex,:);

            data_x_temp=data_x(logical(actual_combos),:);
            data_x_temp=cat(1,data_x_temp,data_y(logical(~actual_combos),:));
            
            data_y_temp=data_y(logical(actual_combos),:);
            data_y_temp=cat(1,data_y_temp,data_x(logical(~actual_combos),:));
            
            %Calculate T value for this grouping (ttest is paired ttest)
            [~,~,~,STATS] = ttest(data_x_temp,data_y_temp,alpha,'right');
            tVals(permIndex,:) = STATS.tstat;
            TFCEdata(permIndex,:) = ClusterEnhancement(STATS.tstat,sparse_channel_adjacency_matrix,E,H);
        end
        clear  p STATS group1 group
        
        %calculate real t
        [~,p.real,~,REALSTATS] = ttest(data_x,data_y,alpha,'right');
        T.tMax = sort(max(tVals,[],2),'descend');
        T.chk_T = mean(tVals);
        T.tMaxTFCE = sort(max(TFCEdata,[],2),'descend');
        
        T.real_T = REALSTATS.tstat;
        T.real_TFCE = ClusterEnhancement(REALSTATS.tstat,sparse_channel_adjacency_matrix,E,H);
          
    case 'pairedTboth'

        nSubj=size(data_x,1);
        requested_permutations = 2^nSubj;
        if nargin == 9
            requested_permutations = permutation_overide;
        end

        % Nichols & Holmes sign-flip labellings (exact if small, else random)
        signflip = snpm_signflip_patterns(nSubj, requested_permutations);
        possible_permutations = size(signflip,1);

        TFCEdata=zeros(possible_permutations,number_of_inputs);
        tVals=zeros(possible_permutations,number_of_inputs);


        for permIndex = 1:possible_permutations
            %print out the status
            if mod(permIndex, 100) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end

            actual_combos=signflip(permIndex,:);

            data_x_temp=data_x(logical(actual_combos),:);
            data_x_temp=cat(1,data_x_temp,data_y(logical(~actual_combos),:));
            
            data_y_temp=data_y(logical(actual_combos),:);
            data_y_temp=cat(1,data_y_temp,data_x(logical(~actual_combos),:));
            
            %Calculate T value for this grouping (ttest is paired ttest)
            [~,~,~,STATS] = ttest(data_x_temp,data_y_temp);
            tVals(permIndex,:) = STATS.tstat;
            TFCEdata(permIndex,:) = ClusterEnhancement(STATS.tstat,sparse_channel_adjacency_matrix,E,H);
        end
        clear  p STATS group1 group
        
        % calculate real t
        [~,p.real,~,REALSTATS] = ttest(data_x,data_y);
        T.tMax = sort(max(abs(tVals),[],2),'descend');
        T.chk_T = mean(abs(tVals));
        T.tMaxTFCE = sort(max(abs(TFCEdata),[],2),'descend');
        
        T.real_T = REALSTATS.tstat;
        T.real_TFCE = ClusterEnhancement(REALSTATS.tstat,sparse_channel_adjacency_matrix,E,H);
    
    case 'onesampleTboth'
        nSubj=size(data_x,1);
        requested_permutations = 2^nSubj;
        if nargin == 9
            requested_permutations = permutation_overide;
        end

        % Nichols & Holmes sign-flip labellings (exact if small, else random)
        signflip = snpm_signflip_patterns(nSubj, requested_permutations);
        possible_permutations = size(signflip,1);

        TFCEdata=zeros(possible_permutations,number_of_inputs);
        tVals=zeros(possible_permutations,number_of_inputs);


        for permIndex = 1:possible_permutations
            %print out the status
            if mod(permIndex, 100) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end

            actual_combos=signflip(permIndex,:);

            data_x_temp=data_x(logical(actual_combos),:);
            data_x_temp=cat(1,data_x_temp,data_y(logical(~actual_combos),:));
            
            data_y_temp=data_y(logical(actual_combos),:);
            data_y_temp=cat(1,data_y_temp,data_x(logical(~actual_combos),:));
            
            %Calculate T value for this grouping (ttest is paired ttest)
            [~,~,~,STATS] = ttest(data_x_temp,data_y_temp);
            tVals(permIndex,:) = STATS.tstat;
            TFCEdata(permIndex,:) = ClusterEnhancement(STATS.tstat,sparse_channel_adjacency_matrix,E,H);
        end
        clear  p STATS group1 group
        
        % calculate real t
        [~,p.real,~,REALSTATS] = ttest(data_x,data_y);
        T.tMax = sort(max(abs(tVals),[],2),'descend');
        T.chk_T = mean(abs(tVals));
        T.tMaxTFCE = sort(max(abs(TFCEdata),[],2),'descend');
        
        T.real_T = REALSTATS.tstat;
        T.real_TFCE = ClusterEnhancement(REALSTATS.tstat,sparse_channel_adjacency_matrix,E,H);


    case 'unpairedTleft'
        nSubj=size(data_x,1)+size(data_y,1);
        nGrp=size(data_y,1);
        G = nchoosek(nSubj,nGrp);

        data =  cat(1,data_y,data_x); clear data_*;

        % Choose how many relabelings to evaluate. An explicit override is
        % honoured; otherwise enumerate exactly when the group is small, else
        % fall back to a sane Monte-Carlo ceiling. snpm_relabel_assignments then
        % samples WITH REPLACEMENT (no uniqueness loop) so the run always
        % terminates, even when requested >= G (the historical hang).
        if nargin == 9
            requested_permutations = permutation_overide;
        elseif G <= 300000
            requested_permutations = G;          % exact enumeration
        else
            requested_permutations = 50000;      % default Monte-Carlo ceiling
        end

        actual_combos = snpm_relabel_assignments(nSubj,nGrp,requested_permutations);
        possible_permutations = size(actual_combos,1);

        TFCEdata=zeros(possible_permutations,number_of_inputs);
        tVals=zeros(possible_permutations,number_of_inputs);

        for permIndex = 1:possible_permutations
            %print out the status
            if mod(permIndex, 100) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end

            group1 = actual_combos(permIndex,1:nGrp);
            group2 = actual_combos(permIndex,nGrp+1:end);

            %Calculate T value for this grouping (ttest2 is unpaired ttest)
            [~,~,~,STATS] = ttest2(data(group1,:),data(group2,:),alpha,'right');
            tVals(permIndex,:) = STATS.tstat;
            TFCEdata(permIndex,:) = ClusterEnhancement(STATS.tstat,sparse_channel_adjacency_matrix,E,H);
        end
        clear  p STATS group1 group

        %calculate real ttest
        [~,p.real,~,REALSTATS] = ttest2(data(1:nGrp,:),data(nGrp+1:end,:),alpha,'right');
        T.tMax = sort(max(tVals,[],2),'descend');
        T.chk_T = mean(tVals);
        T.tMaxTFCE = sort(max(TFCEdata,[],2),'descend');
        
        T.real_T = REALSTATS.tstat;
        T.real_TFCE = ClusterEnhancement(REALSTATS.tstat,sparse_channel_adjacency_matrix,E,H);
        
    case 'unpairedTright'
        nSubj=size(data_x,1)+size(data_y,1);
        nGrp=size(data_x,1);
        G = nchoosek(nSubj,nGrp);
        data =  cat(1,data_x,data_y); clear data_*;

        % See 'unpairedTleft' for the branch logic. snpm_relabel_assignments
        % samples WITH REPLACEMENT so requested >= G can no longer hang.
        if nargin == 9
            requested_permutations = permutation_overide;
        elseif G <= 300000
            requested_permutations = G;          % exact enumeration
        else
            requested_permutations = 50000;      % default Monte-Carlo ceiling
        end

        actual_combos = snpm_relabel_assignments(nSubj,nGrp,requested_permutations);
        possible_permutations = size(actual_combos,1);

        TFCEdata=zeros(possible_permutations,number_of_inputs);
        tVals=zeros(possible_permutations,number_of_inputs);

        for permIndex = 1:possible_permutations
            %print out the status
            if mod(permIndex, 100) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end

            group1 = actual_combos(permIndex,1:nGrp);
            group2 = actual_combos(permIndex,nGrp+1:end);

            %Calculate T value for this grouping (ttest2 is unpaired ttest)
            [~,~,~,STATS] = ttest2(data(group1,:),data(group2,:),alpha,'right');
            tVals(permIndex,:) = STATS.tstat;
            TFCEdata(permIndex,:) = ClusterEnhancement(STATS.tstat,sparse_channel_adjacency_matrix,E,H);
        end
        clear  p STATS group1 group

        %calculate r ttest
        [~,p.real,~,REALSTATS] = ttest2(data(1:nGrp,:),data(nGrp+1:end,:),alpha,'right');
        T.tMax = sort(max(tVals,[],2),'descend');
        T.chk_T = mean(tVals);
        T.tMaxTFCE = sort(max(TFCEdata,[],2),'descend');
        
         T.real_T = REALSTATS.tstat;
         T.real_TFCE = ClusterEnhancement(REALSTATS.tstat,sparse_channel_adjacency_matrix,E,H);
        
    case 'unpairedTboth'
        nSubj=size(data_x,1)+size(data_y,1);
        nGrp=size(data_x,1);
        G = nchoosek(nSubj,nGrp);
        data =  cat(1,data_x,data_y); clear data_*;

        % See 'unpairedTleft' for the branch logic. snpm_relabel_assignments
        % samples WITH REPLACEMENT so requested >= G can no longer hang.
        if nargin == 9
            requested_permutations = permutation_overide;
        elseif G <= 300000
            requested_permutations = G;          % exact enumeration
        else
            requested_permutations = 50000;      % default Monte-Carlo ceiling
        end

        actual_combos = snpm_relabel_assignments(nSubj,nGrp,requested_permutations);
        possible_permutations = size(actual_combos,1);

        TFCEdata=zeros(possible_permutations,number_of_inputs);
        tVals=zeros(possible_permutations,number_of_inputs);

        for permIndex = 1:possible_permutations
            %print out the status
            if mod(permIndex, 100) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end

            group1 = actual_combos(permIndex,1:nGrp);
            group2 = actual_combos(permIndex,nGrp+1:end);

            %Calculate T value for this grouping (ttest2 is unpaired ttest)
            [~,~,~,STATS] = ttest2(data(group1,:),data(group2,:),alpha);
            tVals(permIndex,:) = STATS.tstat;
            TFCEdata(permIndex,:) = ClusterEnhancement(STATS.tstat,sparse_channel_adjacency_matrix,E,H);
        end
        clear  p STATS group1 group

        % calculate r ttest
        [~,p.real,~,REALSTATS] = ttest2(data(1:nGrp,:),data(nGrp+1:end,:));
        [T.tMax ,~] = sort(max(abs(tVals),[],2),'descend');
        T.chk_T = mean(abs(tVals));
        [T.tMaxTFCE ,~] = sort(max(abs(TFCEdata),[],2),'descend');
        
         T.real_T = REALSTATS.tstat;
         T.real_TFCE = ClusterEnhancement(REALSTATS.tstat,sparse_channel_adjacency_matrix,E,H);
        
    case 'correlationPboth'

        nSubj=size(data_x,1);
        % Choose how many row orderings to evaluate. Honour an explicit override;
        % otherwise use the sane default ceiling. snpm_row_permutations enumerates
        % exactly when requested >= nSubj! (nSubj<=10) and otherwise samples WITH
        % REPLACEMENT (no uniqueness loop) so requested >= nSubj! cannot hang.
        if nargin < 9
            requested_permutations = min(factorial(nSubj), 50000);
        else
            requested_permutations = permutation_overide;
        end

        actual_combos = snpm_row_permutations(nSubj, requested_permutations);
        possible_permutations = size(actual_combos,1);

        TFCEdata=zeros(possible_permutations,number_of_inputs);
        tVals=zeros(possible_permutations,number_of_inputs);

        % COMPLETE-COLUMN RULE (see snpm_corr_columns). A channel with any
        % missing cell is NaN in the observed map AND in every permuted map, so
        % the analysed row set cannot depend on the permutation. The old code
        % recomputed a per-channel `idxfinite` INSIDE this loop, which made the
        % number of complete pairs a function of the permutation and mis-scaled
        % the null (raw r has spread ~1/sqrt(n-1)).
        %
        % isfinite, not ~isnan: log10(0) = -Inf under datatype 'logscale' passes
        % ~isnan but still makes snpm_corr_columns return NaN, so an ~isnan mask
        % here would disagree with the statistic AND with the mask
        % core_snpm_analysis reports to the user. Same test in the Spearman arm
        % below, in snpm_cluster_analysis, and in core_snpm_analysis.
        evaluable_col = all(isfinite(data_x) & isfinite(data_y), 1);

        for permIndex = 1:possible_permutations
            %print out the status
            if mod(permIndex, 100) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end
            randompermutations = actual_combos(permIndex,:);

            % One vectorised call for the whole map; NaN columns stay NaN.
            r_corr = snpm_corr_columns(data_x(randompermutations,:), data_y);
            r_corr(~evaluable_col) = NaN;   % defensive: unreachable, the rule above already gives NaN

            tVals(permIndex,:)=r_corr;
            TFCEdata(permIndex,:) = ClusterEnhancement(r_corr,sparse_channel_adjacency_matrix,E,H);
        end
        clear  p STATS group1 group2

        % OBSERVED MAP VIA THE IDENTICAL CODE PATH AS THE PERMUTED MAPS.
        % Same snpm_corr_columns call as inside the loop, so the identity
        % relabelling (always row 1 of actual_combos) reproduces it bitwise and
        % the max-statistic ties are exact. Do not "simplify" this into a
        % per-channel corr() loop: a one-ulp difference between the observed
        % statistic and its own image under the identity relabelling can flip a
        % >= comparison against the max null, changing b by one and moving p by
        % 1/(N+1) -- 0.005 at 200 permutations, enough to cross 0.05. The
        % uncorrected parametric p still comes from corr() per channel, which is
        % a separate quantity and not scored against any permutation null.
         T.real_T = snpm_corr_columns(data_x, data_y);
         T.real_T(~evaluable_col) = NaN;
         p.real=NaN(1,number_of_inputs);
         for i = 1:number_of_inputs
             if ~evaluable_col(i), continue; end
             [~,p_tmp]=corr(data_x(:,i),data_y(:,i));
             p.real(i)=p_tmp;
         end

        T.tMax = sort(max(abs(tVals),[],2),'descend');
        T.chk_T = mean(abs(tVals));
        [T.tMaxTFCE ,~] = sort(max(abs(TFCEdata),[],2),'descend');
        T.real_TFCE = ClusterEnhancement(T.real_T,sparse_channel_adjacency_matrix,E,H);
    
    case 'correlationSboth'

            nSubj=size(data_x,1);
            % See 'correlationPboth'. snpm_row_permutations enumerates exactly when
            % requested >= nSubj! (nSubj<=10) and otherwise samples WITH REPLACEMENT
            % (no uniqueness loop), so requested >= nSubj! cannot hang.
            if nargin < 9
                requested_permutations = min(factorial(nSubj), 50000);
            else
                requested_permutations = permutation_overide;
            end

            actual_combos = snpm_row_permutations(nSubj, requested_permutations);
            possible_permutations = size(actual_combos,1);

            TFCEdata=zeros(possible_permutations,number_of_inputs);
            tVals=zeros(possible_permutations,number_of_inputs);

            % COMPLETE-COLUMN RULE (see the Pearson case and snpm_corr_columns).
            % isfinite, not ~isnan -- see the Pearson case. Read off the RAW
            % matrices, not the ranks: tiedrank maps -Inf to a finite rank, so a
            % mask taken after ranking would silently readmit a channel the
            % transform had destroyed.
            evaluable_col = all(isfinite(data_x) & isfinite(data_y), 1);

            % Spearman = Pearson on ranks. Row permutation commutes with
            % ranking, so the rank transform is hoisted OUT of the permutation
            % loop -- exact, not an approximation, and it removes a tiedrank per
            % permutation.
            rank_x = tiedrank(data_x);
            rank_y = tiedrank(data_y);

            for permIndex = 1:possible_permutations
                %print out the status
                if mod(permIndex, 100) == 0
                    disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
                end
                randompermutations = actual_combos(permIndex,:);

                r_corr = snpm_corr_columns(rank_x(randompermutations,:), rank_y);
                % NOT redundant in the Spearman arm: tiedrank maps -Inf to a
                % finite rank (unlike NaN, which it propagates), so a column made
                % non-finite by log10(0) would otherwise produce a real-looking r.
                r_corr(~evaluable_col) = NaN;

                tVals(permIndex,:)=r_corr;
                clusterTmp = ClusterEnhancement(r_corr,sparse_channel_adjacency_matrix,E,H);
                if length(clusterTmp) <  size(TFCEdata,2) % to solve same data for data_x_temp and data_y_temp corr
                    clusterTmp(end+1) = NaN; %#ok<AGROW>
                end
                TFCEdata(permIndex,:) = clusterTmp;
            end
            clear  p STATS group1 group2

            % OBSERVED MAP VIA THE IDENTICAL CODE PATH AS THE PERMUTED MAPS:
            % same rank_x/rank_y, same snpm_corr_columns call, so the identity
            % relabelling reproduces it bitwise and max-statistic ties are
            % exact. Do not replace this with corr(...,'Type','Spearman'): one
            % ulp of disagreement flips a >= against the max null and moves p by
            % 1/(N+1). The uncorrected parametric p below is a separate quantity
            % that is never scored against the null, so corr() is fine there.
             T.real_T = snpm_corr_columns(rank_x, rank_y);
             T.real_T(~evaluable_col) = NaN;
             p.real=NaN(1,number_of_inputs);
             for i = 1:number_of_inputs
                 if ~evaluable_col(i), continue; end
                 [~,p_tmp]=corr(data_x(:,i),data_y(:,i),'Type','Spearman');
                 p.real(i)=p_tmp;
             end

            T.tMax = sort(max(abs(tVals),[],2),'descend');
            T.chk_T = mean(abs(tVals));
            [T.tMaxTFCE ,~] = sort(max(abs(TFCEdata),[],2),'descend');
            T.real_TFCE = ClusterEnhancement(T.real_T,sparse_channel_adjacency_matrix,E,H);
       
        
    otherwise
        % FATAL. This used to print 'error - improproper comparison' and RETURN,
        % handing the caller a half-filled T (no real_T, no tMax, no TFCE null)
        % and an unset p, which downstream code then indexed into. A comparison
        % this engine has no permutation scheme for is not a warning condition.
        %
        % Nothing legitimate reaches here from core_snpm_analysis: it calls
        % global_stat_test with the same comparison/tail FIRST, and that function
        % errors on exactly the same set of unsupported compstrings, so any
        % combination that would land here has already aborted the run. This
        % branch protects direct callers (export_report, scripts, tests).
        error('snpm_single_threshold_with_TFCE:unsupportedComparison', ...
            ['No permutation scheme for comparison ''%s'' with tail ''%s'' (compstring ''%s''). ' ...
             'Supported here: pairedT and unpairedT (both|left|right), onesampleT (both), ' ...
             'correlationP and correlationS (both). The circular analyses ' ...
             '(circ_phase_group, circ_phase_group_u2, circ_corrAngLinear), the GLM presets ' ...
             '(anova1, ancova, regression, rmanova, mixed2way) and mixedmodel are NOT ' ...
             'compstring cases -- they run in their own pipelines and never reach this engine. ' ...
             'circ_wheeler_watson_Test and circ_WatsonsU2Test are retired.'], ...
            comparison, tail, compstring);
end
clear  random* permIndex 


% single threshold real T
critical_T_indx = floor(alpha*possible_permutations)+1;
T.critical_T = T.tMax(critical_T_indx);
p.corrected=ones(1,number_of_inputs);
p.corrected(isnan(T.real_T))=NaN;
for i = 1:length(T.real_T)
    if isnan(T.real_T(i)); continue; end
    % Minimum-bias FWE p (Nichols & Holmes 2001; Phipson & Smyth 2010):
    %   p = (#{max-null >= obs} + 1)/(N + 1)
    % applied uniformly to every permutation tier so p is never zero and the
    % observed labelling counts as one of the N permutations.
    switch tail
        case 'both'
            obs = abs(T.real_T(i));
        otherwise
            obs = T.real_T(i);
    end
    p.corrected(i) = (sum(T.tMax >= obs) + 1)/(possible_permutations + 1);
end

 % single threshold TFCE
T.critical_T_TFCE = T.tMaxTFCE(critical_T_indx);

p.correctedTFCE=ones(1,number_of_inputs);
p.correctedTFCE(isnan(T.real_TFCE))=NaN;
for i = 1:length(T.real_TFCE)
    if isnan(T.real_TFCE(i)); continue; end
    switch tail
        case 'both'
            obs = abs(T.real_TFCE(i));
        otherwise
            obs = T.real_TFCE(i);
    end
    p.correctedTFCE(i) = (sum(T.tMaxTFCE >= obs) + 1)/(possible_permutations + 1);
end
