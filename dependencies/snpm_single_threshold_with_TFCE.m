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

        for permIndex = 1:possible_permutations
            %print out the status
            if mod(permIndex, 1) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end
            randompermutations = actual_combos(permIndex,:);

            data_x_temp=data_x(randompermutations,:);
            r_corr=NaN(1,number_of_inputs);
            for i = 1:number_of_inputs
               idxfinite=~isnan(data_x_temp(:,i)) & ~isnan(data_y(:,i));
%                 [r_tmp,~]=corrcoef(data_x_temp(idxfinite,i),data_y(idxfinite,i));
%                 r_corr(i)=r_tmp(1,2);
               if sum(idxfinite)~=0
                  r_tmp=corr(data_x_temp(idxfinite,i),data_y(idxfinite,i));
                  r_corr(i)=r_tmp;

               else
                  r_corr(i)=NaN;
               end                
 
            end
            tVals(permIndex,:)=r_corr;
            TFCEdata(permIndex,:) = ClusterEnhancement(r_corr,sparse_channel_adjacency_matrix,E,H);
        end
        clear  p STATS group1 group2
        
        % calculate real correlation values
         T.real_T=NaN(1,number_of_inputs);
         p.real=NaN(1,number_of_inputs);
         for i = 1:number_of_inputs
             idxfinite=~isnan(data_x(:,i)) & ~isnan(data_y(:,i));
             %[r_tmp,p_tmp]=corrcoef(data_x(idxfinite,i),data_y(idxfinite,i));
             if sum(idxfinite)~=0
                 [r_tmp,p_tmp]=corr(data_x(idxfinite,i),data_y(idxfinite,i));           
                 T.real_T(i)=r_tmp;
                 p.real(i)=p_tmp;
             else
                 T.real_T(i)=NaN;
                 p.real(i)=NaN;  
             end
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

            for permIndex = 1:possible_permutations
                %print out the status
                if mod(permIndex, 1) == 0
                    disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
                end
                randompermutations = actual_combos(permIndex,:);

                data_x_temp=data_x(randompermutations,:);
                r_corr=NaN(1,number_of_inputs);
                for i = 1:number_of_inputs
                   idxfinite=~isnan(data_x_temp(:,i)) & ~isnan(data_y(:,i));
    %                 [r_tmp,~]=corrcoef(data_x_temp(idxfinite,i),data_y(idxfinite,i));
    %                 r_corr(i)=r_tmp(1,2);
                   if sum(idxfinite)~=0
                      r_tmp=corr(data_x_temp(idxfinite,i),data_y(idxfinite,i),'Type','Spearman','rows','complete');
                      r_corr(i)=r_tmp;
    
                   else
                      r_corr(i)=NaN;
                   end                
     
                end
                tVals(permIndex,:)=r_corr;
                clusterTmp = ClusterEnhancement(r_corr,sparse_channel_adjacency_matrix,E,H);
                if length(clusterTmp) <  size(TFCEdata,2) % to solve same data for data_x_temp and data_y_temp corr 
                    clusterTmp(end+1) = NaN; 
                end
                TFCEdata(permIndex,:) = clusterTmp;
            end
            clear  p STATS group1 group2
            
            % calculate real correlation values
             T.real_T=NaN(1,number_of_inputs);
             p.real=NaN(1,number_of_inputs);
             for i = 1:number_of_inputs
                 idxfinite=~isnan(data_x(:,i)) & ~isnan(data_y(:,i));
                 %[r_tmp,p_tmp]=corrcoef(data_x(idxfinite,i),data_y(idxfinite,i));
                 if sum(idxfinite)~=0
                     [r_tmp,p_tmp]=corr(data_x(idxfinite,i),data_y(idxfinite,i),'Type','Spearman','rows','complete');           
                     T.real_T(i)=r_tmp;
                     p.real(i)=p_tmp;
                 else
                     T.real_T(i)=NaN;
                     p.real(i)=NaN;  
                 end
             end
             
            T.tMax = sort(max(abs(tVals),[],2),'descend');
            T.chk_T = mean(abs(tVals));
            [T.tMaxTFCE ,~] = sort(max(abs(TFCEdata),[],2),'descend');
            T.real_TFCE = ClusterEnhancement(T.real_T,sparse_channel_adjacency_matrix,E,H);
       
        
    otherwise
        display('error - improproper comparison')
        return;
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
