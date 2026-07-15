function [T,p,E,H]= snpm_single_threshold_with_TFCE_circ(data_x,data_y,neighbors,E,H,alpha,comparison,tail,permutation_overide)

%data_x, data_y = subjects x inputs

%%Optional inputs and defaults
%%alpha=0.05
%%comparison='unpairedT' ('pairedT', potentially add more)

%%Outputs
%%T struct chkT,realT,critT, tMax
%Threshold K means cluser thresholds Hi_Threshold,Med_Threshold,Low_Threshold


% Between-group differences:
% Watson-Wheeler tests: assess differences in mean preferred phases
% Watson’s two-tests: assess differences in distribution of preferred phases
%  
% Within-group correlations of mean preferred phase and memory:
% Spearman’s rank order circular-linear correlations

if nargin<4  % default E H
    E=0.5;
    H=2;
end

if nargin<6   % default alpha
    alpha = 0.05;
    display('using default alpha value of 0.05');
end


if nargin<8 % default tail
    tail = 'both';
    display('using default both tails');
end

sparse_channel_adjacency_matrix = make_neighbors_sparse(neighbors,size(neighbors,1));

number_of_inputs=size(data_x,2);
compstring=[comparison tail];
% this switch sets up the comparisons

% circ_wheeler_watson_Test, circ_cmTest, circ_corrAngLinear,circ_corrAngAng

switch compstring

    case 'circ_wheeler_watson_Testboth'
        nSubj=size(data_x,1)+size(data_y,1);
        nGrp=size(data_x,1);
        possible_permutations = nchoosek(nSubj,nGrp);
        data =  cat(1,data_x,data_y); clear data_*;

        
        if nargin == 9
            possible_permutations = permutation_overide;
        end
         if possible_permutations > 300000 || nargin == 9
            %10,000 permutations is generally considered sufficient
            
            if nargin < 9
            possible_permutations = 50000; %so it doesn't take forever to run
            end
            
            actual_combos = NaN(possible_permutations,nGrp);
            TFCEdata=zeros(possible_permutations,number_of_inputs);
            tVals=zeros(possible_permutations,number_of_inputs);
            
            for permIndex = 1:possible_permutations
                %print out the status
                if mod(permIndex, 100) == 0
                    disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
                end
                randompermutations = randperm(nSubj);
                group1 = sort(randompermutations(1:nGrp));
                
                % test whether grouping has already been used
                while ismember(group1,actual_combos,'rows')
                    randompermutations = randperm(nSubj);
                    group1 = sort(randompermutations(1:nGrp));
                end
                
                actual_combos(permIndex,:)=group1;
                group2=randompermutations(nGrp+1:end);
                
                STATS = circ_wheeler_watson_test(data(group1,:),data(group2,:));
                tVals(permIndex,:) = STATS;
                TFCEdata(permIndex,:) = ClusterEnhancement(STATS,sparse_channel_adjacency_matrix,E,H);
            end
            clear  p STATS group1 group
        else
            actual_combos = snpm_enumerate_combinations(nSubj,nGrp);
            TFCEdata=zeros(possible_permutations,number_of_inputs);
            tVals=zeros(possible_permutations,number_of_inputs);
            
            
            for permIndex = 1:possible_permutations
                %print out the status
                if mod(permIndex, 100) == 0
                    disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
                end
                
                group1 = actual_combos(permIndex,1:nGrp);
                group2 = actual_combos(permIndex,nGrp+1:end);
                
                STATS = circ_wheeler_watson_test(data(group1,:),data(group2,:));
                tVals(permIndex,:) = STATS;
                TFCEdata(permIndex,:) = ClusterEnhancement(STATS,sparse_channel_adjacency_matrix,E,H);
            end
            clear  p STATS group1 group
         end
        % calculate real value
        REALSTATS  = circ_wheeler_watson_test(data(1:nGrp,:),data(nGrp+1:end,:));

        
        T.tMax = sort(max(abs(tVals),[],2),'descend');
        T.chk_T = mean(abs(tVals));
        T.tMaxTFCE = sort(max(abs(TFCEdata),[],2),'descend');
        
        T.real_T = REALSTATS';
        T.real_TFCE = ClusterEnhancement(REALSTATS,sparse_channel_adjacency_matrix,E,H);
        %p.critical_value = chi2inv(1-alpha,2);
        
    case('circ_WatsonsU2Testboth')

        nSubj=size(data_x,1)+size(data_y,1);
        nGrp=size(data_x,1);
        possible_permutations = nchoosek(nSubj,nGrp);
        data =  cat(1,data_x,data_y); clear data_*;

        
        if nargin == 9
            possible_permutations = permutation_overide;
        end
         if possible_permutations > 300000 || nargin == 9
            %10,000 permutations is generally considered sufficient
            
            if nargin < 9
            possible_permutations = 50000; %so it doesn't take forever to run
            end
            
            actual_combos = NaN(possible_permutations,nGrp);
            TFCEdata=zeros(possible_permutations,number_of_inputs);
            tVals=zeros(possible_permutations,number_of_inputs);
            
            for permIndex = 1:possible_permutations
                %print out the status
                if mod(permIndex, 100) == 0
                    disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
                end
                randompermutations = randperm(nSubj);
                group1 = sort(randompermutations(1:nGrp));
                
                % test whether grouping has already been used
                while ismember(group1,actual_combos,'rows')
                    randompermutations = randperm(nSubj);
                    group1 = sort(randompermutations(1:nGrp));
                end
                
                actual_combos(permIndex,:)=group1;
                group2=randompermutations(nGrp+1:end);
                
                [~, STATS] = watsons_U2_approx_p(data(group1,:),data(group2,:));
                tVals(permIndex,:) = STATS;
                TFCEdata(permIndex,:) = ClusterEnhancement(STATS,sparse_channel_adjacency_matrix,E,H);
            end
            clear  STATS group1 group
        else
            actual_combos = snpm_enumerate_combinations(nSubj,nGrp);
            TFCEdata=zeros(possible_permutations,number_of_inputs);
            tVals=zeros(possible_permutations,number_of_inputs);
            
            
            for permIndex = 1:possible_permutations
                %print out the status
                if mod(permIndex, 100) == 0
                    disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
                end
                
                group1 = actual_combos(permIndex,1:nGrp);
                group2 = actual_combos(permIndex,nGrp+1:end);
                
                [~, STATS] = watsons_U2_approx_p(data(group1,:),data(group2,:));
                tVals(permIndex,:) = STATS;
                TFCEdata(permIndex,:) = ClusterEnhancement(STATS,sparse_channel_adjacency_matrix,E,H);
            end
            clear  STATS group1 group
         end
        % calculate real value
        [~, REALSTATS]  = watsons_U2_approx_p(data(1:nGrp,:),data(nGrp+1:end,:));

        
        T.tMax = sort(max(abs(tVals),[],2),'descend');
        T.chk_T = mean(abs(tVals));
        T.tMaxTFCE = sort(max(abs(TFCEdata),[],2),'descend');
        
        T.real_T = REALSTATS';
        T.real_TFCE = ClusterEnhancement(REALSTATS,sparse_channel_adjacency_matrix,E,H);



    case('circ_corrAngLinearboth')
        nSubj=size(data_x,1);
        if nargin < 9
        possible_permutations = gamma(nSubj+1);
        else
            possible_permutations = permutation_overide;
        end
        if possible_permutations > 50000 %if there are more than 9 subjects per group, run a subset (50000)
            %10,000 permutations is generally considered sufficient
            possible_permutations = 50000; %so it doesn't take forever to run
        end
        actual_combos = NaN(possible_permutations,nSubj);
        
        TFCEdata=zeros(possible_permutations,number_of_inputs);
        tVals=zeros(possible_permutations,number_of_inputs);
        
        for permIndex = 1:possible_permutations
            %print out the status
            if mod(permIndex, 1) == 0
                disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
            end
            randompermutations = randperm(nSubj);
            
            % test whether grouping has already been used
            while ismember(randompermutations,actual_combos,'rows')
                randompermutations = randperm(nSubj);
                
            end
            
            actual_combos(permIndex,:)=randompermutations;
            
            data_x_temp=data_x(randompermutations,:);
            r_corr=NaN(1,number_of_inputs);
            for i = 1:number_of_inputs
               idxfinite=~isnan(data_x_temp(:,i)) & ~isnan(data_y(:,i));

               if sum(idxfinite)~=0
                  %r_tmp=corr(data_x_temp(idxfinite,i),data_y(idxfinite,i));
                  [~, r_tmp] = circ_corrcl(data_x_temp(idxfinite,i),data_y(idxfinite,i));
                  
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
             if sum(idxfinite)~=0
                 [p_tmp,r_tmp]=circ_corrcl(data_x(idxfinite,i),data_y(idxfinite,i));           
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


% single threshold real T
critical_T_indx = floor(alpha*possible_permutations)+1;
T.critical_T = T.tMax(critical_T_indx);
p.corrected=ones(1,number_of_inputs);
p.corrected(isnan(T.real_T))=NaN;
for i = 1:length(T.real_T)
    switch tail
        case 'both'
            if find(T.tMax<abs(T.real_T(i)),1,'first')/possible_permutations
            p.corrected(i)=find(T.tMax<abs(T.real_T(i)),1,'first')/possible_permutations;
            end
        otherwise
            if T.real_T(i) >= T.tMax(end)
                p.corrected(i)=find(T.tMax<T.real_T(i),1,'first')/possible_permutations;
            end
    end
end

 % single threshold TFCE 
T.critical_T_TFCE = T.tMaxTFCE(critical_T_indx); 

p.correctedTFCE=ones(1,number_of_inputs);
p.correctedTFCE(isnan(T.real_TFCE))=NaN;
for i = 1:length(T.real_TFCE)
    switch tail
        case 'both'
            if abs(T.real_TFCE(i)) > T.tMaxTFCE(end)
                p.correctedTFCE(i)=find(T.tMaxTFCE<abs(T.real_TFCE(i)),1,'first')/possible_permutations;
            end
        otherwise
            if T.real_TFCE(i) > T.tMaxTFCE(end)
                p.correctedTFCE(i)=find(T.tMaxTFCE<T.real_TFCE(i),1,'first')/possible_permutations;
            end
    end
end
