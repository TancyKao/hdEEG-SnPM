function [Clusters] = snpm_cluster_analysis_circ(data_x,data_y,threshold,neighbors,alpha,comparison,tail,permutation_overide)

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

if nargin < 5
    alpha = 0.05;
    disp('using default alpha value of 0.05');
end


if nargin < 7
    tail = 'both';
    disp('using default both tails');
    
end


thresholds = threshold;

number_of_inputs = size(data_x,2);
compstring=[comparison tail];

switch compstring

    case 'circ_wheeler_watson_Testboth'
        
        nSubj=size(data_x,1)+size(data_y,1);
        nGrp=size(data_x,1);
        possible_permutations = nchoosek(nSubj,nGrp);
        data =  cat(1,data_x,data_y); clear data_*;
        
        if nargin == 8
            possible_permutations = permutation_overide;
        end
        
        if possible_permutations > 300000 || nargin == 8 %if there are more than 9 subjects per group, run a subset (50000)
            
            if nargin < 8
                possible_permutations = 50000;
            end
            
            actual_combos = NaN(possible_permutations,nGrp);
            
            max_cluster_sizes = zeros(possible_permutations,length(thresholds));
            
            for permIndex=1:possible_permutations
                %print out the status
                if mod(permIndex, 1000) == 0
                    disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
                end
                % arranges the subject number in different orders and then takes the first grouping and the second group
                randompermutations = randperm(nSubj);
                group1 = sort(randompermutations(1:nGrp));

                %test whether grouping has already been used
                while ismember(group1,actual_combos,'rows')
                    randompermutations = randperm(nSubj);
                    group1 = sort(randompermutations(1:nGrp));
                end
                
                actual_combos(permIndex,:)=group1;
                group2=randompermutations(nGrp+1:end);
                
                
                STATS = circ_wheeler_watson_test(data(group1,:),data(group2,:));
               

                % find clusters above threshold
                for ti = 1:length(thresholds)
                    temp_clusters = snpm_find_clusters_graphalgs(STATS,thresholds(ti),sparse_channel_adjacency_matrix);
                    if ~isempty(temp_clusters)
                        max_cluster_sizes(permIndex,ti) = max(cellfun('length',temp_clusters));
                    end
                end
            end
        else
            actual_combos = snpm_enumerate_combinations(nSubj,nGrp);
            max_cluster_sizes = zeros(possible_permutations,length(thresholds));
             for permIndex=1:possible_permutations
                %print out the status
                if mod(permIndex, 1000) == 0
                    disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
                end
                
                % arranges the subject number in different orders and then takes the first grouping and the second group
                group1 = actual_combos(permIndex,1:nGrp);
                group2 = actual_combos(permIndex,nGrp+1:end);
                
 
                STATS  = circ_wheeler_watson_test(data(group1,:),data(group2,:));

                % calculate clusters above threshold
                for ti = 1:length(thresholds)
                    temp_clusters = snpm_find_clusters_graphalgs(STATS,thresholds(ti),sparse_channel_adjacency_matrix);
                    if ~isempty(temp_clusters)
                        max_cluster_sizes(permIndex,ti) = max(cellfun('length',temp_clusters));
                    end
                end
            end
        end
        
        % real t
        REALSTATS = circ_wheeler_watson_test(data(1:nGrp,:),data(nGrp+1:end,:));  


    case 'circ_WatsonsU2Testboth'
        
        nSubj=size(data_x,1)+size(data_y,1);
        nGrp=size(data_x,1);
        possible_permutations = nchoosek(nSubj,nGrp);
        data =  cat(1,data_x,data_y); clear data_*;
        
        if nargin == 8
            possible_permutations = permutation_overide;
        end
        
        if possible_permutations > 300000 || nargin == 8 %if there are more than 9 subjects per group, run a subset (50000)
            
            if nargin < 8
                possible_permutations = 50000;
            end
            
            actual_combos = NaN(possible_permutations,nGrp);
            
            max_cluster_sizes = zeros(possible_permutations,length(thresholds));
            
            for permIndex=1:possible_permutations
                %print out the status
                if mod(permIndex, 1000) == 0
                    disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
                end
                % arranges the subject number in different orders and then takes the first grouping and the second group
                randompermutations = randperm(nSubj);
                group1 = sort(randompermutations(1:nGrp));

                %test whether grouping has already been used
                while ismember(group1,actual_combos,'rows')
                    randompermutations = randperm(nSubj);
                    group1 = sort(randompermutations(1:nGrp));
                end
                
                actual_combos(permIndex,:)=group1;
                group2=randompermutations(nGrp+1:end);
                
                [~, STATS] = watsons_U2_approx_p(data(group1,:),data(group2,:));
  
                % find clusters above threshold
                for ti = 1:length(thresholds)
                    temp_clusters = snpm_find_clusters_graphalgs(STATS,thresholds(ti),sparse_channel_adjacency_matrix);
                    if ~isempty(temp_clusters)
                        max_cluster_sizes(permIndex,ti) = max(cellfun('length',temp_clusters));
                    end
                end
            end
        else
            actual_combos = snpm_enumerate_combinations(nSubj,nGrp);
            max_cluster_sizes = zeros(possible_permutations,length(thresholds));
             for permIndex=1:possible_permutations
                %print out the status
                if mod(permIndex, 1000) == 0
                    disp([num2str(permIndex),' out of ',num2str(possible_permutations),' combinations completed...']);
                end
                
                % arranges the subject number in different orders and then takes the first grouping and the second group
                group1 = actual_combos(permIndex,1:nGrp);
                group2 = actual_combos(permIndex,nGrp+1:end);
                
                %Calculate T value for this grouping (ttest2 is unpaired ttest)
                %[~,~,~,REALSTATS] = ttest2(data(group1,:),data(group2,:));
                [~, STATS]  = watsons_U2_approx_p(data(group1,:),data(group2,:));

                % calculate clusters above threshold
                for ti = 1:length(thresholds)
                    temp_clusters = snpm_find_clusters_graphalgs(STATS,thresholds(ti),sparse_channel_adjacency_matrix);
                    if ~isempty(temp_clusters)
                        max_cluster_sizes(permIndex,ti) = max(cellfun('length',temp_clusters));
                    end
                end
            end
        end
        
        % real t
        [~,REALSTATS] = watsons_U2_approx_p(data(1:nGrp,:),data(nGrp+1:end,:));  

 case 'circ_corrAngLinearboth'
        
        nSubj=size(data_x,1);
        possible_permutations = gamma(nSubj+1);
        
        if possible_permutations > 50000 %if there are more than 9 subjects per group, run a subset (50000)
            %10,000 permutations is generally considered sufficient
            possible_permutations = 50000; %so it doesn't take forever to run
        end
        
        if nargin == 8
            possible_permutations = permutation_overide;
        end
        
        actual_combos = NaN(possible_permutations,nSubj);
        
        max_cluster_sizes = zeros(possible_permutations,length(thresholds));
        
        for permIndex=1:possible_permutations
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
            r_corr   = NaN(1,number_of_inputs);
            nr2_corr = NaN(1,number_of_inputs);
            for i = 1:number_of_inputs
               idxfinite=~isnan(data_x_temp(:,i)) & ~isnan(data_y(:,i));

               if sum(idxfinite)~=0
                  
                  [~, r_tmp, nr2_tmp] = circ_corrcl(data_x_temp(idxfinite,i),data_y(idxfinite,i));
                  
                  r_corr(i)   = r_tmp;
                  nr2_corr(i) = nr2_tmp;

               else
                  r_corr(i)   = NaN;
                  nr2_corr(i) = NaN;
               end                
 
            end

            
            for ti = 1:length(thresholds)
                temp_clusters = snpm_find_clusters_graphalgs(nr2_corr,thresholds(ti),sparse_channel_adjacency_matrix);
                if ~isempty(temp_clusters)
                    max_cluster_sizes(permIndex,ti) = max(cellfun('length',temp_clusters));
                end
            end
        end

        % calculate real correlation values
         REALSTATS=NaN(1,number_of_inputs);
         p=NaN(1,number_of_inputs);
         for i = 1:number_of_inputs
             idxfinite=~isnan(data_x(:,i)) & ~isnan(data_y(:,i));
             if sum(idxfinite)~=0
                 [p_tmp,r_tmp, chi2_tmp]=circ_corrcl(data_x(idxfinite,i),data_y(idxfinite,i));           
                 REALSTATS_r(i)=r_tmp;
                 REALSTATS_chi2(i) = chi2_tmp;
                 p(i)=p_tmp;
             else
                 REALSTATS(i)=NaN;
                 REALSTATS_chi2(i) = NaN;
                 p(i)=NaN;  
             end
         end

    otherwise
        disp('ERROR - improproper comparison')
end
clear ti temp_clusters group*

%%%%%%
cl=1;
Clusters=struct([]);
for ti = 1:length(thresholds)
    temp_clusters = snpm_find_clusters_graphalgs(REALSTATS_chi2,thresholds(ti),sparse_channel_adjacency_matrix);
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
if strcmp(tail,'both') && length(threshold) == 1
    max_cluster_sizes =  sort(max(max_cluster_sizes,[],2),1,'descend');
    for cli = 1 : length(Clusters)
        potential_clusters=Clusters(cli).channels;
        Clusters(cli).p=find(max_cluster_sizes<length(potential_clusters),1,'first')/possible_permutations;
        if isempty(Clusters(cli).p)
            Clusters(cli).p = 1;
        end
        Clusters(cli).permutations=possible_permutations;
    end
    
else
    max_cluster_sizes =  sort(max_cluster_sizes,1,'descend');
    for cli = 1 : length(Clusters)
        potential_clusters=Clusters(cli).channels;
        ti = Clusters(cli).threshold==thresholds;
        Clusters(cli).p=find(max_cluster_sizes(:,ti)<length(potential_clusters),1,'first')/possible_permutations;
        if isempty(Clusters(cli).p)
            Clusters(cli).p = 1;
        end
        Clusters(cli).permutations=possible_permutations;
    end
    
end




