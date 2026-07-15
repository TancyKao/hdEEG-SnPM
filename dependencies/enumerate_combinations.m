% written by Brian Hamm 4/11/2014
% works for symmetric, asymmetric - Brian Hamm

function combinations = enumerate_combinations(num_subjects, group_1_size)
    group_2_size      = num_subjects - group_1_size;
    num_combinations  = nchoosek(num_subjects, group_1_size);
    combinations      = zeros(num_combinations, num_subjects);
    combination       = 1:num_subjects;
    
    for combination_num = 1:num_combinations
        group_2_index = group_1_size + 1;
        subject_num   = 1;
        
        for group_1_index = 1:group_1_size
            while subject_num < combination(group_1_index)
                combination(group_2_index) = subject_num;
                group_2_index              = group_2_index + 1;
                subject_num                = subject_num + 1;
            end
            
            subject_num = combination(group_1_index) + 1;
        end
        
        for group_2_index = group_2_index:num_subjects
            combination(group_2_index) = group_2_index;
        end
        
        combinations(combination_num, :) = combination;
        
        for index = group_1_size:-1:1
            combination(index) = combination(index) + 1;
            
            if combination(index) <= group_2_size + index
                for broken_index = (index + 1):group_1_size
                    combination(broken_index) = combination(broken_index - 1) + 1;
                end
                
                break;
            end
        end
    end
end
