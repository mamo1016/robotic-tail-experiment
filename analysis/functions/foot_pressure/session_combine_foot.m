function foot_data_boot_combined = session_combine_foot(foot_cleaned)
name_list = fieldnames(foot_cleaned);
subject_number_list = zeros(length(name_list),1);

for struct_num = 1:length(name_list)
        current_name = name_list{struct_num};           % Get the name (string)    
        trial_type = fieldnames(foot_cleaned.(current_name)); 
        tokens = regexp(current_name, 'subject_(\d+)', 'tokens'); % get sub num and ses num    
        subject_number_list(struct_num) = str2double(tokens{1,1}{1});
end

for loop = 1:1
    unique_num = unique(subject_number_list);
    for struct_num = 1:length(unique_num)
        
        for trial_type_num = 1: length(trial_type)
            trial_type_name = trial_type{trial_type_num};
            footA = [];
            footB = [];
            footC = [];
            
            for i = 1:length(find(subject_number_list==unique_num(struct_num)))
                current_sub_num = find(subject_number_list==unique_num(struct_num));
                current_name = name_list{current_sub_num(i)};           % Get the name (string)  
                
                % ersp_session_combined.(string("subject"+tokens{1,1}{1})) = 
                foot_data_boot = foot_cleaned.(current_name).(trial_type_name);
                if i == 1
                    footA = foot_data_boot;
                elseif i == 2
                    footB = foot_data_boot;
                elseif i == 3
                    footC = foot_data_boot;
                end
            end 
            tokens = regexp(current_name, 'subject_(\d+)', 'tokens'); % get sub num and ses num
            foot_data_boot_combined.(string("subject"+tokens{1,1}{1})).(trial_type_name) = cat(1, footA, footB, footC);
            length_store.(string("subject"+tokens{1,1}{1})).(trial_type_name) = [size(footA); size(footB); size(footC)];
            % foot_data_boot = foot_data_boot_combined.(string("subject"+tokens{1,1}{1})).(trial_type_name);
            % rng('shuffle');

        end
    end
    
    % loop_3dig = sprintf('%03d', loop); 
end
filename = "output\foot\foot_cleaned_session_combined.mat";        
save(filename, 'foot_data_boot_combined', '-v7.3');

filename = "output\foot\length_store.mat";        
save(filename, 'length_store', '-v7.3');