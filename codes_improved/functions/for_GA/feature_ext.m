function [eeg_feature, foot_feature] = feature_ext(eeg_ersp_data, foot_ersp_data)

    sub_name = "sub_1";
    % --- 1. Mock Data Generation (Replace this with your real X and Y) ---
    % X: 25 Participants x 1200 "Tiles" (Time-Freq bins)
    % Y: 25 Participants x 1200 "Tiles" (Time-Freq bins)
    sub_list = fieldnames(eeg_ersp_data); % subject name list
    num_subs = length(sub_list); % number of subjects
    
    types_eeg = fieldnames(eeg_ersp_data.(sub_name)); 
    types_foot = fieldnames(foot_ersp_data.(sub_name));
    
    % type_eeg = types_eeg{1}; % here change
    % type_foot = name_finder(types_foot, type_eeg);
    
    eeg_feature = zeros(num_subs,numel(eeg_ersp_data.(sub_list{1}).expected.vector)*3);
    foot_feature = zeros(num_subs,numel(eeg_ersp_data.(sub_list{1}).expected.vector)*3);
    
    for sub_num = 1:num_subs 
        for i = 1:3
            switch i
                case 1
                    eeg_str = ["expected_duplicated", "expected"]; % 1 - 2
                case 2
                    eeg_str = ["unexpected", "expected"]; % 1 - 2
                case 3
                    eeg_str = ["right_after", "right_before"]; % 1 - 2
            end
            foot_str = [string(name_finder(eeg_str(1))), string(name_finder(eeg_str(2)))];
            eeg_feature(sub_num,1+(i-1)*400:400*i)  = reshape(eeg_ersp_data.(sub_list{sub_num}).(eeg_str(1)).vector   -  eeg_ersp_data.(sub_list{sub_num}).(eeg_str(2)).vector, 1, []);
            foot_feature(sub_num,1+(i-1)*400:400*i) = reshape(foot_ersp_data.(sub_list{sub_num}).(foot_str(1)).vector -  foot_ersp_data.(sub_list{sub_num}).(foot_str(2)).vector, 1, []);
        end
    end
end