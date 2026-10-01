function out_put = foot_ersp_t_test(config_para)


num_subs = length(config_para.name_list);
trial_size = config_para.trial_size;
out_put = struct();
for sub_num=1:num_subs
    for i = 1:3
        if trial_size(i,sub_num)==0
            continue;
        end
        filename = string("output\foot\foot_data_dif_bet_each_ex_sub" + sub_num +"_session" + i + "_.mat");        
        dif_bet_each_ex = load(filename).dif_bet_each_ex;
    
        if sub_num == 1 && i == 1
            num_tiles = length(dif_bet_each_ex.t)*length(dif_bet_each_ex.f);
            config_para.tile_t = dif_bet_each_ex.t;
            config_para.tile_f = dif_bet_each_ex.f;
            save("config_para.mat", 'config_para', '-v7.3');
            common_size = size(dif_bet_each_ex.matrix,1);
            all_sub_diffs_mat = zeros(num_subs*common_size*3, num_tiles);
            % all_sub_diffs_mat_store = zeros(3, num_subs, common_size, dif_bet_each_ex.f, dif_bet_each_ex.t);
            idx=1;
        end
        % Initialize Matrix: [25 Subjects x Num_Tiles]
        
                
        for trial_no = 1:common_size
            all_sub_diffs_mat(idx, :) = reshape(dif_bet_each_ex.matrix(trial_no,:,:), [1, num_tiles]);
            % all_sub_diffs_mat_store(i,sub_num,trial_no,:,:) = dif_bet_each_ex.matrix(trial_no,:,:);
            idx = idx + 1;
        end
    end
end
all_sub_diffs_mat(idx:end, :) = [];
%% 2. Perform One-Sample T-Test (The "Drill")
% We test every tile (column) against 0.
% Alpha = 0.05 means 95% confidence.
[h, p_values, ci, stats] = ttest(all_sub_diffs_mat, 0, 'Alpha', 0.05);
%% 3. Results & Filtering
sig_indices = h == 0; % Get the index numbers of not significant tiles
% sig_t_vals  = stats.tstat(find(h == 1)); % Get the T-values for those tiles
all_sub_diffs_mat(:,sig_indices) = 0;
out_put.all_sub_diffs_mat = all_sub_diffs_mat;
out_put.sig_indices = sig_indices;
% out_put.all_sub_diffs_mat_store;