function out_put = tile_t_test(all_sub_diffs,sub_list,types,type_name_order,str)
    num_tiles = size(all_sub_diffs.(sub_list{1}).(str).(string(types{type_name_order(2*1-1)} + "_vs_" + types{type_name_order(2*1)})));
    % Initialize Matrix: [25 Subjects x Num_Tiles]
    num_subs = length(sub_list);
    all_sub_diffs_mat = zeros(num_subs, num_tiles(1)*num_tiles(2));
    out_put = struct();
    for j = 1:3
            
        for i = 1:num_subs
    
            sub_name = sub_list{i};
            eeg_type_A = types{type_name_order(2*j-1)};
            eeg_type_B = types{type_name_order(2*j)};
    
            all_sub_diffs_mat(i,:) = reshape(all_sub_diffs.(sub_name).(str).(string(eeg_type_A + "_vs_" + eeg_type_B)),size(all_sub_diffs_mat(i,:)));
        end
        %% 2. Perform One-Sample T-Test (The "Drill")
        % We test every tile (column) against 0.
        % Alpha = 0.05 means 95% confidence.
        [h, p_values, ci, stats] = ttest(all_sub_diffs_mat, 0, 'Alpha', 0.05);
        %% 3. Results & Filtering
        sig_indices = h == 0; % Get the index numbers of not significant tiles
        % sig_t_vals  = stats.tstat(find(h == 1)); % Get the T-values for those tiles
    
        for i = 1:num_subs
            
            sub_name = sub_list{i};
            eeg_type_A = types{type_name_order(2*j-1)};
            eeg_type_B = types{type_name_order(2*j)};
            all_sub_diffs_mat(i,sig_indices) = 0;
            tile_size = size(all_sub_diffs.(sub_name).(str).(string(eeg_type_A + "_vs_" + eeg_type_B)));
            out_put.(sub_name).(string(str+"_significant")).(string(eeg_type_A + "_vs_" + eeg_type_B)) = reshape(all_sub_diffs_mat(i,:), tile_size);
            out_put.(sub_name).(string(str+"_significant")).(string(eeg_type_A + "_vs_" + eeg_type_B + "__p_val")) = reshape(p_values, tile_size);
            out_put.(sub_name).(string(str+"_significant")).(string(eeg_type_A + "_vs_" + eeg_type_B + "__stats")) = stats;
            
        end
    end
end