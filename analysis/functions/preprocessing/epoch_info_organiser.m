function [trial_data,config,single_session_data] = epoch_info_organiser(trial_data,config,single_session_dataset,adjusted_times)
    for loop_num = 1:2 % loop 1: expected, loop 2: unexpected
        if loop_num == 1
            indices_tail_emg = trial_data.expected_indices_tail_emg;
        elseif loop_num == 2
            indices_tail_emg = trial_data.unexpected_indices_tail_emg;
        end
        epochs_tail_emg=zeros(length(indices_tail_emg), config.emg_fs*config.epoch_end-config.emg_fs*config.epoch_start+1);
        for i = 1:length(indices_tail_emg)
            % tail plot
            idx_st=indices_tail_emg(i)+config.emg_fs*config.epoch_start;
            idx_ed=indices_tail_emg(i)+config.emg_fs*config.epoch_end;    
            epochs_tail_emg(i,:) = single_session_dataset.tail_data_contami_split(idx_st:idx_ed);   
        
        
        end
        
        
        if loop_num == 1
            trial_data.expected_epoch_times_tail_emg_st = single_session_dataset.tail_emg_time(indices_tail_emg);
            trial_data.expected_epochs_tail_emg = epochs_tail_emg;
            trial_data.expected_epoch_times_eeg_st = single_session_dataset.eeg_time(trial_data.expected_indices_eeg)-adjusted_times.eeg_to_emg_offset;

            trial_data.expected_epoch_times_foot_st = single_session_dataset.foot_time_dist_Both(trial_data.expected_indices_foot)-adjusted_times.eeg_to_emg_offset;
        elseif loop_num == 2
            trial_data.unexpected_epoch_times_tail_emg_st = single_session_dataset.tail_emg_time(indices_tail_emg);
            trial_data.unexpected_epochs_tail_emg = epochs_tail_emg;
            trial_data.unexpected_epoch_times_eeg_st = single_session_dataset.eeg_time(trial_data.unexpected_indices_eeg)-adjusted_times.eeg_to_emg_offset;

            trial_data.expected_epoch_times_foot_st = single_session_dataset.foot_time_dist_Both(trial_data.expected_indices_foot)-adjusted_times.eeg_to_emg_offset;
        end    
    end

    config = orderfields(config);
    single_session_data = orderfields(single_session_dataset);
    trial_data = orderfields(trial_data);


    
end