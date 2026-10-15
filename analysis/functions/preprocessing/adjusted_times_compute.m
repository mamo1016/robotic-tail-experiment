function adjusted_times = adjusted_times_compute(single_session_dataset, adjusted_times)
    if ~isempty(single_session_dataset.foot_time_left)
        adjusted_times.adjusted_foot_time_left = single_session_dataset.foot_time_left + adjusted_times.eeg_to_emg_offset;
    else
        adjusted_times.adjusted_foot_time_left = [];
    end
     % ADD THIS: Adjust distributed foot data times
    if ~isempty(single_session_dataset.foot_time_dist_L)
        adjusted_times.adjusted_foot_time_dist_L = single_session_dataset.foot_time_dist_L + adjusted_times.eeg_to_emg_offset;
    else
        adjusted_times.adjusted_foot_time_dist_L = [];
    end
    
    if ~isempty(single_session_dataset.foot_time_dist_R)
        adjusted_times.adjusted_foot_time_dist_R = single_session_dataset.foot_time_dist_R + adjusted_times.eeg_to_emg_offset;
    else
        adjusted_times.adjusted_foot_time_dist_R = [];
    end
    
    adjusted_times.adjusted_tail_emg_time = single_session_dataset.tail_emg_time + adjusted_times.eeg_to_emg_offset;