function epoch_data = epoch_organise(config, single_session_dataset, adjusted_times, epoch_data, plot_switch)
    time_eeg = single_session_dataset.eeg_time-adjusted_times.eeg_to_emg_offset;
    data_eeg = single_session_dataset.eeg_filtered;
    
    time_tail_emg = single_session_dataset.tail_emg_time;
    data_emg = single_session_dataset.emg_raw;
    data_tail = single_session_dataset.tail_data_contami_split;
    
    time_foot = single_session_dataset.foot_time_dist_Both;
    data_foot = single_session_dataset.load_average;
    
    time_master = time_eeg;
    data_resampled_emg = interp1(time_tail_emg, data_emg, time_master, 'linear', 'extrap');
    data_resampled_tail = interp1(time_tail_emg, data_tail, time_master, 'linear', 'extrap');
    data_resampled_foot = interp1(time_foot, data_foot, time_master, 'linear', 'extrap');
    
    xlims = [20 100];
    if plot_switch
        figure
        subplot(4,1,1)
        plot(time_master, data_eeg)
        xlim(xlims)
        
        subplot(4,1,2)
        plot(time_master, data_resampled_emg)
        xlim(xlims)
        
        subplot(4,1,3)
        plot(time_master, data_resampled_tail)
        xlim(xlims)
        
        subplot(4,1,4)
        plot(time_master, data_resampled_foot)
        xlim(xlims)
        
        subplot(4,1,3)
        hold on
    end

    % expected
    epoch_times_eeg_st = epoch_data.times_eeg_st_expected;
    color_char = 'r';
    idxs_expected = plot_each_trials(config, epoch_times_eeg_st, time_master, color_char, plot_switch);

    % unexpected
    epoch_times_eeg_st = epoch_data.times_eeg_st_unexpected;
    color_char = 'b';
    idxs_unexpected = plot_each_trials(config, epoch_times_eeg_st, time_master, color_char, plot_switch);

    % right before
    epoch_times_eeg_st = epoch_data.times_eeg_st_right_before;
    color_char = 'g';
    idxs_right_before = plot_each_trials(config, epoch_times_eeg_st, time_master, color_char, plot_switch);

    % right after
    epoch_times_eeg_st = epoch_data.times_eeg_st_right_after;
    color_char = 'y';
    idxs_right_after = plot_each_trials(config, epoch_times_eeg_st, time_master, color_char, plot_switch);

    epoch_data.time_master = time_master;
    epoch_data.data_eeg = data_eeg;
    epoch_data.data_resampled_emg = data_resampled_emg;
    epoch_data.data_resampled_tail = data_resampled_tail;
    epoch_data.data_resampled_foot = data_resampled_foot;
    
    epoch_data.idxs_expected = idxs_expected;
    epoch_data.idxs_unexpected = idxs_unexpected;
    epoch_data.idxs_right_before = idxs_right_before;
    epoch_data.idxs_right_after = idxs_right_after;


end