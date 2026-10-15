function idx = plot_each_trials(config, epoch_times_eeg_st, time_master, color, plot_switch)

    idx_a = zeros(length(epoch_times_eeg_st),1);
    idx_b = zeros(length(epoch_times_eeg_st),1);

    for i = 1:length(epoch_times_eeg_st)
        a = epoch_times_eeg_st(i)+config.epoch_start;
        % b = epoch_times_eeg_st(i)+config.epoch_end;
        [~, idx_a(i)] = min(abs(time_master-a));
        idx_b(i)= idx_a(i) + config.eeg_fs*(config.epoch_end-config.epoch_start);
        % [~, idx_b] = min(abs(time_master-b));
        if idx_b(i)<length(time_master)
            y = zeros(length(time_master(idx_a(i):idx_b(i))),1);
            if rem(i,2) == 0
                y = y + 3.5;
            else
                y = y + 3.3;
            end
    
            if plot_switch == 1
                plot(time_master(idx_a(i):idx_b(i)), y,color)
                ylim([-1.5 4])
            end
        else
            break;
        end
    end

    if idx_b(i)>length(time_master)
        idx = [idx_a(1:i-1), idx_b(1:i-1)];
    else
        idx = [idx_a, idx_b];
    end