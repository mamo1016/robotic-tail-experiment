function epochs = epochs_compute(epoch_data, idxs, plot_switch,fig_epoch, plot_color)

    a=50; b = 100;
    if plot_switch && strcmp(plot_color, 'blue')

        figure(fig_epoch)
        t = tiledlayout(4,1, 'TileSpacing', 'compact', 'Padding', 'compact');

        % Plot 1
        nexttile;
        plot(epoch_data.time_master, epoch_data.data_eeg, "Color",'black');
        ylabel('EEG'); title('Surprise ERD Analysis'); xlim([a b]);
        
        % Plot 2
        nexttile;
        plot(epoch_data.time_master, epoch_data.data_resampled_emg, "Color",'black');
        ylabel('EMG'); xlim([a b]);
        
        % Plot 3
        nexttile;
        plot(epoch_data.time_master, epoch_data.data_resampled_foot, "Color",'black');
        ylabel('Foot CoP'); xlim([a b]);
        
        % Plot 4
        nexttile;
        plot(epoch_data.time_master, epoch_data.data_resampled_tail, "Color",'black');
        ylabel('Tail Pos'); xlabel('Time (s)'); xlim([a b]);
        
        linkaxes(t.Children, 'x'); % Sync zoom
    end

    % eeg---------------------------
    epochs_data_eeg = zeros(length(idxs),idxs(1,2)-idxs(1,1)+1,size(epoch_data.data_eeg,2));
    for loop_epoch = 1:length(idxs)
        idx_a = idxs(loop_epoch,1);
        idx_b = idxs(loop_epoch,2);
        epoch_eeg = epoch_data.data_eeg(idx_a:idx_b,:);
        if plot_switch
            figure(fig_epoch)
            nexttile(1); 
            hold on;
            plot(epoch_data.time_master(idx_a:idx_b), epoch_eeg+2, "Color",plot_color)
            ylim([-2 4])
            xlim([a b])
        end
        epochs_data_eeg(loop_epoch,:,:) = epoch_eeg;
    end
    epochs = struct();
    epochs.eeg = epochs_data_eeg;


    % emg---------------------------
    epochs_data_emg = zeros(length(idxs),idxs(1,2)-idxs(1,1)+1);
    for loop_epoch = 1:length(idxs)
        idx_a = idxs(loop_epoch,1);
        idx_b = idxs(loop_epoch,2);
        epoch_emg = epoch_data.data_resampled_emg(idx_a:idx_b,:);
        if plot_switch
            nexttile(2); 
            hold on;
            plot(epoch_data.time_master(idx_a:idx_b), epoch_emg+200, "Color", plot_color)
            xlim([a b])
        end
        epochs_data_emg(loop_epoch,:) = epoch_emg;
    end
    epochs.emg = epochs_data_emg;

    % foot---------------------------
    epochs_data_foot = zeros(length(idxs),idxs(1,2)-idxs(1,1)+1);
    for loop_epoch = 1:length(idxs)
        idx_a = idxs(loop_epoch,1);
        idx_b = idxs(loop_epoch,2);
        epoch_foot = epoch_data.data_resampled_foot(idx_a:idx_b,:);
        if plot_switch
            nexttile(3); 
            hold on;
            plot(epoch_data.time_master(idx_a:idx_b), epoch_foot+2, "Color", plot_color)
            xlim([a b])
        end
        epochs_data_foot(loop_epoch,:) = epoch_foot;
    end
    epochs.foot = epochs_data_foot;

    % tail---------------------------
    epochs_data_tail = zeros(length(idxs),idxs(1,2)-idxs(1,1)+1);
    for loop_epoch = 1:length(idxs)
        idx_a = idxs(loop_epoch,1);
        idx_b = idxs(loop_epoch,2);
        epoch_tail = epoch_data.data_resampled_tail(idx_a:idx_b,:);
        if plot_switch
            nexttile(4); 
            hold on;
            plot(epoch_data.time_master(idx_a:idx_b), epoch_tail+2, "Color", plot_color)
            xlim([a b])
        end
        epochs_data_tail(loop_epoch,:) = epoch_tail;
    end
    epochs.tail = epochs_data_tail;
end