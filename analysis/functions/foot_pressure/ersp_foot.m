function [data_ersp] = ersp_foot(config, current_data)

    [num_trials, ~] = size(current_data);
    % ersp_parameters = struct();
    % Define time vector for the raw data (-2s to +4s)
    % We create a vector from -2 to +4 spaced by 1/config.eeg_fs
    % time_axis = linspace(config.epoch_start, config.epoch_end, num_samples); 
    
    %% 3. Time-Frequency Decomposition (STFT)
    % We calculate the power spectrum for each trial individually
    % Parameters for the windowing (tweak these to trade off time vs freq resolution)
    
    % Initialize accumulators
    % total_power = [];
    
    % fprintf('Processing trials...\n');
    for i = 1:num_trials
        % Extract single trial
        trial_signal = current_data(i, :);
        
        % Calculate Spectrogram
        % [s, f, t_spec] = spectrogram(x, window, noverlap, nfft, fs)
        [s, freqs, t_spec] = spectrogram(trial_signal, config.window_size, config.overlap, config.nfft, config.eeg_fs);
        
        % Calculate Power (Magnitude Squared)
        trial_power = abs(s).^2;
        
        % Accumulate sum of power (we will average later)
        % if i == 1
        %     sum_power = trial_power;
        % else
        %     sum_power = sum_power + trial_power;
        % end
    
        % Compute Mean Power across all trials
        % mean_power = sum_power / num_trials;
        
        %% 4. Baseline Normalization (The "ERSP" Step)
        % We need to convert raw power to Decibels (dB) relative to the baseline period.
        
        % Adjust spectrogram time axis to match real time (-2s start)
        % t_spec starts at roughly config.window_size/2. We shift it to align with -2s.
        % Note: The spectrogram cuts off edges, so t_spec will be shorter than 6s.
        t_real = t_spec + config.epoch_start; 
        
        % Find indices for the baseline period (-1.5s to -0.5s usually safe)
        baseline_window = config.baseline_window;
        
        base_idx = t_real >= baseline_window(1) & t_real <= baseline_window(2);
        
        % Calculate mean baseline power for EACH frequency row
        baseline_vector = mean(trial_power(:, base_idx), 2); 
        
        cmap = [linspace(0, 1, 128)', linspace(0, 1, 128)', ones(128, 1); % Blue to White
                ones(128, 1), linspace(1, 0, 128)', linspace(1, 0, 128)']; % White to Red
         
        % Apply formula: 10 * log10(Signal / Baseline)
        % We replicate the baseline vector to match the size of trial_power matrix
        single_trial_ersp_db = 10 * log10(trial_power ./ repmat(baseline_vector, 1, size(trial_power, 2)));
        if i == 1
            all_trial_ersp_db = zeros(num_trials, size(single_trial_ersp_db,1), size(single_trial_ersp_db,2));
        end
        all_trial_ersp_db(i,:,:) = single_trial_ersp_db;
    end
    data_ersp.all_trial_ersp_db = all_trial_ersp_db;
    data_ersp.time = t_real;
    data_ersp.freq = freqs;
    data_ersp.cmap = cmap;
    %% 5. Visualization
    % figure;
    % clims = [-5 5]; % Color limits in dB (Adjust based on your signal strength)
    % imagesc(t_real, freqs, all_trial_ersp_db, clims); 
    % axis xy; % Flip Y-axis so low freq is at bottom
    % colormap(jet); % 'jet' or 'parula' are standard
    % colorbar;
    % title(['ERSP Channel ' num2str(channel_to_plot)]);
    % xlabel('Time (s)');
    % ylabel('Frequency (Hz)');
    % 
    % % Add a vertical line at Time 0 (Event)
    % hold on; xline(0, 'k--', 'LineWidth', 2); hold off;
    % 
    % % Limit view to relevant frequencies (e.g., 1-40 Hz)
    % ylim([5 30]);
    % colormap(cmap);
end