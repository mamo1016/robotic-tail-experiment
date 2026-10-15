function [expected_foot_epochs, unexpected_foot_epochs] = extract_foot_epochs_each(foot_center, adjusted_foot_time_left, eeg_event_indices, expected_trials, unexpected_trials, eeg_time, epoch_samples, config, foot_raw_dist_L,foot_time_dist_L,foot_raw_dist_R,foot_time_dist_R)
    % Function to extract foot center epochs from continuous foot pressure data
    %
    % Inputs:
    %   foot_center: Processed foot center data (left-right difference)
    %   adjusted_foot_time_left: Time vector for foot data (already adjusted for offset)
    %   eeg_event_indices: Indices of events in EEG time
    %   expected_trials: Indices of expected trials
    %   unexpected_trials: Indices of unexpected trials
    %   eeg_time: Time vector for EEG data
    %   epoch_samples: Two-element vector with [start, end] samples relative to event
    %   config.foot_fs: Foot sensor sampling frequency
    %   config.eeg_fs: EEG sampling frequency
    %
    % Outputs:
    %   expected_foot_epochs: Foot center epochs for expected trials
    %   unexpected_foot_epochs: Foot center epochs for unexpected trials
    
    % Check if foot data is available
    if isempty(foot_raw_dist_L)
        fprintf('    Warning: No foot center data available for epoching\n');
        expected_foot_epochs = [];
        unexpected_foot_epochs = [];
        return;
    end
    
    % Initialize
    n_expected = length(expected_trials);
    n_unexpected = length(unexpected_trials);
    n_samples = epoch_samples(2) - epoch_samples(1) + 1;
    samples_foot = round(n_samples * config.foot_fs / config.eeg_fs); % Adjust for different sampling rates
    
    % Calculate epoch in seconds from samples (using EEG sampling rate)
    epoch_start = epoch_samples(1) / config.eeg_fs;
    epoch_end = epoch_samples(2) / config.eeg_fs;
    
    % Adjust epoch start/end for foot sampling rate
    foot_epoch_samples = round([epoch_start, epoch_end] * config.foot_fs); % Convert from seconds to samples
    
    % Pre-allocate arrays
    expected_foot_epochs = zeros(n_expected, samples_foot, 10, 'single');
    unexpected_foot_epochs = zeros(n_unexpected, samples_foot, 10, 'single');
    
    % Process expected trials
    valid_expected_epochs = 0;
    for i = 1:n_expected
        e = expected_trials(i);
        % Get event time in EEG time domain
        event_time_eeg = eeg_time(eeg_event_indices(e));
        
        % Find closest time point in foot time (already adjusted for offset)
        [~, event_idx_foot] = min(abs(adjusted_foot_time_left - event_time_eeg));
        
        % Check if we have enough data before and after the event
        if (event_idx_foot + foot_epoch_samples(1) >= 1) && (event_idx_foot + foot_epoch_samples(2) <= length(foot_center))
            % Create time vector around the event (in seconds)
            epoch_time_range = foot_epoch_samples(1):foot_epoch_samples(2);
            
            % Extract epoch - with bounds checking
            start_idx = max(1, event_idx_foot + foot_epoch_samples(1));
            end_idx = min(length(foot_center), event_idx_foot + foot_epoch_samples(2));
            
            if end_idx > start_idx
                % epoch_data = foot_center(start_idx:end_idx);
                epoch_data = [foot_raw_dist_L(start_idx:end_idx,:), foot_raw_dist_R(start_idx:end_idx,:)];
                
                % Resample to expected length if necessary
                if length(epoch_data) ~= samples_foot
                    epoch_data = resample(double(epoch_data), samples_foot, length(epoch_data));
                end
                
                % Store epoch
                valid_expected_epochs = valid_expected_epochs + 1;
                expected_foot_epochs(valid_expected_epochs, :,:) = single(epoch_data);
            end
        end
    end
    
    % Trim array if needed
    if valid_expected_epochs < n_expected
        expected_foot_epochs = expected_foot_epochs(1:valid_expected_epochs, :, :);
    end
    
    % Process unexpected trials
    valid_unexpected_epochs = 0;
    for i = 1:n_unexpected
        e = unexpected_trials(i);
        % Get event time in EEG time domain
        event_time_eeg = eeg_time(eeg_event_indices(e));
        
        % Find closest time point in foot time (already adjusted for offset)
        [~, event_idx_foot] = min(abs(adjusted_foot_time_left - event_time_eeg));
        
        % Check if we have enough data before and after the event
        if (event_idx_foot + foot_epoch_samples(1) >= 1) && (event_idx_foot + foot_epoch_samples(2) <= length(foot_center))
            % Create time vector around the event (in seconds)
            epoch_time_range = foot_epoch_samples(1):foot_epoch_samples(2);
            
            % Extract epoch - with bounds checking
            start_idx = max(1, event_idx_foot + foot_epoch_samples(1));
            end_idx = min(length(foot_center), event_idx_foot + foot_epoch_samples(2));
            
            if end_idx > start_idx
                % epoch_data = foot_center(start_idx:end_idx);
                epoch_data = [foot_raw_dist_L(start_idx:end_idx,:), foot_raw_dist_R(start_idx:end_idx,:)];

                % Resample to expected length if necessary
                if length(epoch_data) ~= samples_foot
                    epoch_data = resample(double(epoch_data), samples_foot, length(epoch_data));
                end
                
                % Store epoch
                valid_unexpected_epochs = valid_unexpected_epochs + 1;
                unexpected_foot_epochs(valid_unexpected_epochs, :, :) = single(epoch_data);
            end
        end
    end
    
    % Trim array if needed
    if valid_unexpected_epochs < n_unexpected
        unexpected_foot_epochs = unexpected_foot_epochs(1:valid_unexpected_epochs, :, :);
    end
    
    fprintf('    Extracted %d expected and %d unexpected foot center epochs\n', ...
        size(expected_foot_epochs, 1), size(unexpected_foot_epochs, 1));
end