
function [expected_emg_epochs, unexpected_emg_epochs] = extract_emg_epochs(emg_raw, emg_time, eeg_event_indices, expected_trials, unexpected_trials, eeg_time, epoch_samples, eeg_to_emg_offset, config)
    % Function to extract EMG epochs from continuous EMG data
    %
    % Inputs:
    %   emg_raw: Raw EMG data
    %   emg_time: Time vector for EMG data
    %   eeg_event_indices: Indices of events in EEG time
    %   expected_trials: Indices of expected trials
    %   unexpected_trials: Indices of unexpected trials
    %   eeg_time: Time vector for EEG data
    %   epoch_samples: Two-element vector with [start, end] samples relative to event
    %   eeg_to_emg_offset: Time offset between EEG and EMG recordings
    %   config.emg_fs: EMG sampling frequency
    %   config.eeg_fs: EEG sampling frequency
    %
    % Outputs:
    %   expected_emg_epochs: EMG epochs for expected trials
    %   unexpected_emg_epochs: EMG epochs for unexpected trials
    
    % Initialize
    n_expected = length(expected_trials);
    n_unexpected = length(unexpected_trials);
    n_samples = epoch_samples(2) - epoch_samples(1) + 1;
    samples_emg = round(n_samples * config.emg_fs / config.eeg_fs); % Adjust for different sampling rates
    
    % Adjust epoch start/end for EMG sampling rate
    emg_epoch_samples = round(epoch_samples * config.emg_fs / config.eeg_fs);
    
    % Pre-allocate arrays
    expected_emg_epochs = zeros(n_expected, samples_emg, 'single');
    unexpected_emg_epochs = zeros(n_unexpected, samples_emg, 'single');
    
    % Convert EEG event times to EMG event times using offset
    valid_expected_epochs = 0;
    for i = 1:n_expected
        e = expected_trials(i);
        % Get event time in EEG time domain
        event_time_eeg = eeg_time(eeg_event_indices(e));
        
        % Convert to EMG time domain (accounting for offset)
        event_time_emg = event_time_eeg - eeg_to_emg_offset;
        
        % Find closest time point in EMG time
        [~, event_idx_emg] = min(abs(emg_time - event_time_emg));
        
        % Check if epoch is within boundaries
        if (event_idx_emg + emg_epoch_samples(1) >= 1) && (event_idx_emg + emg_epoch_samples(2) <= length(emg_raw))
            % Extract epoch
            epoch_range = (event_idx_emg + emg_epoch_samples(1)):(event_idx_emg + emg_epoch_samples(2));
            
            % Resample to match expected output size if necessary
            emg_epoch = emg_raw(epoch_range);
            if length(emg_epoch) ~= samples_emg
                emg_epoch = resample(double(emg_epoch), samples_emg, length(emg_epoch));
            end
            
            % Store epoch
            valid_expected_epochs = valid_expected_epochs + 1;
            expected_emg_epochs(valid_expected_epochs, :) = single(emg_epoch);
        end
    end
    
    % Trim array if needed
    if valid_expected_epochs < n_expected
        expected_emg_epochs = expected_emg_epochs(1:valid_expected_epochs, :);
    end
    
    % Do the same for unexpected trials
    valid_unexpected_epochs = 0;
    for i = 1:n_unexpected
        e = unexpected_trials(i);
        % Get event time in EEG time domain
        event_time_eeg = eeg_time(eeg_event_indices(e));
        
        % Convert to EMG time domain (accounting for offset)
        event_time_emg = event_time_eeg - eeg_to_emg_offset;
        
        % Find closest time point in EMG time
        [~, event_idx_emg] = min(abs(emg_time - event_time_emg));
        
        % Check if epoch is within boundaries
        if (event_idx_emg + emg_epoch_samples(1) >= 1) && (event_idx_emg + emg_epoch_samples(2) <= length(emg_raw))
            % Extract epoch
            epoch_range = (event_idx_emg + emg_epoch_samples(1)):(event_idx_emg + emg_epoch_samples(2));
            
            % Resample to match expected output size if necessary
            emg_epoch = emg_raw(epoch_range);
            if length(emg_epoch) ~= samples_emg
                emg_epoch = resample(double(emg_epoch), samples_emg, length(emg_epoch));
            end
            
            % Store epoch
            valid_unexpected_epochs = valid_unexpected_epochs + 1;
            unexpected_emg_epochs(valid_unexpected_epochs, :) = single(emg_epoch);
        end
    end
    
    % Trim array if needed
    if valid_unexpected_epochs < n_unexpected
        unexpected_emg_epochs = unexpected_emg_epochs(1:valid_unexpected_epochs, :);
    end
    
    fprintf('    Extracted %d expected and %d unexpected EMG epochs\n', ...
        size(expected_emg_epochs, 1), size(unexpected_emg_epochs, 1));
end