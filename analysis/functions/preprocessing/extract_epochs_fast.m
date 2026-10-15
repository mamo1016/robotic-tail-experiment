
%% 4. FAST EPOCHING WITH PRE-ALLOCATION
function [expected_epochs, unexpected_epochs] = extract_epochs_fast(single_session_dataset, trial_data, phase)
        %[expected_epochs, unexpected_epochs] = extract_epochs_fast(eeg_filtered, eeg_event_indices, expected_trials, unexpected_trials, epoch_samples, baseline_idx, artifact_mask)
        
        %[expected_epochs, unexpected_epochs] = extract_epochs_fast(eeg_filtered, eeg_event_indices, expected_trials, unexpected_trials, epoch_samples, baseline_idx, artifact_mask);
        %[bef_epochs, aft_epochs] =             extract_epochs_fast(eeg_filtered, eeg_event_indices_contami, bef_trials, aft_trials, epoch_samples, baseline_idx, artifact_mask);
        
    if phase == 1
        eeg_event_indices = trial_data.eeg_event_indices;
        expected_trials = trial_data.expected_trials;
        unexpected_trials = trial_data.unexpected_trials;
    elseif phase == 2
        eeg_event_indices = trial_data.eeg_event_indices_contami;
        expected_trials = trial_data.bef_trials;
        unexpected_trials = trial_data.aft_trials;
    elseif phase == 3
        eeg_event_indices = trial_data.eeg_event_indices_contami;
        expected_trials = trial_data.expected_first_trials;
        unexpected_trials = trial_data.exsecondpected_first_trials;
    end
    fprintf('    Fast epoching with pre-allocation...\n');
    
    n_samples = trial_data.epoch_samples(2) - trial_data.epoch_samples(1) + 1;
    n_channels = size(single_session_dataset.eeg_filtered, 2);
    
    % Pre-allocate maximum possible size
    max_expected = length(expected_trials);
    max_unexpected = length(unexpected_trials);
    
    expected_epochs = zeros(max_expected, n_channels, n_samples, 'single');
    unexpected_epochs = zeros(max_unexpected, n_channels, n_samples, 'single');
    
    % Vectorized epoch extraction for expected trials
    valid_expected = 0;
    for i = 1:max_expected
        e = expected_trials(i);
        event_idx = eeg_event_indices(e);
        
        % Boundary check
        if event_idx + trial_data.epoch_samples(1) >= 1 && event_idx + trial_data.epoch_samples(2) <= size(single_session_dataset.eeg_filtered, 1)
            epoch_range = (event_idx + trial_data.epoch_samples(1)):(event_idx + trial_data.epoch_samples(2));
            
            % Quick artifact check
            if ~any(trial_data.artifact_mask(epoch_range))
                % Extract epoch (vectorized)
                epoch_data = single_session_dataset.eeg_filtered(epoch_range, :);
                
                % Vectorized baseline correction
                baseline_mean = mean(epoch_data(trial_data.baseline_idx, :), 1);
                epoch_data = epoch_data - baseline_mean;
                
                valid_expected = valid_expected + 1;
                expected_epochs(valid_expected, :, :) = epoch_data';
            end
        end
    end
    
    % Trim to actual size
    expected_epochs = expected_epochs(1:valid_expected, :, :);
    
    % Same for unexpected trials
    valid_unexpected = 0;
    for i = 1:max_unexpected
        e = unexpected_trials(i);
        event_idx = eeg_event_indices(e);
        
        if event_idx + trial_data.epoch_samples(1) >= 1 && event_idx + trial_data.epoch_samples(2) <= size(single_session_dataset.eeg_filtered, 1)
            epoch_range = (event_idx + trial_data.epoch_samples(1)):(event_idx + trial_data.epoch_samples(2));
            
            if ~any(trial_data.artifact_mask(epoch_range))
                epoch_data = single_session_dataset.eeg_filtered(epoch_range, :);
                baseline_mean = mean(epoch_data(trial_data.baseline_idx, :), 1);
                epoch_data = epoch_data - baseline_mean;
                
                valid_unexpected = valid_unexpected + 1;
                unexpected_epochs(valid_unexpected, :, :) = epoch_data';
            end
        end
    end
    
    unexpected_epochs = unexpected_epochs(1:valid_unexpected, :, :);
    
    fprintf('    Extracted %d expected and %d unexpected epochs\n', size(expected_epochs, 1), size(unexpected_epochs, 1));
end