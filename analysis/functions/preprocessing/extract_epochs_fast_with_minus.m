
%% 4. FAST EPOCHING WITH PRE-ALLOCATION
function [minus_epochs] = extract_epochs_fast_with_minus(eeg_filtered, eeg_event_indices, epoch_samples, baseline_idx, artifact_mask)
    
    fprintf('    Fast epoching with pre-allocation...\n');
    
    n_samples = epoch_samples(2) - epoch_samples(1) + 1;
    n_channels = size(eeg_filtered, 2);
    
    % Pre-allocate maximum possible size
    max_expected = length(eeg_event_indices);
    % max_unexpected = length(unexpected_trials);
    
    minus_epochs = zeros(max_expected, n_channels, n_samples, 'single');
    % unexpected_epochs = zeros(max_unexpected, n_channels, n_samples, 'single');
    
    % Vectorized epoch extraction for expected trials
    valid_expected = 0;
    for i = 1:max_expected
        % e = eeg_event_indices(i);
        event_idx = eeg_event_indices(i);
        
        % Boundary check
        if event_idx + epoch_samples(1) >= 1 && event_idx + epoch_samples(2) <= size(eeg_filtered, 1)
            epoch_range = (event_idx + epoch_samples(1)):(event_idx + epoch_samples(2));
            
            % Quick artifact check
            if ~any(artifact_mask(epoch_range))
                % Extract epoch (vectorized)
                epoch_data = eeg_filtered(epoch_range, :);
                
                % Vectorized baseline correction
                baseline_mean = mean(epoch_data(baseline_idx, :), 1);
                epoch_data = epoch_data - baseline_mean;
                
                valid_expected = valid_expected + 1;
                minus_epochs(valid_expected, :, :) = epoch_data';
            end
        end
    end
    
    % Trim to actual size
    minus_epochs = minus_epochs(1:valid_expected, :, :);
    
    
    
    fprintf('    Extracted %d minus epochs\n', size(eeg_event_indices, 1));
end