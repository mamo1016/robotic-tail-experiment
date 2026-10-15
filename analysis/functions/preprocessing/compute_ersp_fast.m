%% 5. MUCH FASTER TIME-FREQUENCY ANALYSIS
function [expected_ersp, unexpected_ersp, ersp_times, freqs] = compute_ersp_fast(first_group_epochs, second_group_epochs, config, eeg_epoch_time, boot_strap_len)
    fprintf('    Computing time-frequency analysis (optimized)...\n');
    
    % Use shorter time-frequency analysis focused on region of interest
    roi_time_idx = find(eeg_epoch_time >= -1.0 & eeg_epoch_time <= 3.0); % Focus on relevant time
    roi_freqs = config.freq_range(config.freq_range >= 10 & config.freq_range <= 35); % Focus on relevant frequencies
    
    % Reduce bootstrap iterations for speed (still statistically valid)
    fast_bootstrap = min(boot_strap_len, 1000); % Cap at 1000 for speed
    
    % Use optimized STFT parameters for speed
    if config.useGPU
        [expected_ersp, ~, ~, ersp_times, freqs] = compute_ersp_stft_fast_gpu(first_group_epochs(:, :, roi_time_idx), config, roi_freqs, eeg_epoch_time(roi_time_idx), fast_bootstrap);
        [unexpected_ersp, ~, ~, ~, ~] = compute_ersp_stft_fast_gpu(second_group_epochs(:, :, roi_time_idx), config, roi_freqs, eeg_epoch_time(roi_time_idx), fast_bootstrap);
                                        % compute_ersp_stft_fast_gpu(epochs, fs, config, times, boot_strap_len)
    else
        [expected_ersp, ~, ~, ersp_times, freqs] = compute_ersp_stft_fast_cpu(first_group_epochs(:, :, roi_time_idx), config, roi_freqs, eeg_epoch_time(roi_time_idx), fast_bootstrap);
        [unexpected_ersp, ~, ~, ~, ~] = compute_ersp_stft_fast_cpu(second_group_epochs(:, :, roi_time_idx), config, roi_freqs, eeg_epoch_time(roi_time_idx), fast_bootstrap);
    end
end