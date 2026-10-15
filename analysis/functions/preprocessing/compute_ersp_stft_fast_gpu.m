%% 8. SIMPLIFIED FAST GPU STFT (Missing Function)
function [ersp, itc, powbase, freqs] = compute_ersp_stft_fast_gpu(epochs, config, freq_range, boot_strap_len)
    
    % % Optimized parameters for speed
    % window_length = 512;
    % overlap = 0.90;
    % nfft = 1024;
    hop = round(config.window_length * (1 - config.overlap));
    config.eeg_epoch_time = linspace(config.epoch_start, config.epoch_end, n_samples);

    [n_epochs, n_channels, n_eeg_epoch_time] = size(epochs);
    
    % Find frequency indices
    freq_indices = round(freq_range * config.nfft / config.eeg_fs) + 1;
    freq_indices = freq_indices(freq_indices <= config.nfft/2);
    n_freqs = length(freq_indices);
    freqs = (freq_indices - 1) * config.eeg_fs / config.nfft;
    
    % Calculate time points
    n_windows = floor((n_eeg_epoch_time - config.window_length) / hop) + 1;
    times_stft = config.eeg_epoch_time(1) + (0:n_windows-1) * hop / config.eeg_fs;
    
    % Initialize results
    ersp = zeros(n_channels, n_freqs, n_windows, 'single');
    itc = zeros(n_channels, n_freqs, n_windows, 'single');
    powbase = zeros(n_channels, n_freqs, 'single');
    
    n_bootstraps = min(boot_strap_len, 500);
    window_func = hamming(config.window_length);

    % Process each channel
    for ch = 1:n_channels
        channel_data = squeeze(epochs(:, ch, :));
        
        % Pre-allocate
        all_tf = zeros(n_epochs, n_freqs, n_windows, 'single');
        
        % Compute STFT for all trials on GPU
        for e = 1:n_epochs
            signal_gpu = gpuArray(channel_data(e, :));
            [~, ~, ~, S] = spectrogram(signal_gpu, window_func, config.window_length-hop, config.nfft, config.eeg_fs);
            all_tf(e, :, :) = gather(abs(S(freq_indices, :)).^2);
        end
        
        % Fast bootstrapping
        bootstrap_indices = randi(n_epochs, n_bootstraps, boot_strap_len);
        boot_power = zeros(n_bootstraps, n_freqs, n_windows, 'single');
        
        for b = 1:n_bootstraps
            boot_power(b, :, :) = squeeze(mean(all_tf(bootstrap_indices(b, :), :, :), 1));
        end
        
        % Calculate baseline and ERSP
        baseline_windows = times_stft < 0;
        for f = 1:n_freqs
            baseline_power = mean(boot_power(:, f, baseline_windows), [1, 3]);
            powbase(ch, f) = baseline_power;
            
            mean_power = squeeze(mean(boot_power(:, f, :), 1));
            ersp(ch, f, :) = 10 * log10(mean_power / max(baseline_power, eps));
        end
    end
    
    config.eeg_epoch_time = times_stft;
end