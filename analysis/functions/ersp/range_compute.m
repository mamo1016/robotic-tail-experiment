function [time_range, freq_range] = range_compute(config)
    % 1. Define your parameters (from your config)
    L = config.window_size;
    noverlap = config.overlap;
    fs = config.eeg_fs;
    
    N = fs*(config.epoch_end - config.epoch_start)+1;;
    
    % 2. Calculate the 'Hop Size' (Step size between windows)
    hop_size = L - noverlap;
    
    % 3. Calculate the number of segments (k)
    % MATLAB drops the last segment if it doesn't fit a full window size
    k = fix((N - noverlap) / hop_size);
    
    % 4. generate the indices for the start of each window (0-based for math)
    window_indices = 0 : (k - 1);
    
    % 5. Calculate t_spec
    % Formula: (Start_Index * Hop) + (Half_Window_Duration)
    % We divide by fs to convert samples to seconds.
    t_spec = (window_indices * hop_size + (L - 1) / 2) / fs;
    time_range = t_spec + config.epoch_start; 
    
    nfft = config.nfft;
    
    % Determine if nfft is even or odd to set the correct range
    if mod(nfft, 2) == 0
        % Case 1: Even nfft (Most common, e.g., 256, 512, 1024)
        % Range: 0 to fs/2
        % Number of points: (nfft / 2) + 1
        freq_range = (0 : nfft/2)' * fs / nfft;
    else
        % Case 2: Odd nfft
        % Range: 0 to slightly less than fs/2
        % Number of points: (nfft + 1) / 2
        freq_range = (0 : (nfft-1)/2)' * fs / nfft;
    end