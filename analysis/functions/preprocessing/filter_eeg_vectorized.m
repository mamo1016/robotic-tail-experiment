%% 2. VECTORIZED FILTERING (MUCH FASTER)
function eeg_filtered = filter_eeg_vectorized(eeg_raw, filter_para, config)
    
    fprintf('    Applying vectorized filtering...\n');
    
    if config.useGPU && gpuDeviceCount > 0
        % GPU-based filtering - process all channels at once
        eeg_gpu = gpuArray(double(eeg_raw));
        
        % Apply all filters in sequence to entire matrix
        eeg_gpu = filtfilt(filter_para.b_high, filter_para.a_high, eeg_gpu);
        eeg_gpu = filtfilt(filter_para.b_low, filter_para.a_low, eeg_gpu);
        eeg_gpu = filtfilt(filter_para.b_notch, filter_para.a_notch, eeg_gpu);
        
        eeg_filtered = gather(single(eeg_gpu));
        clear eeg_gpu;
    else
        % CPU vectorized filtering - much faster than parfor loop
        eeg_double = double(eeg_raw);
        
        % Apply filters to all channels simultaneously
        eeg_double = filtfilt(filter_para.b_high, filter_para.a_high, eeg_double);
        eeg_double = filtfilt(filter_para.b_low, filter_para.a_low, eeg_double);
        eeg_filtered = single(filtfilt(filter_para.b_notch, filter_para.a_notch, eeg_double));
        
        clear eeg_double;
    end

   
end