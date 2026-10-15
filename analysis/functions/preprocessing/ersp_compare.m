function ersps = ersp_compare(trial_data,config)            
    % Define needed variables
    % n_samples = size(trial_data.expected_epochs, 3);
    
    
    boot_strap_len = min(size(trial_data.expected_epochs, 1), size(trial_data.unexpected_epochs, 1));
    boot_strap_len_contami = min(size(trial_data.bef_epochs, 1), size(trial_data.aft_epochs, 1));
    
    % Use the much faster, optimized time-frequency function
    ersps = compute_ersp_fast(trial_data.expected_epochs, trial_data.unexpected_epochs, config, boot_strap_len);
    
    [ersps.bef_ersp, ersps.aft_ersp, ~, ~] = ...
        compute_ersp_fast(trial_data.bef_epochs, trial_data.aft_epochs, config, eeg_epoch_time, boot_strap_len_contami);
    
    [~, ersps.minus_ersp, ~, ~] = ...
        compute_ersp_fast(trial_data.bef_epochs, trial_data.minus_bef_epochs, config, eeg_epoch_time, boot_strap_len_contami);
end