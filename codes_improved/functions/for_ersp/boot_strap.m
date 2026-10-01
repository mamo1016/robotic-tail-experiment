function [final_ersp, final_std] = boot_strap(config, ersp_boot, boot_strap_pick_num, select_idx)
    % ersp_data_boot = ersp_epoch.subject01_session01.expected.data_ersp;
    rng('shuffle');

    %% 1. Setup Parameters
    % 'ersp_boot' is your variable: [Trials x Freq x Time]
    % Example: 50 x 257 x 281
    % [n_trials, ~, ~] = size(ersp_boot);
    % select_idx = randperm(n_trials, boot_strap_pick_num);
    ersp_boot = ersp_boot(select_idx,:,:);
    [n_trials, n_freqs, n_times] = size(ersp_boot);

    % Pre-allocate memory to store the sum of all bootstraps
    % (Accumulating the sum is faster/lighter than storing 1000 separate matrices)
    sum_ersp = zeros(n_freqs, n_times);
    sum_sq_ersp = zeros(n_freqs, n_times); % Optional: to calculate SD later
    %% 2. The Bootstrap Loop
    for i = 1:config.n_boots
        % A. Randomly select indices
        % randperm(N, K) returns K unique integers from 1 to N (sampling WITHOUT replacement)
        % This simulates "what if this subject only had these 12 trials?"
        rng('shuffle');
        idx = randi([1 n_trials], 1, boot_strap_pick_num);
        
        % B. Extract the subset of data
        subset_data = ersp_boot(idx, :, :);
        
        % C. Calculate the mean for this specific iteration
        % mean(..., 1) collapses the trial dimension
        iteration_avg = squeeze(mean(subset_data, 1)); 
        
        % D. Accumulate (Add to running total)
        sum_ersp = sum_ersp + iteration_avg;
        
        % Optional: Accumulate square for Standard Deviation calculation
        sum_sq_ersp = sum_sq_ersp + (iteration_avg .^ 2);
    end
    
    %% 3. Final Calculation
    % The "Grand Average" of your bootstraps
    final_ersp = sum_ersp / config.n_boots;
    
    % Optional: Calculate Bootstrap Standard Deviation (Estimate of Variability)
    % Formula: sqrt( E[x^2] - (E[x])^2 )
    final_std = sqrt( (sum_sq_ersp / config.n_boots) - (final_ersp .^ 2) );
    
    % fprintf('Done. Final ERSP size: %s\n', mat2str(size(final_ersp)));
end