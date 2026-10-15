function samples = subtract_idx(trial_size, min_trial_num)

    % 1. Create your original array
    data = 1:1:trial_size-1;
    
    % 2. Generate 100 equally spaced points between index 1 and 150
    % linspace(x1, x2, n) generates n points between x1 and x2
    ideal_indices = linspace(1, trial_size-1, min_trial_num-1);
    
    % 3. Round them to the nearest integers to get valid array indices
    valid_indices = round(ideal_indices);
    
    % 4. Extract the samples using these indices
    samples = data(valid_indices);
end