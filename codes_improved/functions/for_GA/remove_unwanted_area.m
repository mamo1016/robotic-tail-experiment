function clean_data = remove_unwanted_area(all_sub_diffs,data, freq_range, f_range)

names = fieldnames(all_sub_diffs);
[n_subs, n_features] = size(data);
clean_data = zeros(n_subs, n_features);

rmv_idx = freq_range<f_range(1) | freq_range>f_range(2);

% Handle multiple blocks of 400 features
num_blocks = n_features / 400;
if mod(n_features, 400) ~= 0
    error('Feature count must be a multiple of 400');
end

for i = 1:n_subs
    full_row = data(i, :);
    cleaned_row = [];
    
    for b = 1:num_blocks
        % Extract block
        idx_start = (b-1)*400 + 1;
        idx_end = b*400;
        block = full_row(idx_start:idx_end);
        
        % Reshape and Mask
        % DYNAMICALLY calculate dimensions based on block length and freq_range length
        n_total = length(block);
        n_freqs = length(freq_range);
        n_times = n_total / n_freqs;
        
        if mod(n_total, n_freqs) ~= 0
           error('Block size (%d) is not divisible by frequency range length (%d)', n_total, n_freqs);
        end
        
        block_mat = reshape(block, n_freqs, n_times);
        block_mat(rmv_idx, :) = 0;
        
        % Flatten and Append
        cleaned_row = [cleaned_row, block_mat(:)'];
    end
    
    clean_data(i, :) = cleaned_row;
end