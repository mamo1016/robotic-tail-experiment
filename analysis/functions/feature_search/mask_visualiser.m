function [selected_t, selected_f] = mask_visualiser(best_mask_x,f,t)
    % 2. Calculate Grid Dimensions
    num_freqs = length(f);
    num_times = length(t);
    
    % 3. Find the "Winner" Indices
    selected_indices = find(best_mask_x == 1);
    
    % 4. Convert Linear Index -> 2D Grid Coordinates
    % MATLAB fills matrices column by column (Freq first, then Time)
    [f_idx, t_idx] = ind2sub([num_freqs, num_times], selected_indices);
    % 5. Map Indices to Real Values (Seconds and Hz)
    selected_t = t(t_idx);
    selected_f = f(f_idx);
end