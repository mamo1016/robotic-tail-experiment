function tiled = ersp_tile_foot_ex(ersp_matrix, time_range, freq_range)

    
    % Define the boundaries of your grid (ROI)
    % You can adjust these to focus only on 0s-2s and 8Hz-30Hz
    t_start = -2;       t_end = 3.0; 
    f_start = 0;       f_end = 50;
    
    tile_size.f_step = 3;
    tile_size.t_step = 0.2;
    
    f_vec = f_start : tile_size.f_step : (f_end - tile_size.f_step);
    t_vec = t_start : tile_size.t_step : (t_end - tile_size.t_step);

    loop_len = size(ersp_matrix,1);
    tiled_mat = zeros(loop_len, length(f_vec),length(t_vec));
    % tile_labels = cell(length(f_vec),length(t_vec));
    

    for loop = 1:loop_len
        for f_num = 1:length(f_vec)
            f = f_vec(f_num);
            % --- Loop through Time Windows ---
            for t_num = 1:length(t_vec)
                t = t_vec(t_num);
                % 1. Find the Indices (The Mapping Step)
                % "Which rows correspond to 8Hz - 13Hz?"
                f_indices = freq_range >= f & freq_range < (f + tile_size.f_step);
                
                % "Which columns correspond to 0.5s - 0.75s?"
                t_indices = time_range >= t & time_range < (t + tile_size.t_step);
                
                % 2. Extract the Sub-Matrix (The "Tile")
                sub_block = ersp_matrix(loop, f_indices, t_indices);
                
                % 3. Calculate Average Power (Downsampling)
                % If the block is empty (indices not found), handle safely
                if isempty(sub_block)
                    avg_val = 0; 
                else
                    avg_val = mean(sub_block(:)); 
                end
                
                % 4. Append to your Feature Vector
                tiled_mat(loop, f_num, t_num) = avg_val;
                
                % 5. Create a Label (Useful for debugging/plotting later)
                % tile_labels{f_num,t_num} = sprintf('%.2fs_%.0fHz', t, f);
            end
        end
    end
    % visualisation_ersp(tiled_mat, figure, t_vec, f_vec); title("subA unexpected - expected ERSP");
    tiled.matrix = tiled_mat;
    % tiled.label = tile_labels;
    tiled.f = f_vec;
    tiled.t = t_vec;
    
end