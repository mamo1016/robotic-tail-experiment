%% Helper function to find a matching file based on timestamp - Optimized
function matching_file = find_matching_file(file_list, reference_timestamp)
    % Pre-allocate memory for time differences
    n_files = length(file_list);
    if n_files == 0
        matching_file = '';
        return;
    end
    
    timeDiff = zeros(1, n_files);
    
    % Define the format for datetime parsing
    formatSpec = 'dd-MMM-yyyy_HH-mm-ss-SSS';
    
    % Convert reference timestamp to datetime only once
    try
        dt_ref = datetime(reference_timestamp, 'InputFormat', formatSpec);
    catch
        matching_file = '';
        return;
    end
    
    % Process all files with vectorization where possible
    for i = 1:n_files
        file_timestamp = extract_timestamp(file_list(i).name);
        try
            dt_file = datetime(file_timestamp, 'InputFormat', formatSpec);
            timeDiff(i) = abs(seconds(dt_file - dt_ref));
        catch
            timeDiff(i) = Inf; % Mark as invalid
        end
    end
    
    % Find minimum time difference
    [~, idx] = min(timeDiff);
    matching_file = file_list(idx).name;
end