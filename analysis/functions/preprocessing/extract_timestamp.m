%% Helper function to extract timestamp from filename
function timestamp = extract_timestamp(filename)
    parts = strsplit(filename, '-st-');
    if length(parts) > 1
        timestamp = strtrim(parts{2});
        timestamp = strrep(timestamp, '.mat', '');
    else
        timestamp = '';
    end
end