
%% Helper function for foot data computation - Optimized version
function datasets = footCompute_optimized(datasets)
    left = datasets.foot_raw(:,1);
    right = datasets.foot_raw(:,2);
    foot_time = datasets.foot_time;

    % Find positions of 255
    left_idx = find(left == 255);
    right_idx = find(right == 255);
    
    % Ensure we have valid data
    if isempty(left_idx) || isempty(right_idx)
        datasets.leftSums = [];
        datasets.rightSums = [];
        datasets.foot_time_left = [];
        return;
    end
    
    % Use foot time at marker points
    datasets.foot_time_left = foot_time(left_idx);
    
    % Process left foot data using vectorized approach
    datasets.leftSums = zeros(length(left_idx), 1, 'single');
    for i = 1:length(left_idx)-1
        startIdx = left_idx(i) + 1;
        endIdx = left_idx(i+1) - 1;
        if startIdx <= endIdx
            datasets.leftSums(i) = sum(left(startIdx:endIdx));
        end
    end
    
    % Process the last segment
    if ~isempty(left_idx) && left_idx(end) < length(left)
        startIdx = left_idx(end) + 1;
        datasets.leftSums(end) = sum(left(startIdx:end));
    end
    
    % Process right foot data using the same approach
    datasets.rightSums = zeros(length(right_idx), 1, 'single');
    for i = 1:length(right_idx)-1
        startIdx = right_idx(i) + 1;
        endIdx = right_idx(i+1) - 1;
        if startIdx <= endIdx
            datasets.rightSums(i) = sum(right(startIdx:endIdx));
        end
    end
    
    % Process the last segment
    if ~isempty(right_idx) && right_idx(end) < length(right)
        startIdx = right_idx(end) + 1;
        datasets.rightSums(end) = sum(right(startIdx:end));
    end

    if ~isempty(datasets.leftSums) && ~isempty(datasets.rightSums) && ~isempty(datasets.foot_time_left)
        left_foot_len = min(length(datasets.leftSums), length(datasets.foot_time_left));
        datasets.leftSums = datasets.leftSums(1:left_foot_len);
        datasets.rightSums = [datasets.rightSums; 0];
        datasets.rightSums = datasets.rightSums(1:left_foot_len);
        datasets.foot_center = datasets.leftSums - datasets.rightSums;
        datasets.foot_time_left = datasets.foot_time_left(1:left_foot_len);
    else
        datasets.foot_center = [];
        fprintf('  Warning: Could not process foot data\n');
    end

end