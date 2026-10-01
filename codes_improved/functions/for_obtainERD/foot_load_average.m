%% 3. OPTIMIZED EVENT DETECTION
function [single_session_dataset] = foot_load_average(single_session_dataset, config)
    fprintf('foot load average...\n');
    
    %foot parameter obtain
    load_average = [single_session_dataset.foot_raw_dist_L, single_session_dataset.foot_raw_dist_R]*config.foot_corrdination;
    single_session_dataset.load_average = load_average(2:end-1,:);
    single_session_dataset.foot_time_dist_Both = (single_session_dataset.foot_time_dist_L+single_session_dataset.foot_time_dist_R)/2;
    single_session_dataset.foot_time_dist_Both = single_session_dataset.foot_time_dist_Both(2:end-1,:);
end