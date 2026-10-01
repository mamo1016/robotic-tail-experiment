function [foot_feature_matrix, feature_labels] = generate_foot_features(foot_data_all_subs, config)
    sub_list = fieldnames(foot_data_all_subs);
    n_subs = length(sub_list);

    % Define Windows and Features
    time_windows = {[-0.5, 0], [0, 0.5], [0.5, 1.0], [1.0, 1.5], [0, 1.5]};
    window_names = {'Baseline(-0.5_to_0s)', 'Reflex(0_to_0.5s)', 'EarlyComp(0.5_to_1s)', 'LateComp(1_to_1.5s)', 'TotalPost(0_to_1.5s)'};

    conditions = {'epochs_expected', 'epochs_unexpected', 'epochs_right_before', 'epochs_right_after'};

    % Initialize
    foot_feature_matrix = [];
    feature_labels = {};
    
    % Generate Time Vector based on config
    times = config.epoch_start:1/config.eeg_fs:config.epoch_end;

    for s = 1:n_subs
        sub_name = sub_list{s};
        sub_data = foot_data_all_subs.(sub_name);
        
        all_cond_features = [];

        for c_idx = 1:length(conditions)
            cond_name = conditions{c_idx};
            has_data = isfield(sub_data, cond_name);
            
            if has_data
                % data is [trials x time] for a single physiological channel
                data = sub_data.(cond_name); 
            end

            cond_features = [];
            temp_labels = {};
            chan_label = 'FootCOP'; % Single signal representation
            
            for w_idx = 1:length(time_windows)
                win = time_windows{w_idx};
                win_name = window_names{w_idx};

                if has_data
                    % Extract window time points
                    t_idx = (times >= win(1)) & (times <= win(2));
                    
                    if sum(t_idx) > 0
                        % win_data is [trials x time_window]
                        win_data = data(:, t_idx);

                        % 1. RMS -> Calculate for each trial, then average across trials
                        rms_val = mean(sqrt(mean(win_data.^2, 2)), 1);

                        % 2. Area (Absolute Integral) -> Calculate for each trial, then average
                        area_val = mean(trapz(abs(win_data), 2), 1);
                    else
                        rms_val = 0;
                        area_val = 0;
                    end
                else
                    % Zero pad missing conditions ensuring stable array dims
                    rms_val = 0;
                    area_val = 0;
                end

                cond_features = [cond_features, rms_val, area_val];
                
                if s == 1
                    temp_labels{end+1} = [cond_name, '_', chan_label, '_RMS_', win_name];
                    temp_labels{end+1} = [cond_name, '_', chan_label, '_Area_', win_name];
                end
            end
            
            all_cond_features = [all_cond_features, cond_features];
            if s == 1
                feature_labels = [feature_labels, temp_labels];
            end
        end
        foot_feature_matrix = [foot_feature_matrix; all_cond_features];
    end
end
