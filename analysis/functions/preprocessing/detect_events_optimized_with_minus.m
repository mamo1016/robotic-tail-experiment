%% 3. OPTIMIZED EVENT DETECTION
function [trial_data, single_session_dataset] = detect_events_optimized_with_minus(single_session_dataset, eeg_to_emg_offset, config)
    trial_data = struct();
    fprintf('    Detecting events with vectorized approach...\n');
    

    %foot parameter obtain

    single_session_dataset.load_average = [single_session_dataset.foot_raw_dist_L, single_session_dataset.foot_raw_dist_R]*config.foot_corrdination;
    single_session_dataset.load_average = single_session_dataset.load_average(2:end-1,:);
    single_session_dataset.foot_time_dist_Both = (single_session_dataset.foot_time_dist_L+single_session_dataset.foot_time_dist_R)/2;
    single_session_dataset.foot_time_dist_Both = single_session_dataset.foot_time_dist_Both(2:end-1,:);

    % Vectorized event detection
    tail_diff = diff(single_session_dataset.tail_data);
    transitions = find(tail_diff ~= 0);

    
    tail_diff_contami_split = diff(single_session_dataset.tail_data_contami_split);
    transitions_contami_split = find(tail_diff_contami_split >= 2);

    % here check indices
    % Find expected (0->1) and unexpected (1->-1) transitions
    expected_indices_tail_emg = transitions(single_session_dataset.tail_data(transitions + 1) == 1 & transitions > config.restDur) + 1;
    unexpected_indices_tail_emg = transitions(single_session_dataset.tail_data(transitions + 1) == -1 & transitions > config.restDur) + 1;
    trial_data.expected_indices_tail_emg = expected_indices_tail_emg;
    trial_data.unexpected_indices_tail_emg = unexpected_indices_tail_emg;

    % Find expected (0->2) and unexpected (0->3) transitions
    bef_indices_contami_split = transitions_contami_split(single_session_dataset.tail_data_contami_split(transitions_contami_split + 1) == 2 & transitions_contami_split > config.restDur) + 1;
    aft_indices_contami_split = transitions_contami_split(single_session_dataset.tail_data_contami_split(transitions_contami_split + 1) == 3 & transitions_contami_split > config.restDur) + 1;
    
    bef_indices_contami_split_2 = find(tail_diff_contami_split == 2);    
    aft_indices_contami_split_2 = find(tail_diff_contami_split == 3); 

    % get minus one epoch samples (tail time domain)
    bef_indices_contami_split_minus = zeros(size(bef_indices_contami_split_2));
    for i=1:length(bef_indices_contami_split_2)
        temp=expected_indices_tail_emg-bef_indices_contami_split_2(i); 
        [~,idx]=min(abs(temp));
        bef_indices_contami_split_minus(i)=expected_indices_tail_emg(idx-1)-2;
    end

    closest_finder=zeros(length(expected_indices_tail_emg),length(bef_indices_contami_split));
    % bef_indices_contami_split=zeros(size(bef_indices_contami_split));
    for i=1:length(bef_indices_contami_split)
        closest_finder(:,i)=expected_indices_tail_emg-bef_indices_contami_split(i);
        [~,idx]=min(abs(closest_finder(:,i)));
        if idx>2
            bef_indices_contami_split(i)=expected_indices_tail_emg(idx-1);
        end
    end

    % Combine and sort events
    all_event_indices = [expected_indices_tail_emg; unexpected_indices_tail_emg];
    event_labels = [ones(size(expected_indices_tail_emg)); -ones(size(unexpected_indices_tail_emg))];

    all_event_indices_contami = [bef_indices_contami_split; aft_indices_contami_split];
    all_event_indices_contami_2 = [bef_indices_contami_split_2; aft_indices_contami_split_2];
    
    event_labels_contami = [ones(size(bef_indices_contami_split)); -ones(size(aft_indices_contami_split))];
    event_labels_contami_2 = [ones(size(bef_indices_contami_split_2)); -ones(size(aft_indices_contami_split_2))];
    
    
    % Sort by time
    [all_event_indices, sort_idx] = sort(all_event_indices);
    event_labels = event_labels(sort_idx);

    [all_event_indices_contami, sort_idx] = sort(all_event_indices_contami);
    event_labels_contami = event_labels_contami(sort_idx);

    [all_event_indices_contami_2, sort_idx] = sort(all_event_indices_contami_2);
    event_labels_contami_2 = event_labels_contami_2(sort_idx);

    % Convert to EEG time domain efficiently
    adjusted_tail_emg_time =  single_session_dataset.tail_emg_time + eeg_to_emg_offset;
    emg_event_times =         adjusted_tail_emg_time(all_event_indices);    

    expected_emg_event_times = adjusted_tail_emg_time(expected_indices_tail_emg);    
    unexpected_emg_event_times = adjusted_tail_emg_time(unexpected_indices_tail_emg);    

    trial_data.expected_emg_event_times = expected_emg_event_times;
    trial_data.unexpected_emg_event_times = unexpected_emg_event_times;

    %old code
    % emg_event_times_contami = adjusted_tail_emg_time(all_event_indices_contami); 
    
    %new code
    emg_event_times_contami = adjusted_tail_emg_time(all_event_indices_contami_2);
    emg_event_times_contami_minus = adjusted_tail_emg_time(bef_indices_contami_split_minus);

    % Vectorized time matching
    [~, trial_data.eeg_event_indices] = arrayfun(@(x) min(abs(single_session_dataset.eeg_time - x)), emg_event_times);
    [~, trial_data.eeg_event_indices_contami] = arrayfun(@(x) min(abs(single_session_dataset.eeg_time - x)), emg_event_times_contami);
    [~, trial_data.eeg_event_indices_contami_minus] = arrayfun(@(x) min(abs(single_session_dataset.eeg_time - x)), emg_event_times_contami_minus);

    [~, trial_data.expected_indices_eeg] = arrayfun(@(x) min(abs(single_session_dataset.eeg_time - x)), adjusted_tail_emg_time(expected_indices_tail_emg));
    [~, trial_data.unexpected_indices_eeg] = arrayfun(@(x) min(abs(single_session_dataset.eeg_time - x)), adjusted_tail_emg_time(unexpected_indices_tail_emg));


    [~, trial_data.expected_indices_foot] = arrayfun(@(x) min(abs(single_session_dataset.foot_time_dist_Both - x)), adjusted_tail_emg_time(expected_indices_tail_emg));
    [~, trial_data.unexpected_indices_foot] = arrayfun(@(x) min(abs(single_session_dataset.foot_time_dist_Both - x)), adjusted_tail_emg_time(unexpected_indices_tail_emg));

    % Separate expected and unexpected trials
    trial_data.expected_trials = find(event_labels == 1);
    trial_data.unexpected_trials = find(event_labels == -1);



    % before and after
    trial_data.bef_trials = find(event_labels_contami_2 == 1); %think how to get minus before trials 10/12
    trial_data.aft_trials = find(event_labels_contami_2 == -1);


    % split all expected index between each unexpected trial
    first_idx_temp=1;
    second_idx_temp=trial_data.unexpected_trials(1)-1;
    expected_idx=zeros(2,length(trial_data.unexpected_trials));
    for i = 1:length(trial_data.unexpected_trials)
        expected_idx(1,i)=first_idx_temp;
        expected_idx(2,i)=second_idx_temp;
        if i ~= length(trial_data.unexpected_trials)
            first_idx_temp=trial_data.unexpected_trials(i)+1;
            second_idx_temp=trial_data.unexpected_trials(i+1)-1;
        end
    end
    
    random_expected=zeros(size(expected_idx));
    for i = 1:size(expected_idx,2)
        random_expected(1,i) = randi([expected_idx(1,i), expected_idx(2,i)]);
        random_expected(2,i) = random_expected(1,i);
        while random_expected(1,i) == random_expected(2,i)
            random_expected(2,i) = randi([expected_idx(1,i), expected_idx(2,i)]);
        end
    end    
    % Sort each column in ascending order
    sorted_random_expected = sort(random_expected);

    % random expected 1 and random expected 2
    trial_data.expected_first_trials = sorted_random_expected(1,:); 
    trial_data.expected_second_trials = sorted_random_expected(2,:);

    % eeg_event_indices_minus=zeros(size(eeg_event_indices_contami));
    % for i=1:length(eeg_event_indices_contami)
    %     eeg_event_indices_minus(i)=eeg_event_indices((eeg_event_indices==eeg_event_indices_contami(i))-1);
    % end

    % bef_minus_indices_contami_split=

    fprintf('Found %d expected and %d unexpected events\n', length(trial_data.expected_trials), length(trial_data.unexpected_trials));

    fprintf(config.log_file, '    Detected %d expected trials and %d unexpected trials\n', ...
            length(trial_data.expected_trials), length(trial_data.unexpected_trials));
        
    % Verify we have enough trials to proceed
    if length(trial_data.expected_trials) < 3 || length(trial_data.unexpected_trials) < 3 || length(trial_data.bef_trials) < 1 || length(trial_data.aft_trials) < 1
        fprintf('  Warning: Not enough trials detected. Expected: %d, Unexpected: %d\n', ...
            length(trial_data.expected_trials), length(trial_data.unexpected_trials));
        fprintf(config.log_file, '  Warning: Not enough trials detected. Expected: %d, Unexpected: %d\n', ...
            length(trial_data.expected_trials), length(trial_data.unexpected_trials));
        % continue; % Skip this session
        trial_data.skip = 1;
    else
        trial_data.skip = 0;
    end

    fprintf('  Detecting artifacts...\n');
    trial_data.artifact_mask = any(abs(single_session_dataset.eeg_filtered) > config.artifact_threshold, 2);         % Vectorized artifact detection is much cleaner and faster        
    artifact_percentage = sum(trial_data.artifact_mask) / length(trial_data.artifact_mask) * 100; % Store artifact information
    fprintf(config.log_file, '    Artifacts detected in %.2f%% of samples\n', artifact_percentage);

    % Define epoch and baseline samples (these lines are still needed)
    trial_data.epoch_samples = round([config.epoch_start * config.eeg_fs, config.epoch_end * config.eeg_fs]);
    baseline_samples = round([config.baseline_start * config.eeg_fs, config.baseline_end * config.eeg_fs]);
    trial_data.baseline_idx = (-trial_data.epoch_samples(1) + baseline_samples(1) + 1):(-trial_data.epoch_samples(1) + baseline_samples(2) + 1);


end