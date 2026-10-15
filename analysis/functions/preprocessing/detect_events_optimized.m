%% 3. OPTIMIZED EVENT DETECTION
function [expected_trials, unexpected_trials, eeg_event_indices, bef_trials, aft_trials, eeg_event_indices_contami] = detect_events_optimized(tail_data, tail_data_contami_split, tail_emg_time, eeg_time, eeg_to_emg_offset)
    
    fprintf('    Detecting events with vectorized approach...\n');
    
    restDur = 250*20;
    
    % Vectorized event detection
    tail_diff = diff(tail_data);
    transitions = find(tail_diff ~= 0);

    
    tail_diff_contami_split = diff(tail_data_contami_split);
    transitions_contami_split = find(tail_diff_contami_split >= 2);
    
    % Find expected (0->1) and unexpected (1->-1) transitions
    expected_indices = transitions(tail_data(transitions + 1) == 1 & transitions > restDur) + 1;
    unexpected_indices = transitions(tail_data(transitions + 1) == -1 & transitions > restDur) + 1;
    
    % Find expected (0->2) and unexpected (0->3) transitions
    bef_indices_contami_split = transitions_contami_split(tail_data_contami_split(transitions_contami_split + 1) == 2 & transitions_contami_split > restDur) + 1;
    aft_indices_contami_split = transitions_contami_split(tail_data_contami_split(transitions_contami_split + 1) == 3 & transitions_contami_split > restDur) + 1;
    
    % closest_finder=zeros(length(expected_indices),length(bef_indices_contami_split));
    % eeg_event_indices_contami_minus=zeros(size(bef_indices_contami_split));
    % for i=1:length(bef_indices_contami_split)
    %     closest_finder(:,i)=expected_indices-bef_indices_contami_split(i);
    %     [~,idx]=min(abs(closest_finder(:,i)));
    %     if idx>2
    %         eeg_event_indices_contami_minus(i)=expected_indices(idx-1);
    %     end
    % end

    % Combine and sort events
    all_event_indices = [expected_indices; unexpected_indices];
    event_labels = [ones(size(expected_indices)); -ones(size(unexpected_indices))];

    all_event_indices_contami = [bef_indices_contami_split; aft_indices_contami_split];
    event_labels_contami = [ones(size(bef_indices_contami_split)); -ones(size(aft_indices_contami_split))];
    
    
    % Sort by time
    [all_event_indices, sort_idx] = sort(all_event_indices);
    event_labels = event_labels(sort_idx);

    [all_event_indices_contami, sort_idx] = sort(all_event_indices_contami);
    event_labels_contami = event_labels_contami(sort_idx);


    % Convert to EEG time domain efficiently
    adjusted_tail_emg_time =  tail_emg_time + eeg_to_emg_offset;
    emg_event_times =         adjusted_tail_emg_time(all_event_indices);    
    emg_event_times_contami = adjusted_tail_emg_time(all_event_indices_contami);


    % Vectorized time matching
    [~, eeg_event_indices] = arrayfun(@(x) min(abs(eeg_time - x)), emg_event_times);
    [~, eeg_event_indices_contami] = arrayfun(@(x) min(abs(eeg_time - x)), emg_event_times_contami);
    
    % Separate expected and unexpected trials
    expected_trials = find(event_labels == 1);
    unexpected_trials = find(event_labels == -1);

    bef_trials = find(event_labels_contami == 1);
    aft_trials = find(event_labels_contami == -1);

    fprintf('    Found %d expected and %d unexpected events\n', length(expected_trials), length(unexpected_trials));
end