function epoch_times_eeg_st = epoch_index(single_session_dataset, eeg_to_emg_offset, config, phase)
    % Vectorized event detection
    tail_diff = diff(single_session_dataset.tail_data_contami_split);
    transitions= find(tail_diff ~= 0);

    % Find expected (0->1) and unexpected (1->-1) transitions     % Find expected (0->2) and unexpected (0->3) transitions
    if phase == 1
        % expected
        epoch_indices_tail_emg = transitions(single_session_dataset.tail_data_contami_split(transitions + 1) == 1 & transitions > config.restDur) + 1;
    elseif phase == 2
        % unexpected
        epoch_indices_tail_emg = transitions(single_session_dataset.tail_data_contami_split(transitions + 1) == -1 & transitions > config.restDur) + 1;
    elseif phase == 3
        % right before
        epoch_indices_tail_emg = transitions(single_session_dataset.tail_data_contami_split(transitions + 1) == 2 & transitions > config.restDur) + 1;
    elseif phase == 4
        % right after
        epoch_indices_tail_emg = transitions(single_session_dataset.tail_data_contami_split(transitions + 1) == 3 & transitions > config.restDur) + 1;
    end

    adjusted_tail_emg_time =  single_session_dataset.tail_emg_time + eeg_to_emg_offset;
    [~, indices_eeg] = arrayfun(@(x) min(abs(single_session_dataset.eeg_time - x)), adjusted_tail_emg_time(epoch_indices_tail_emg));
    epoch_times_eeg_st = single_session_dataset.eeg_time(indices_eeg)-eeg_to_emg_offset;