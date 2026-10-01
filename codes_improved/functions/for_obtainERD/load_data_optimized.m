%% 1. OPTIMIZED FILE LOADING (BIGGEST BOTTLENECK)
% Replace the inefficient matfile approach with direct loading
function datasets = load_data_optimized(subject_folder, matchedfiles, config)
    
    datasets = struct();
    % fprintf('    Loading EEG file...\n');
    % Load EEG data directly - much faster than matfile
    eeg_struct = load(fullfile(subject_folder, matchedfiles.matching_eeg), 'ans');
    datasets.eeg_raw = single(eeg_struct.ans.Data(:, config.coi)); % Extract only channels of interest
    datasets.eeg_time = eeg_struct.ans.Time;
    clear eeg_struct; % Free memory immediately
    
    % fprintf('    Loading EMG files...\n');
    % Load EMG data
    emg_struct = load(fullfile(subject_folder, matchedfiles.matching_tail), 'TailsEMG');
    datasets.emg_raw = single(emg_struct.TailsEMG.Data(:, 3));
    datasets.tail_emg_time = emg_struct.TailsEMG.Time;
    clear emg_struct;
    
    % Load tail output data
    tail_struct = load(fullfile(subject_folder, matchedfiles.matching_emg), 'data');
    datasets.tail_data = tail_struct.data;
    clear tail_struct;

    tail_struct = load(fullfile(subject_folder, matchedfiles.matching_emg_contami_split), 'tail_move_data_contami_split');
    datasets.tail_data_contami_split = tail_struct.tail_move_data_contami_split;
    clear tail_struct;

    % fprintf('    Loading foot data...\n');
    % Load foot data
    foot_struct = load(fullfile(subject_folder, matchedfiles.matching_foot), 'Foot');
    datasets.foot_raw = single(foot_struct.Foot.Data);
    datasets.foot_time = foot_struct.Foot.Time;
    clear foot_struct;
end