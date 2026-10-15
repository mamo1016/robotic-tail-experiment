function files=all_files(subject_folder)
    % Find all sessions for this subject
    files=struct();
    files.eeg_files = dir(fullfile(subject_folder, 'EEG-st-*.mat'));
    files.cleaned_emg_files = dir(fullfile(subject_folder, 'cleaned_Emg-st-*.mat'));
    files.foot_files = dir(fullfile(subject_folder, 'Foot-st-*.mat'));
    files.tail_emg_files = dir(fullfile(subject_folder, 'Emg-st-*.mat'));
    files.contami_split_emg_files = dir(fullfile(subject_folder, 'contamination_splitcleaned_Emg-st-*.mat'));