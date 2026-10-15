function all_epoch = all_epoch_organise(config)
% List of subjects to process
epoch_list = dir(fullfile(config.epochDir, 'Sub*'));

% Create a log file
fprintf(config.log_file, 'EEG Processing Log - %s\n\n', datetime("now"));

all_epoch = struct();
file_num = length(epoch_list);
% fig_epoch = figure('Color', 'w'); % 'Name' sets the window title

for num = 1:file_num
    epoch_file = fullfile(config.epochDir, epoch_list(num).name);
    epoch = load(epoch_file).epoch_data;

    filename = epoch_list(num).name;
    tokens = regexp(filename, 'subject_(\d+)___*session(\d+)', 'tokens');
    
    if ~isempty(tokens)
        sub_id_str = pad(string(tokens{1}{1}), 2, 'left', '0');
        sess_id_str = pad(string(tokens{1}{2}), 2, 'left', '0');

        field_name = string("subject"+sub_id_str+"_session"+sess_id_str); 
        
        all_epoch.(field_name).expected = epoch.epochs_expected;
        all_epoch.(field_name).unexpected = epoch.epochs_unexpected;
        all_epoch.(field_name).right_before = epoch.epochs_right_before;
        all_epoch.(field_name).right_after = epoch.epochs_right_after;
    end
end
all_epoch = orderfields(all_epoch);
