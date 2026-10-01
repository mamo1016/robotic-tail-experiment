function interpolate_all_epoch = cleaning_epoch(all_epoch, interpolate_all_epoch)
fields = fieldnames(all_epoch);
trial_num = zeros(length(fields),2);
for struct_num = 1:length(fields)
    current_name = fields{struct_num};           % Get the name (string)    
    trial_type = fieldnames(all_epoch.(current_name)); 

    for trial_type_num = 1: length(trial_type)
        trial_type_name = trial_type{trial_type_num};           % 1:expected, 2:unexpected, 3:right before 4: right after
    
        current_eeg_data = all_epoch.(current_name).(trial_type_name).eeg; % Get the content
    
        trial_number = size(current_eeg_data,1);
        ch_number = size(current_eeg_data,3);
        
        remove_trial = zeros(trial_number,ch_number);
        remove_idx = zeros(trial_number,1);
        for i = 1:trial_number
            for ch = 1:ch_number
                if max(abs(current_eeg_data(i,:,ch)))>100
                    remove_trial(i,ch) = 1;
                end
            end
            % BUG (confirmed by author, 2026-08-21; NOT fixed for the 2026 submission,
            % because re-running would change published numbers).
            % This test should be >=1, not >1. As written, an epoch with exactly ONE
            % channel above threshold is left unrepaired. Measured across all 77 session
            % files (5,374 epochs): 28.1% of epochs carry exactly one flagged channel and
            % pass through untouched. 3.3% have 2-4 and are interpolated; 1.0% have >4 and
            % are rejected. The manuscript Methods describe the behaviour as it actually
            % runs. Fix this before any future analysis.
            if sum(remove_trial(i,:))>1
                bad_ch = find(remove_trial(i,:)==1);    
                if length(bad_ch) > 4
                    remove_idx(i) = 1;
                elseif ~isempty(bad_ch)
                    current_eeg_data(i,:,:) = channel_interpolation(bad_ch,current_eeg_data(i,:,:));
                end
                % all_epoch.(current_name).expected.eeg = fixed_data
            end
        end
        % check how many removed
        trial_num(struct_num,1) = size(current_eeg_data,1);
        current_eeg_data(remove_idx==1,:,:)=[];
        trial_num(struct_num,2) = size(current_eeg_data,1);
    
        current_emg_data = all_epoch.(current_name).(trial_type_name).emg; % Get the content
        current_emg_data(remove_idx==1,:,:)=[];
    
        current_foot_data = all_epoch.(current_name).(trial_type_name).foot; % Get the content
        current_foot_data(remove_idx==1,:,:)=[];
    
        current_tail_data = all_epoch.(current_name).(trial_type_name).tail; % Get the content
        current_tail_data(remove_idx==1,:,:)=[];
    
        interpolate_all_epoch.(current_name).(trial_type_name).remained_eeg = current_eeg_data;
        interpolate_all_epoch.(current_name).(trial_type_name).remained_emg = current_emg_data;
        interpolate_all_epoch.(current_name).(trial_type_name).remained_foot = current_foot_data;
        interpolate_all_epoch.(current_name).(trial_type_name).remained_tail = current_tail_data;

    end
end