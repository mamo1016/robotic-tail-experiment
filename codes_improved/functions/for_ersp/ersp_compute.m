function ersp_compute(fields,interp_epoch, config,out_dir)
    % compute ersp and save
    for struct_num = 1:length(fields)
        ersp_epoch = struct();
        current_name = fields{struct_num};           % Get the name (string)    
        trial_type = fieldnames(interp_epoch.(current_name)); 
    
        for trial_type_num = 1: length(trial_type)
            trial_type_name = trial_type{trial_type_num};           % 1:expected, 2:unexpected, 3:right before 4: right after
    
            current_eeg_data = interp_epoch.(current_name).(trial_type_name).remained_eeg; % Get the content
            % ersp_epoch.(current_name).(trial_type_name).data_ersp = ersp(config, current_eeg_data);  
            ersp_epoch.(current_name).(trial_type_name).data_ersp = ersp(config, current_eeg_data);  

            % interp_epoch.(current_name).(trial_type_name).ersp_parameter = ersp_parameter;       
    
        end
        filename = string(out_dir+"\ersp_ave\ersp_epoch_ave_sub_"+struct_num+".mat");        
        save(filename, 'ersp_epoch', '-v7.3');
        clear ersp_epoch
    end