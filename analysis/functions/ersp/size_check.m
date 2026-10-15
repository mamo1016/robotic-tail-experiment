function each_size = size_check(fields,interp_epoch)
    each_size = struct();
    size_check = zeros(4,26,3); % trial type, sub num, session num
    for struct_num = 1:length(fields)
            current_name = fields{struct_num};           % Get the name (string)    
            trial_type = fieldnames(interp_epoch.(current_name)); 
            tokens = regexp(current_name, 'subject(\d+)_*session(\d+)', 'tokens'); % get sub num and ses num
        
            for trial_type_num = 1: length(trial_type)
                trial_type_name = trial_type{trial_type_num};
                size_num = [str2double(tokens{1,1}{1}), str2double(tokens{1,1}{2})];
                size_check(trial_type_num, size_num(1), size_num(2)) = size(interp_epoch.(current_name).(trial_type_name).remained_eeg,1); 
            end
    end
    ex=squeeze(size_check(1,:,:));
    uex=squeeze(size_check(2,:,:));
    bex=squeeze(size_check(3,:,:));
    aex=squeeze(size_check(4,:,:));
    
    each_size.expected = ex;
    each_size.unexpected = uex;
    each_size.right_before = bex;
    each_size.right_after = aex;