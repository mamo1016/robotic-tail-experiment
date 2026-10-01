function boot_strap_pick_num = min_epoch(interp_epoch)
    each_size = size_check(fieldnames(interp_epoch),interp_epoch);
    fields = fieldnames(each_size);
    mini_val = zeros(length(fields),1);
    for i = 1:length(fields)
        current_name = fields{i};           % Get the name (string)    
        size_list = each_size.(current_name); 
        size_sum = sum(size_list,2);
        mini_val(i) = min(size_sum(size_sum~=0));   
    end
    boot_strap_pick_num = min(mini_val); % Number of trials to pick per iteration (to match bad subjects)
end