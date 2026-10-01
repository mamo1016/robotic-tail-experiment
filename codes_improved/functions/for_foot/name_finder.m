function current_type_name_foot = name_finder(current_type_name)
    
    if     current_type_name == "expected_duplicated"
        current_type_name_foot = "expected_duplicated";
    
    elseif current_type_name == "expected"
        current_type_name_foot = "epochs_expected";
    
    elseif current_type_name == "unexpected"
        current_type_name_foot = "epochs_unexpected";

    elseif current_type_name == "right_before"
        current_type_name_foot = "epochs_right_before";

    elseif current_type_name == "right_after"
        current_type_name_foot = "epochs_right_after";
    
    end
    
end