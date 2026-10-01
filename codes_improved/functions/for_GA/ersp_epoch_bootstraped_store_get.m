function ersp_epoch_bootstraped_store = ersp_epoch_bootstraped_store_get()
    ersp_epoch_bootstraped_store = struct();
    for loop = 1:211
        loop_3dig = sprintf('%03d', loop); 
        ersp_epoch_bootstraped = load(string("output\ersp\ersp_epoch_bootstraped"+loop_3dig+".mat")).ersp_epoch_bootstraped;
        fields = fieldnames(ersp_epoch_bootstraped);           
    
        for struct_num = 1:length(fields)
            current_name = fields{struct_num};           % Get the name (string)  
            if loop == 1
                ersp_epoch_bootstraped_store.(current_name).expected_ersp_ave =       zeros(size(ersp_epoch_bootstraped.(current_name).expected.final_ersp));
                ersp_epoch_bootstraped_store.(current_name).unexpected_ersp_ave =     zeros(size(ersp_epoch_bootstraped.(current_name).expected.final_ersp));
                ersp_epoch_bootstraped_store.(current_name).right_before_ersp_ave =   zeros(size(ersp_epoch_bootstraped.(current_name).expected.final_ersp));
                ersp_epoch_bootstraped_store.(current_name).right_after_ersp_ave =    zeros(size(ersp_epoch_bootstraped.(current_name).expected.final_ersp));
                ersp_epoch_bootstraped_store.(current_name).expected_duplicated_ave = zeros(size(ersp_epoch_bootstraped.(current_name).expected.final_ersp));
    
            end
            ersp_epoch_bootstraped_store.(current_name).expected_ersp_ave =       ersp_epoch_bootstraped_store.(current_name).expected_ersp_ave +                   ersp_epoch_bootstraped.(current_name).expected.final_ersp;
            ersp_epoch_bootstraped_store.(current_name).unexpected_ersp_ave =     ersp_epoch_bootstraped_store.(current_name).unexpected_ersp_ave +                 ersp_epoch_bootstraped.(current_name).unexpected.final_ersp;
            ersp_epoch_bootstraped_store.(current_name).right_before_ersp_ave =   ersp_epoch_bootstraped_store.(current_name).right_before_ersp_ave +               ersp_epoch_bootstraped.(current_name).right_before.final_ersp;
            ersp_epoch_bootstraped_store.(current_name).right_after_ersp_ave =    ersp_epoch_bootstraped_store.(current_name).right_after_ersp_ave +                ersp_epoch_bootstraped.(current_name).right_after.final_ersp;
            ersp_epoch_bootstraped_store.(current_name).expected_duplicated_ave = ersp_epoch_bootstraped_store.(current_name).expected_duplicated_ave +             ersp_epoch_bootstraped.(current_name).expected_duplicated.final_ersp;
        end
    end
    
    for struct_num = 1:length(fields)
        current_name = fields{struct_num};           % Get the name (string)  
    
        ersp_epoch_bootstraped_store.(current_name).expected_ersp_ave =       ersp_epoch_bootstraped_store.(current_name).expected_ersp_ave / loop;
        ersp_epoch_bootstraped_store.(current_name).unexpected_ersp_ave =     ersp_epoch_bootstraped_store.(current_name).unexpected_ersp_ave / loop;
        ersp_epoch_bootstraped_store.(current_name).right_before_ersp_ave =   ersp_epoch_bootstraped_store.(current_name).right_before_ersp_ave / loop;
        ersp_epoch_bootstraped_store.(current_name).right_after_ersp_ave =    ersp_epoch_bootstraped_store.(current_name).right_after_ersp_ave / loop;
        ersp_epoch_bootstraped_store.(current_name).expected_duplicated_ave = ersp_epoch_bootstraped_store.(current_name).expected_duplicated_ave / loop;
    end
end