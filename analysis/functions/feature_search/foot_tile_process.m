function foot_ave = foot_tile_process(foot)
    foot_sub_list = fieldnames(foot);
    foot_ave = struct();
    for i = 1:length(foot_sub_list)
        if i == 1 
            types = fieldnames(foot.(foot_sub_list{1}));
        end
    
        for j = 1:length(types)
            if i == 1
                foot_ave.(types{j}) = foot.(foot_sub_list{i}).(types{j}).vector;
            elseif i == 25
                foot_ave.(types{j}) = foot_ave.(types{j}) + foot.(foot_sub_list{i}).(types{j}).vector;
                foot_ave.(types{j}) = foot_ave.(types{j})/i;
            else
                foot_ave.(types{j}) = foot_ave.(types{j}) + foot.(foot_sub_list{i}).(types{j}).vector;
            end
        end
    
    end
end