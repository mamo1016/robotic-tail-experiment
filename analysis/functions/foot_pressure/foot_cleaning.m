function [foot_cleaned] = foot_cleaning(config)
% List of subjects to process
epoch_list = dir(fullfile(config.epochDir, 'Sub*'));

file_num = length(epoch_list);
% fig_epoch_foot = figure('Color', 'w'); % 'Name' sets the window title
% figure(fig_epoch_foot);
% hold on

types = ["epochs_expected", "epochs_unexpected", "epochs_right_before", "epochs_right_after"];
foot_cleaned = struct();

for num = 1:file_num
    epoch_file = fullfile(config.epochDir, epoch_list(num).name);
    epoch = load(epoch_file).epoch_data;
    for type_num = 1:4
        data_store_idx=1;
        data=epoch.(types(type_num)).foot;
        data_store = zeros(size(data));
        for i = 1:size(data,1)
            df = diff(data(i,:));
            ddf = diff(df);
            ddf_remove = 0;
            max_ddf_remove = 0;
            for ddf_idx=1:length(ddf)
                if abs(ddf(1,ddf_idx))<=0.001
                    ddf_remove = ddf_remove + 1;
                    if ddf_remove>max_ddf_remove
                        max_ddf_remove = ddf_remove;
                    end
                else
                    ddf_remove = 0;                    
                end
            end
            % 1. Find the starting positions of each new run
            run_starts = [1, find(diff(ddf) ~= 0) + 1];
            % 2. Calculate the length of each run
            run_lengths = diff([run_starts, numel(ddf) + 1]);

            % if max(run_lengths) > 1000
            %     plot(data(i,:),'b');
            % elseif max_ddf_remove>1000
            %     plot(data(i,:),'r');
            % end

            if max(run_lengths) < 1000 && max_ddf_remove < 1000
                data_store(data_store_idx,:)=data(i,:);
                data_store_idx = data_store_idx + 1;     
                % plot(data(i,:));
            end
        end
        data_store(data_store_idx:end,:)=[];
        foot_cleaned.(erase(epoch_list(num).name, '.mat')).(types(type_num))= data_store;
    end    
end