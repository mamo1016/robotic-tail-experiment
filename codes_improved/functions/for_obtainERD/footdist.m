function datasets = footdist(datasets)
    for i=1:2
        foot_raw = datasets.foot_raw(:,i);
        foot_time = datasets.foot_time; 
        
        foot_raw_dist=zeros(length(foot_raw),5);
        foot_time_dist=zeros(length(foot_raw),1);
        row_num=1;
        %1-5,                       1-... skip6x
        foot_indent=find(foot_raw(1:6)==255);
        for foot_raw_loop=1:length(foot_raw)
            if foot_raw(foot_raw_loop)==255
                foot_time_dist(row_num,1)=foot_time(foot_raw_loop,1);
                row_num=row_num+1;
            else
                % foot_raw_dist(row_num,rem(foot_raw_loop+foot_indent,6))=foot_raw(foot_raw_loop);
                foot_raw_dist(row_num,rem(foot_raw_loop+(6-foot_indent),6))=foot_raw(foot_raw_loop);
            end
        end
        foot_raw_dist(row_num+1:end,:)=[];
        foot_time_dist(row_num+1:end,:)=[];
        
        if i==1
            datasets.foot_raw_dist_L = foot_raw_dist;
            datasets.foot_time_dist_L = foot_time_dist;
        elseif i==2
            datasets.foot_raw_dist_R = foot_raw_dist;
            datasets.foot_time_dist_R = foot_time_dist;
        end
    end

    % Ensure foot_time_dist arrays are properly sized and aligned
    if ~isempty(datasets.foot_time_dist_L) && ~isempty(datasets.foot_time_dist_R)
        % Make sure both have the same length (use the shorter one)
        min_length = min(length(datasets.foot_time_dist_L), length(datasets.foot_time_dist_R));
        datasets.foot_time_dist_L = datasets.foot_time_dist_L(1:min_length);
        datasets.foot_time_dist_R = datasets.foot_time_dist_R(1:min_length);
        datasets.foot_raw_dist_L = datasets.foot_raw_dist_L(1:min_length, :);
        datasets.foot_raw_dist_R = datasets.foot_raw_dist_R(1:min_length, :);
    end
%%
% a=zeros(20,1);
% for i=1:20
%     a(i)=rem(i+(6-5),6);
% end
% 4
% 5
% 255
% 1
% 2
% 3
% 4
% 5
% 255
% 1
% 2
% 3
% 4
% 5
% 255