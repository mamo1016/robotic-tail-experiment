% 4. Visualization Check
function average = average_ersp(ersp_epoch_bootstraped, plot_switch,time_range, freq_range)
fields = fieldnames(ersp_epoch_bootstraped);

if plot_switch
    fig_ersp = figure('Color', 'w'); % 'Name' sets the window title
end

for struct_num = 1:length(fields)
    current_name = fields{struct_num};           % Get the name (string)  
    % trial_type = fieldnames(ersp_epoch_bootstraped.(current_name)); 
    ersp_vis_expected = ersp_epoch_bootstraped.(current_name).expected.final_ersp;
    ersp_vis_unexpected = ersp_epoch_bootstraped.(current_name).unexpected.final_ersp;
    ersp_vis_right_before = ersp_epoch_bootstraped.(current_name).right_before.final_ersp;
    ersp_vis_right_after = ersp_epoch_bootstraped.(current_name).right_after.final_ersp;
    ersp_vis_expected_duplicated = ersp_epoch_bootstraped.(current_name).expected_duplicated.final_ersp;
    if struct_num == 1
        ersp_vis_unexpected_expected_average = zeros(size(ersp_vis_unexpected));
        ersp_vis_expected_expected_average = zeros(size(ersp_vis_unexpected));
        ersp_vis_after_before_average = zeros(size(ersp_vis_unexpected));
    end
    if plot_switch
        tiledlayout(4,1, 'TileSpacing', 'compact', 'Padding', 'compact');
        visualisation_ersp(ersp_vis_expected, fig_ersp, time_range, freq_range); title("expected ERSP");
        visualisation_ersp(ersp_vis_unexpected, fig_ersp, time_range, freq_range); title("unexpected ERSP");
        visualisation_ersp(ersp_vis_unexpected-ersp_vis_expected, fig_ersp, time_range, freq_range); title("unexpected - expected ERSP");
        visualisation_ersp(ersp_vis_expected_duplicated-ersp_vis_expected, fig_ersp, time_range, freq_range); title("expected - expected ERSP");
        fig = gcf; % Get current figure handle
        exportgraphics(fig, string("output\unexpected_expected_pic\"+current_name+".png"), 'Resolution', 300);
    end
    ersp_vis_unexpected_expected_average = ersp_vis_unexpected_expected_average + (ersp_vis_unexpected-ersp_vis_expected);
    ersp_vis_expected_expected_average = ersp_vis_expected_expected_average + (ersp_vis_expected_duplicated-ersp_vis_expected);

    if plot_switch
        tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');   
        visualisation_ersp(ersp_vis_right_before, fig_ersp, time_range, freq_range); title("right before ERSP");
        visualisation_ersp(ersp_vis_right_after, fig_ersp, time_range, freq_range); title("right after ERSP");
        visualisation_ersp(ersp_vis_right_after-ersp_vis_right_before, fig_ersp, time_range, freq_range); title("right after - right before ERSP");
        fig = gcf; % Get current figure handle
        exportgraphics(fig, string("output\rightafter_rightbefore_pic\"+current_name+".png"), 'Resolution', 300);
    end
    ersp_vis_after_before_average = ersp_vis_after_before_average + (ersp_vis_right_after-ersp_vis_right_before);
           
% check all trial types and dir=fference in each subject, and summarise for all subject
end
average = struct();
average.ersp_vis_unexpected_expected_average = ersp_vis_unexpected_expected_average/struct_num;
average.ersp_vis_expected_expected_average = ersp_vis_expected_expected_average/struct_num;
average.ersp_vis_after_before_average = ersp_vis_after_before_average/struct_num;