% 4. Visualization Check 
function ave_ersp_single_sub_optimized(ave_single_sub, plot_switch, time_range, freq_range, sub_num, fig_ersp)
    % [OPTIMIZATION] The loop and redundant load logic has been removed.
    % We process the pre-calculated ave_single_sub passed in memory.

    types = fieldnames(ave_single_sub);

    if plot_switch
        tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        visualisation_ersp(ave_single_sub.(types{2}).ersp - ave_single_sub.(types{1}).ersp, fig_ersp, time_range, freq_range);
        title("expected - expected duplicated ERSP");
        
        visualisation_ersp(ave_single_sub.(types{3}).ersp - ave_single_sub.(types{2}).ersp, fig_ersp, time_range, freq_range);
        title("unexpected - expected ERSP");
        
        visualisation_ersp(ave_single_sub.(types{5}).ersp - ave_single_sub.(types{4}).ersp, fig_ersp, time_range, freq_range);
        title("right after - right before ERSP");
        
        fig = gcf; % Get current figure handle
        exportgraphics(fig, string("output\all_in_one_pic\boot_average_sub_" + sub_num + ".png"), 'Resolution', 300);
    end 
end
