function visualisation_ersp(ersp_vis, fig_ersp, time_range, freq_range, clims)
    arguments
        ersp_vis
        fig_ersp
        time_range
        freq_range
        clims = [-1 1]   % is OPTIONAL. If missing, it becomes [-1 1].
    end
%% 5. Visualization
    % figure(fig_ersp);
    nexttile;
    % clims = [-1 1]; % Color limits in dB (Adjust based on your signal strength)
    imagesc(time_range, freq_range, ersp_vis, clims); 
    axis xy; % Flip Y-axis so low freq is at bottom
    cmap = [linspace(0, 1, 128)', linspace(0, 1, 128)', ones(128, 1); % Blue to White
        ones(128, 1), linspace(1, 0, 128)', linspace(1, 0, 128)']; % White to Red
    colormap(cmap); % 'jet' or 'parula' are standard
    colorbar;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');

    % Add a vertical line at Time 0 (Event)
    hold on; xline(0, 'k--', 'LineWidth', 0.4); hold off;

    % Limit view to relevant frequencies (e.g., 1-40 Hz)
    % ylim([5 30]);