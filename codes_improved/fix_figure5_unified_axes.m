%% fix_figure5_unified_axes.m  (mn-4)
% Regenerate the LOOCV 7-method scatter plot (Figure 5) with two-tier
% unified axis ranges: one for the top row (unregularized methods) and
% one for the bottom row (regularized methods).
%
% Run from: codes_improved

clear; clc;
addpath('functions\for_GA');
addpath('functions\for_obtainERD');
config = configure_parameters();
out_dir = config.out_dir;
save_dir = fullfile(out_dir, 'ga_feature');
paper_fig_dir = 'figures';

%% --- Load the CORRECT LOOCV results (Feb 24 run with 7 methods) ---
% NOTE: The most recent .mat files (Mar 5-6) have broken regression scores
% due to a code change. The Feb 24 15:32 run matches the paper values:
% LASSO r=0.995, RIDGE r=0.993, SVR r=0.989.
loocv_path = fullfile(save_dir, 'results_LOOCV_2026-02-24_15-32.mat');
if ~exist(loocv_path, 'file')
    error('Target LOOCV file not found: %s', loocv_path);
end
fprintf('Loading: %s\n', loocv_path);
load(loocv_path, 'results');

methods = fieldnames(results);
N_methods = length(methods);

%% --- Define two-tier grouping ---
% Top row: unregularized (indices 1-4), Bottom row: regularized (5-7)
n_cols = 4;
top_idx = 1:min(4, N_methods);
bot_idx = 5:N_methods;

%% --- Compute axis range per tier ---
top_vals = [];
for m = top_idx
    pred   = results.(methods{m}).predicted;
    actual = results.(methods{m}).actual;
    top_vals = [top_vals; pred(:); actual(:)];
end
top_range = [min(top_vals) - 0.05*range(top_vals), ...
             max(top_vals) + 0.05*range(top_vals)];

bot_vals = [];
for m = bot_idx
    pred   = results.(methods{m}).predicted;
    actual = results.(methods{m}).actual;
    bot_vals = [bot_vals; pred(:); actual(:)];
end
bot_range = [min(bot_vals) - 0.05*range(bot_vals), ...
             max(bot_vals) + 0.05*range(bot_vals)];

fprintf('Top row axis range: [%.3f, %.3f]\n', top_range(1), top_range(2));
fprintf('Bottom row axis range: [%.3f, %.3f]\n', bot_range(1), bot_range(2));

%% --- Plot with two-tier axes ---
fig = figure('Color', 'w', 'Position', [50 50 300*n_cols 600]);
tl = tiledlayout(2, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

for m = 1:N_methods
    % Determine row and axis range
    if m <= 4
        ax_range = top_range;
    else
        ax_range = bot_range;
    end

    nexttile;

    pred   = results.(methods{m}).predicted;
    actual = results.(methods{m}).actual;
    r_val  = corr(pred, actual, 'Type', 'Spearman');
    [~, p_val] = corr(pred, actual, 'Type', 'Spearman');

    scatter(actual, pred, 50, 'filled', 'MarkerFaceAlpha', 0.7);
    hold on;

    % Unity line
    plot(ax_range, ax_range, 'k--', 'LineWidth', 1);

    % Regression line
    p_fit = polyfit(actual, pred, 1);
    x_fit = linspace(ax_range(1), ax_range(2), 100);
    plot(x_fit, polyval(p_fit, x_fit), 'r-', 'LineWidth', 1.5);
    hold off;

    xlim(ax_range);
    ylim(ax_range);
    axis square;

    % Title with stats
    if p_val < 0.001
        p_str = 'p<0.001';
    else
        p_str = sprintf('p=%.3f', p_val);
    end
    title(sprintf('%s\nr=%.3f, %s', upper(methods{m}), r_val, p_str), ...
        'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Actual', 'FontSize', 9);
    ylabel('Predicted', 'FontSize', 9);
    set(gca, 'FontSize', 9);
end

title(tl, 'Cross-Validated Brain--Body Predictions', 'FontSize', 12);

%% --- Save ---
fig_path = fullfile(paper_fig_dir, 'loocv_comparison_7methods.png');
exportgraphics(fig, fig_path, 'Resolution', 300);
fprintf('Figure saved to: %s\n', fig_path);
fprintf('=== DONE ===\n');
