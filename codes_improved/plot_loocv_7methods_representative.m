%% plot_loocv_7methods_representative.m
% Generate the 7-method LOOCV scatter figure (Supplementary Fig. S9) from the
% CANONICAL run, so the panel r values match the supplementary caption and the
% headline numbers in main.tex (LASSO 0.995, Ridge 0.993, SVR 0.989).
%
% NOTE ON THE FILE NAME: this script no longer picks a "representative" run.
% It previously searched ga_feature_30runs/bin/ for the run whose per-method
% r-values sat closest to the batch mean, which made the published panels
% (LASSO 0.960, Ridge 0.980, SVR 0.957) silently disagree with every number
% quoted in supplementary.tex. That is the same defect that was fixed in
% plot_loocv_4methods_trimmed.m. Fixed 2026-07-25. The 102-run archive is the
% robustness check only and is reported separately as mean+/-SD.
%
% Output: <repo>/paper/figures/loocv_comparison_7methods.png

clear; clc;
addpath('functions\for_GA');
addpath('functions\for_obtainERD');
config = configure_parameters();

paper_fig_dir = fullfile(fileparts(pwd), 'paper', 'figures');

methods   = {'spearman','pca','robust','pls','lasso','ridge','svr'};
n_methods = numel(methods);

%% 1. Resolve the canonical run
CANONICAL_FILE = 'results_LOOCV_2026-02-24_15-32.mat';
canonical_path = fullfile(config.out_dir, 'ga_feature', CANONICAL_FILE);

if ~isfile(canonical_path)
    alt_path = fullfile(config.out_dir, 'ga_feature_30runs', 'toResult', CANONICAL_FILE);
    if isfile(alt_path)
        canonical_path = alt_path;
    else
        error('plot_loocv_7methods:canonicalMissing', ...
              ['Canonical results file not found. Looked in:\n  %s\n  %s\n' ...
               'Set canonical_path manually if the file has moved.'], ...
              canonical_path, alt_path);
    end
end

%% 2. Load the canonical run
fprintf('Loading canonical run: %s\n', canonical_path);
results = load(canonical_path, 'results').results;

% Report the r/p values actually plotted, so any drift from the supplementary
% caption is visible in the console at generation time.
fprintf('\nPanel values (Spearman, predicted vs actual):\n');
for mi = 1:n_methods
    m = methods{mi};
    [r_c, p_c] = corr(results.(m).predicted(:), results.(m).actual(:), 'Type', 'Spearman');
    fprintf('  %-9s r = %+.4f   p = %.4f\n', m, r_c, p_c);
end
fprintf('\n');

%% 3. Two-tier unified axis ranges
top_idx = 1:4;   % unregularised
bot_idx = 5:7;   % regularised

top_vals = [];
for m = top_idx
    top_vals = [top_vals; ...
        results.(methods{m}).predicted(:); ...
        results.(methods{m}).actual(:)];  %#ok<AGROW>
end
pad = 0.05 * range(top_vals);
top_range = [min(top_vals) - pad, max(top_vals) + pad];

bot_vals = [];
for m = bot_idx
    bot_vals = [bot_vals; ...
        results.(methods{m}).predicted(:); ...
        results.(methods{m}).actual(:)];  %#ok<AGROW>
end
pad = 0.05 * range(bot_vals);
bot_range = [min(bot_vals) - pad, max(bot_vals) + pad];

fprintf('Top row axis range : [%.3f, %.3f]\n', top_range(1), top_range(2));
fprintf('Bottom row axis range: [%.3f, %.3f]\n', bot_range(1), bot_range(2));

%% 4. Plot
n_cols = 4;
fig = figure('Color', 'w', 'Position', [50 50 300*n_cols 600]);
tl  = tiledlayout(2, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

panel_labels = 'abcdefg';

for m = 1:n_methods
    ax_range = bot_range;
    if m <= 4
        ax_range = top_range;
    end

    nexttile;

    pred   = results.(methods{m}).predicted(:);
    actual = results.(methods{m}).actual(:);
    [r_val, p_val] = corr(pred, actual, 'Type', 'Spearman');

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

    % p-value string
    if p_val < 0.001
        p_str = 'p<0.001';
    else
        p_str = sprintf('p=%.3f', p_val);
    end

    title(sprintf('%s\nr=%.3f, %s', upper(methods{m}), r_val, p_str), ...
          'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Actual CoP response',    'FontSize', 9);
    ylabel('Predicted CoP response', 'FontSize', 9);
    set(gca, 'FontSize', 9);

    % Panel label (a)–(g) in top-left corner
    rng_x = range(ax_range);
    rng_y = range(ax_range);
    text(ax_range(1) + 0.03*rng_x, ax_range(2) - 0.06*rng_y, ...
         sprintf('(%s)', panel_labels(m)), ...
         'FontSize', 11, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
end

% Plain hyphen: the LaTeX en-dash form "--" rendered literally as two hyphens
% in the exported PNG.
title(tl, 'Cross-validated brain-body predictions', 'FontSize', 12);

%% 5. Save
fig_path = fullfile(paper_fig_dir, 'loocv_comparison_7methods.png');
exportgraphics(fig, fig_path, 'Resolution', 300);
fprintf('\nFigure saved to: %s\n', fig_path);
fprintf('=== DONE ===\n');
