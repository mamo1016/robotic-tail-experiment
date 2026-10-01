%% plot_loocv_4methods_trimmed.m
% Generate a trimmed 4-panel LOOCV scatter figure for main text Figure 3A.
% Shows: Spearman (failed exemplar) + LASSO + Ridge + SVR (successes).
% Uses the CANONICAL run, so the panel r values match the Figure 3 caption and
% the headline numbers in main.tex (LASSO 0.995, Ridge 0.993, SVR 0.989).
%
% Output: paper/figures/loocv_comparison_4methods.png

clear; clc;
addpath('functions\for_GA');
addpath('functions\for_obtainERD');
config = configure_parameters();

to_dir = fullfile(config.out_dir, 'ga_feature_30runs', 'toResult');
paper_fig_dir = fullfile(fileparts(pwd), 'paper', 'figures');

show_methods = {'spearman', 'lasso', 'ridge', 'svr'};
all_methods  = {'spearman','pca','robust','pls','lasso','ridge','svr'};
panel_titles = {'SPEARMAN', 'LASSO', 'RIDGE', 'SVR'};
panel_labels = {'(a)', '(b)', '(c)', '(d)'};

%% 1. Load the CANONICAL run
% The canonical run is the paper's primary analysis. The 102-run archive in
% ga_feature_30runs/toResult/ is the robustness check only, and is reported
% separately as mean+/-SD in main.tex.
%
% Do NOT swap in a "most representative" run here. This script previously
% searched toResult/ for the run closest to the batch mean, which made the
% panel r values (LASSO 0.976, Ridge 0.948, SVR 0.991) silently disagree with
% the Figure 3 caption and the headline numbers in main.tex. Fixed 2026-07-21.
CANONICAL_FILE = 'results_LOOCV_2026-02-24_15-32.mat';
canonical_path = fullfile(config.out_dir, 'ga_feature', CANONICAL_FILE);

if ~isfile(canonical_path)
    alt_path = fullfile(to_dir, CANONICAL_FILE);
    if isfile(alt_path)
        canonical_path = alt_path;
    else
        error('plot_loocv_4methods:canonicalMissing', ...
              ['Canonical results file not found. Looked in:\n  %s\n  %s\n' ...
               'Set canonical_path manually if the file has moved.'], ...
              canonical_path, alt_path);
    end
end

%% 2. Load
fprintf('Loading canonical run: %s\n', canonical_path);
results = load(canonical_path, 'results').results;

% Report the r values actually plotted, so any drift from the caption is
% visible in the console at generation time.
fprintf('Panel r values (Spearman correlation vs actual):\n');
for mi = 1:numel(show_methods)
    m = show_methods{mi};
    fprintf('  %-9s r = %+.3f\n', m, ...
            corr(results.(m).predicted(:), results.(m).actual(:), 'Type', 'Spearman'));
end

%% 3. Compute axis ranges
fail_vals = [results.spearman.predicted(:); results.spearman.actual(:)];
pad = 0.05 * range(fail_vals);
fail_range = [min(fail_vals) - pad, max(fail_vals) + pad];

succ_vals = [];
for m = {'lasso','ridge','svr'}
    succ_vals = [succ_vals; results.(m{1}).predicted(:); results.(m{1}).actual(:)];
end
pad = 0.05 * range(succ_vals);
succ_range = [min(succ_vals) - pad, max(succ_vals) + pad];

%% 4. Plot 2x2 grid
fig = figure('Color', 'w', 'Position', [50 50 900 750]);
tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:4
    m = show_methods{i};
    if i == 1
        ax_range = fail_range;
    else
        ax_range = succ_range;
    end

    nexttile;

    pred   = results.(m).predicted(:);
    actual = results.(m).actual(:);
    [r_val, p_val] = corr(pred, actual, 'Type', 'Spearman');

    scatter(actual, pred, 60, 'filled', 'MarkerFaceAlpha', 0.7);
    hold on;
    plot(ax_range, ax_range, 'k--', 'LineWidth', 1);
    p_fit = polyfit(actual, pred, 1);
    x_fit = linspace(ax_range(1), ax_range(2), 100);
    plot(x_fit, polyval(p_fit, x_fit), 'r-', 'LineWidth', 1.5);
    hold off;

    xlim(ax_range); ylim(ax_range);
    axis square;

    if p_val < 0.001
        p_str = 'p<0.001';
    else
        p_str = sprintf('p=%.3f', p_val);
    end

    title(sprintf('%s\nr=%.3f, %s', panel_titles{i}, r_val, p_str), ...
          'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Actual CoP response', 'FontSize', 10);
    ylabel('Predicted CoP response', 'FontSize', 10);
    set(gca, 'FontSize', 10, 'LineWidth', 1);

    rng_x = range(ax_range);
    rng_y = range(ax_range);
    text(ax_range(1) + 0.03*rng_x, ax_range(2) - 0.06*rng_y, ...
         panel_labels{i}, 'FontSize', 12, 'FontWeight', 'bold', ...
         'VerticalAlignment', 'top');
end

%% 5. Save
out_path = fullfile(paper_fig_dir, 'loocv_comparison_4methods.png');
exportgraphics(fig, out_path, 'Resolution', 300);
close(fig);
fprintf('Saved: %s\n', out_path);
