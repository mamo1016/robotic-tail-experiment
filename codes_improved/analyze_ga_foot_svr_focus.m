% analyze_ga_foot_svr_focus.m
% Deep dive into SVR foot feature selection across 102 LOOCV runs.
%
% SVR is the best performer in LOOCV (r=0.978±0.012).
% This script:
%   1. Shows which foot features SVR selected (with clean labels)
%   2. Compares SVR selections across 102 runs
%   3. Identifies consensus SVR features (selected in most runs)
%   4. Links back to trial-by-trial trends (for next analysis)

clear; close all; clc;
addpath('functions\for_obtainERD');
config = configure_parameters();
out_dir = config.out_dir;

%% ---- Build foot feature labels (structured for readability) ----
conditions   = {'Expected', 'Unexpected', 'RightBefore', 'RightAfter'};
cond_short   = {'Exp', 'Unexp', 'RBefore', 'RAfter'};

time_windows = {'Baseline', 'Reflex', 'EarlyComp', 'LateComp', 'TotalPost'};
win_short    = {'Base', 'Refl', 'Early', 'Late', 'Total'};

metric_names = {'RMS', 'Area'};
metric_short = {'RMS', 'Area'};

% Build feature labels (full and short)
foot_feature_labels_full = {};
foot_feature_labels_short = {};
feat_idx = 0;
for c = 1:length(conditions)
    for w = 1:length(time_windows)
        for m = 1:length(metric_names)
            feat_idx = feat_idx + 1;
            foot_feature_labels_full{feat_idx} = sprintf('%s_%s_%s', ...
                conditions{c}, time_windows{w}, metric_names{m});
            foot_feature_labels_short{feat_idx} = sprintf('%s_%s_%s', ...
                cond_short{c}, win_short{w}, metric_short{m});
        end
    end
end

n_features = length(foot_feature_labels_full);
fprintf('Loaded %d foot feature labels (full and short).\n', n_features);

%% ---- Find all LOOCV result files ----
result_dir = fullfile(out_dir, 'ga_feature_30runs', 'toResult');
d = dir(fullfile(result_dir, 'results_LOOCV_*.mat'));
n_results = length(d);
fprintf('Found %d result files.\n', n_results);

%% ---- Extract SVR foot feature selection ----
svr_selection = zeros(n_results, n_features);
svr_corr = zeros(n_results, 1);

fprintf('Extracting SVR foot feature selection from %d runs...\n', n_results);

for r = 1:n_results
    result_file = fullfile(result_dir, d(r).name);

    try
        load(result_file, 'results');
    catch
        continue;
    end

    if isfield(results, 'svr')
        if isfield(results.svr, 'mask_y')
            svr_selection(r, :) = results.svr.mask_y;
        end
        % Also extract correlation for this run
        if isfield(results.svr, 'predicted') && isfield(results.svr, 'actual')
            pred = results.svr.predicted(:);
            actual = results.svr.actual(:);
            c = corr(pred, actual);
            if isscalar(c)
                svr_corr(r) = c;
            else
                svr_corr(r) = c(1, 1);  % take first element if matrix returned
            end
        end
    end

    if mod(r, 25) == 0
        fprintf('  Processed %d/%d runs...\n', r, n_results);
    end
end

fprintf('Extraction complete.\n\n');

%% ---- Compute SVR feature selection statistics ----
svr_freq = sum(svr_selection, 1) / n_results * 100;  % % of runs where selected
[svr_freq_sorted, sort_idx] = sort(svr_freq, 'descend');

%% ---- Create output folder ----
fig_dir = fullfile(out_dir, 'ga_foot_feature_analysis');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% ---- Figure 1: SVR feature selection frequency (all 40 features) ----
fig = figure('Position', [100 100 1400 600], 'Color', 'w');

% All 40 features, sorted by frequency
bar(1:n_features, svr_freq_sorted, 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none');
hold on;
yline(50, 'r--', 'LineWidth', 2, 'Label', '50% threshold (consensus)');
hold off;

xlabel('Foot Feature (ranked by selection frequency)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('SVR Selection Frequency (% of 102 runs)', 'FontSize', 12, 'FontWeight', 'bold');
title('SVR Foot Feature Selection Across 102 Runs', 'FontSize', 13, 'FontWeight', 'bold');

set(gca, 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 11);
grid on; grid minor;
xlim([0 n_features+1]);
ylim([0 105]);

svg_freq_path = fullfile(fig_dir, '04_svr_all_features_ranked.png');
exportgraphics(fig, svg_freq_path, 'Resolution', 300);
close(fig);

fprintf('Saved: %s\n', svg_freq_path);

%% ---- Figure 2: SVR feature selection — structured by condition/window ----
% Reshape into a 3D array for better visualization: [4 conditions x 5 windows x 2 metrics]
svr_freq_3d = reshape(svr_freq, [4, 5, 2]);  % [condition x window x metric]

fig = figure('Position', [100 100 1200 800], 'Color', 'w');

for c = 1:4
    for m = 1:2
        subplot(4, 2, (c-1)*2 + m);

        % Heatmap: windows (5) as x-axis, selection frequency as color
        data_slice = squeeze(svr_freq_3d(c, :, m));
        imagesc(data_slice);
        colormap('hot');
        caxis([0 100]);

        if m == 2
            cbar = colorbar;
            cbar.Label.String = '% Selected';
        end

        set(gca, 'XTick', 1:5, 'XTickLabel', win_short, 'FontSize', 10);
        set(gca, 'YTick', [1], 'YTickLabel', {metric_short{m}}, 'FontSize', 10);

        title(sprintf('%s (%s)', conditions{c}, metric_names{m}), ...
            'FontSize', 11, 'FontWeight', 'bold');

        % Add percentage text
        hold on;
        for w = 1:5
            pct = data_slice(w);
            if pct > 0
                text(w, 1, sprintf('%.0f%%', pct), 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', 'Color', 'w', 'FontSize', 11, ...
                    'FontWeight', 'bold');
            end
        end
        hold off;
    end
end

sgtitle('SVR Feature Selection: Structured View (Condition × Window × Metric)', ...
    'FontSize', 13, 'FontWeight', 'bold');

svg_struct_path = fullfile(fig_dir, '05_svr_structured_heatmap.png');
exportgraphics(fig, svg_struct_path, 'Resolution', 300);
close(fig);

fprintf('Saved: %s\n', svg_struct_path);

%% ---- Figure 3: Top 15 SVR features (clean labels) ----
fig = figure('Position', [100 100 1000 600], 'Color', 'w');

n_top = 15;
top_features_idx = sort_idx(1:n_top);
top_features_labels = foot_feature_labels_short(top_features_idx);
top_features_freq = svr_freq(top_features_idx);

% Sort by frequency within top 15
[top_freq_sorted, top_sort_idx] = sort(top_features_freq, 'descend');
top_labels_sorted = top_features_labels(top_sort_idx);

bar(1:n_top, top_freq_sorted, 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'k', 'LineWidth', 1);
set(gca, 'XTick', 1:n_top, 'XTickLabel', top_labels_sorted, ...
    'XTickLabelRotation', 45, 'FontSize', 11);
ylabel('Selection Frequency (% of 102 runs)', 'FontSize', 12, 'FontWeight', 'bold');
title('SVR: Top 15 Selected Foot Features', 'FontSize', 13, 'FontWeight', 'bold');

grid on; grid minor;
set(gca, 'Box', 'off', 'LineWidth', 1.2);
ylim([0 105]);

svg_top_path = fullfile(fig_dir, '06_svr_top15_features.png');
exportgraphics(fig, svg_top_path, 'Resolution', 300);
close(fig);

fprintf('Saved: %s\n', svg_top_path);

%% ---- Figure 4: SVR correlation stability across runs ----
fig = figure('Position', [100 100 900 500], 'Color', 'w');

scatter(1:n_results, svr_corr, 50, [0.3 0.5 0.9], 'filled');
hold on;
mean_corr = mean(svr_corr);
yline(mean_corr, 'r-', 'LineWidth', 2, 'Label', sprintf('Mean r=%.4f', mean_corr));
std_corr = std(svr_corr);
yline(mean_corr + std_corr, 'r--', 'LineWidth', 1);
yline(mean_corr - std_corr, 'r--', 'LineWidth', 1);
hold off;

xlabel('Run Number (1–102)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Brain–Body Correlation (r)', 'FontSize', 12, 'FontWeight', 'bold');
title('SVR LOOCV Performance Stability Across 102 Runs', ...
    'FontSize', 13, 'FontWeight', 'bold');

set(gca, 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 11);
grid on; grid minor;
ylim([0.95 1.0]);

svg_corr_path = fullfile(fig_dir, '07_svr_correlation_stability.png');
exportgraphics(fig, svg_corr_path, 'Resolution', 300);
close(fig);

fprintf('Saved: %s\n', svg_corr_path);

%% ---- Summary table: SVR feature selection ----
svr_summary = table();
svr_summary.Feature_Full = foot_feature_labels_full';
svr_summary.Feature_Short = foot_feature_labels_short';
svr_summary.SVR_Selection_Percent = svr_freq';

% Add condition/window/metric breakdown (matching feature loop structure)
condition_list = {};
window_list = {};
metric_list = {};
for c = 1:length(conditions)
    for w = 1:length(time_windows)
        for m = 1:length(metric_names)
            condition_list = [condition_list, conditions{c}];
            window_list = [window_list, time_windows{w}];
            metric_list = [metric_list, metric_names{m}];
        end
    end
end

svr_summary.Condition = condition_list';
svr_summary.TimeWindow = window_list';
svr_summary.Metric = metric_list';

% Sort by selection frequency
[~, sort_order] = sort(svr_summary.SVR_Selection_Percent, 'descend');
svr_summary = svr_summary(sort_order, :);

csv_path = fullfile(fig_dir, 'svr_foot_feature_selection_summary.csv');
writetable(svr_summary, csv_path);

mat_path = fullfile(fig_dir, 'svr_foot_feature_selection_summary.mat');
save(mat_path, 'svr_summary', 'svr_freq', 'svr_corr', 'foot_feature_labels_full', ...
    'foot_feature_labels_short');

fprintf('Saved: %s\n', csv_path);
fprintf('Saved: %s\n', mat_path);

%% ---- Console report ----
fprintf('\n========================================\n');
fprintf('  SVR FOOT FEATURE SELECTION ANALYSIS\n');
fprintf('========================================\n\n');

fprintf('Runs analysed: %d\n', n_results);
fprintf('SVR correlation: %.4f ± %.4f (range: %.4f–%.4f)\n', ...
    mean(svr_corr), std(svr_corr), min(svr_corr), max(svr_corr));

fprintf('\n--- CONSENSUS FEATURES (>50%% of runs) ---\n');
consensus_idx = find(svr_freq > 50);
if ~isempty(consensus_idx)
    [consensus_freq, consensus_sort] = sort(svr_freq(consensus_idx), 'descend');
    consensus_features = foot_feature_labels_short(consensus_idx(consensus_sort));

    for i = 1:length(consensus_idx)
        fprintf('[%2d] %30s  %.1f%%\n', i, consensus_features{i}, consensus_freq(i));
    end
else
    fprintf('(None)\n');
end

fprintf('\n--- TOP 20 SVR SELECTED FEATURES ---\n');
for i = 1:min(20, n_features)
    feat_idx = sort_idx(i);
    fprintf('[%2d] %30s  %.1f%%\n', i, foot_feature_labels_short{feat_idx}, svr_freq(feat_idx));
end

fprintf('\n--- TOP 20 NON-SELECTED FEATURES ---\n');
for i = 1:min(20, n_features)
    feat_idx = sort_idx(n_features - i + 1);
    fprintf('[%2d] %30s  %.1f%%\n', i, foot_feature_labels_short{feat_idx}, svr_freq(feat_idx));
end

fprintf('\n========================================\n');
fprintf('  Output folder: %s\n', fig_dir);
fprintf('========================================\n');
