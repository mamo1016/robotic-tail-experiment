% analyze_ga_foot_feature_selection.m
% Analyze which foot CoP features are selected by the GA across all 7 methods
% and across all 102 independent LOOCV runs.
%
% Output:
%   1. Heatmap: 7 methods x 40 features — selection frequency (% of 102 runs)
%   2. Consensus features: which foot features are consistently selected
%   3. Feature ranking: top features by selection frequency per method
%   4. Comparison with trial-by-trial trends (optional, for next step)

clear; close all; clc;
addpath('functions\for_obtainERD');
config = configure_parameters();
out_dir = config.out_dir;

%% ---- Load foot feature labels ----
load(fullfile(out_dir, 'foot', 'foot_cleaned_session_combined_filtered.mat'), ...
     'foot_filtered');
sub_names = fieldnames(foot_filtered);
num_subs = length(sub_names);

% Regenerate foot feature labels (same as in generate_foot_features.m)
conditions   = {'epochs_expected', 'epochs_unexpected', ...
                'epochs_right_before', 'epochs_right_after'};
cond_labels  = {'Expected', 'Unexpected', 'RightBefore', 'RightAfter'};

time_windows = {[-0.5, 0], [0, 0.5], [0.5, 1.0], [1.0, 1.5], [0, 1.5]};
win_labels   = {'Baseline', 'Reflex', 'EarlyComp', 'LateComp', 'TotalPost'};

metric_names = {'RMS', 'Area'};

% Build feature labels (short and long versions)
cond_short = {'Exp', 'Unexp', 'RBefore', 'RAfter'};
win_short  = {'Base', 'Refl', 'Early', 'Late', 'Total'};

foot_feature_labels_long = {};
foot_feature_labels_short = {};
feat_idx = 0;
for c = 1:length(conditions)
    for w = 1:length(time_windows)
        for m = 1:length(metric_names)
            feat_idx = feat_idx + 1;
            foot_feature_labels_long{feat_idx} = sprintf('%s_%s_%s', ...
                cond_labels{c}, win_labels{w}, metric_names{m});
            foot_feature_labels_short{feat_idx} = sprintf('[%d] %s_%s_%s', ...
                feat_idx, cond_short{c}, win_short{w}, metric_names{m});
        end
    end
end

fprintf('Loaded %d foot feature labels.\n', length(foot_feature_labels_long));

%% ---- Find all LOOCV result files ----
result_dir = fullfile(out_dir, 'ga_feature_30runs', 'toResult');
if ~exist(result_dir, 'dir')
    error('Result directory not found: %s', result_dir);
end

d = dir(fullfile(result_dir, 'results_LOOCV_*.mat'));
n_results = length(d);
fprintf('Found %d result files in %s\n', n_results, result_dir);

if n_results == 0
    error('No result files found!');
end

%% ---- Load all results and aggregate foot feature selection ----
methods = {'spearman', 'pca', 'robust', 'pls', 'lasso', 'ridge', 'svr'};
n_methods = length(methods);
n_features = 40;

% Accumulator: [n_methods x n_features] — counts across all runs
selection_count = zeros(n_methods, n_features);
selection_details = cell(n_methods, n_features);  % store details per (method, feature)

fprintf('\nProcessing %d result files...\n', n_results);

for r = 1:n_results
    result_file = fullfile(result_dir, d(r).name);

    try
        load(result_file, 'results');
    catch ME
        fprintf('  [%3d] ERROR loading %s: %s\n', r, d(r).name, ME.message);
        continue;
    end

    % Extract foot feature selection for each method
    for m = 1:n_methods
        method = methods{m};

        if ~isfield(results, method)
            fprintf('  [%3d] WARNING: method %s not in results\n', r, method);
            continue;
        end

        if isfield(results.(method), 'mask_y')
            mask = results.(method).mask_y;
            selection_count(m, :) = selection_count(m, :) + mask;
        end

        if isfield(results.(method), 'mask_y_counts')
            counts = results.(method).mask_y_counts;
            % Store per-fold counts for later analysis
            if isempty(selection_details{m, 1})
                selection_details{m, 1} = counts;
            else
                % Accumulate across runs (we'll average later)
                selection_details{m, 1} = selection_details{m, 1} + counts;
            end
        end
    end

    if mod(r, 20) == 0
        fprintf('  Processed %d/%d files...\n', r, n_results);
    end
end

% Convert counts to percentages
selection_pct = (selection_count / n_results) * 100;

fprintf('Aggregation complete.\n');

%% ---- Create output folder ----
fig_dir = fullfile(out_dir, 'ga_foot_feature_analysis');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% ---- Figure 1: Heatmap of foot feature selection (7 methods x 40 features) ----
fig = figure('Position', [100 100 1200 600], 'Color', 'w');

% Heatmap
imagesc(selection_pct);
colormap('hot');
cbar = colorbar;
cbar.Label.String = 'Selection Frequency (%)';
cbar.Label.FontSize = 12;

% Method labels (y-axis)
set(gca, 'YTick', 1:n_methods, 'YTickLabel', methods, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'XTick', 1:n_features, 'XTickLabel', foot_feature_labels_short, ...
    'XTickLabelRotation', 90, 'FontSize', 7);

xlabel('Foot Feature (40 features: Cond_Window_Metric)', ...
    'FontSize', 12, 'FontWeight', 'bold');
ylabel('Regression Method', 'FontSize', 12, 'FontWeight', 'bold');
title('GA Foot Feature Selection Frequency Across 102 LOOCV Runs', ...
    'FontSize', 13, 'FontWeight', 'bold');

% Add percentage text on each cell
hold on;
for m = 1:n_methods
    for f = 1:n_features
        pct = selection_pct(m, f);
        if pct > 0
            text(f, m, sprintf('%.0f%%', pct), 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'Color', 'w', 'FontSize', 7);
        end
    end
end
hold off;

set(gca, 'Box', 'off', 'LineWidth', 1.5);

heatmap_path = fullfile(fig_dir, '01_heatmap_selection_all_methods.png');
exportgraphics(fig, heatmap_path, 'Resolution', 300);
close(fig);

fprintf('Saved: %s\n', heatmap_path);

%% ---- Figure 2: Selection frequency per method (bar plot) ----
fig = figure('Position', [100 100 1000 600], 'Color', 'w');

for m = 1:n_methods
    subplot(2, 4, m);
    pct = selection_pct(m, :);
    [pct_sorted, idx] = sort(pct, 'descend');
    labels_sorted = foot_feature_labels_short(idx);

    % Top 10 features
    n_top = min(10, n_features);
    bar(1:n_top, pct_sorted(1:n_top), 'FaceColor', [0.2 0.4 0.8]);
    set(gca, 'XTick', 1:n_top, 'XTickLabel', labels_sorted(1:n_top), ...
        'XTickLabelRotation', 45, 'FontSize', 9);
    ylabel('Selection Freq (%)', 'FontSize', 10);
    title(methods{m}, 'FontSize', 11, 'FontWeight', 'bold');
    grid on; grid minor;
    set(gca, 'Box', 'off', 'LineWidth', 1);
end

sgtitle('Top 10 Foot Features Per Method (102 runs)', ...
    'FontSize', 13, 'FontWeight', 'bold');

top10_path = fullfile(fig_dir, '02_top10_features_per_method.png');
exportgraphics(fig, top10_path, 'Resolution', 300);
close(fig);

fprintf('Saved: %s\n', top10_path);

%% ---- Figure 3: Consensus features — selected across multiple methods ----
fig = figure('Position', [100 100 1000 600], 'Color', 'w');

% Average selection frequency across methods
avg_selection = mean(selection_pct, 1);
[avg_sorted, idx] = sort(avg_selection, 'descend');
labels_sorted = foot_feature_labels_short(idx);

% Highlight by method agreement
n_show = 20;  % top 20 consensus features
colors = [];
for i = 1:n_show
    feat = idx(i);
    % Count how many methods selected this feature (threshold: >50%)
    n_methods_selected = sum(selection_pct(:, feat) > 50);
    if n_methods_selected >= 6
        colors = [colors; 0.8 0 0];  % red: 6-7 methods
    elseif n_methods_selected >= 4
        colors = [colors; 1 0.6 0];  % orange: 4-5 methods
    else
        colors = [colors; 0.2 0.6 0.8];  % blue: <4 methods
    end
end

bar(1:n_show, avg_sorted(1:n_show), 'FaceColor', 'flat', 'CData', colors);
set(gca, 'XTick', 1:n_show, 'XTickLabel', labels_sorted(1:n_show), ...
    'XTickLabelRotation', 45, 'FontSize', 10);
ylabel('Mean Selection Frequency (%) Across Methods', 'FontSize', 11);
title('Consensus Foot Features: Selected Across Multiple Methods', ...
    'FontSize', 13, 'FontWeight', 'bold');
grid on; grid minor;
set(gca, 'Box', 'off', 'LineWidth', 1.2);

% Legend
legend_handle = patch([0 0 0], [0 0 0], 'red', 'FaceAlpha', 0.7);
hold on;
patch([0 0 0], [0 0 0], [1 0.6 0], 'FaceAlpha', 0.7);
patch([0 0 0], [0 0 0], [0.2 0.6 0.8], 'FaceAlpha', 0.7);
legend('6-7 methods', '4-5 methods', '<4 methods', 'Location', 'northeast');
hold off;

consensus_path = fullfile(fig_dir, '03_consensus_features.png');
exportgraphics(fig, consensus_path, 'Resolution', 300);
close(fig);

fprintf('Saved: %s\n', consensus_path);

%% ---- Summary table ----
summary_table = table();
summary_table.Feature_Long = foot_feature_labels_long';
summary_table.Feature_Short = foot_feature_labels_short';
summary_table.MeanSelection = mean(selection_pct, 1)';
for m = 1:n_methods
    summary_table.(methods{m}) = selection_pct(m, :)';
end

% Sort by mean selection
[~, sort_idx] = sort(summary_table.MeanSelection, 'descend');
summary_table = summary_table(sort_idx, :);

% Save to CSV
csv_path = fullfile(fig_dir, 'foot_feature_selection_summary.csv');
writetable(summary_table, csv_path);

% Also save as .mat
mat_path = fullfile(fig_dir, 'foot_feature_selection_summary.mat');
save(mat_path, 'summary_table', 'selection_pct', 'foot_feature_labels_long', ...
    'foot_feature_labels_short', 'methods');

fprintf('Saved: %s\n', csv_path);
fprintf('Saved: %s\n', mat_path);

%% ---- Console report ----
fprintf('\n========================================\n');
fprintf('  GA FOOT FEATURE SELECTION ANALYSIS\n');
fprintf('========================================\n\n');

fprintf('Total LOOCV runs analysed: %d\n', n_results);
fprintf('Methods: %s\n', strjoin(methods, ', '));
fprintf('Foot features: %d (4 conditions × 5 windows × 2 metrics)\n\n', n_features);

fprintf('--- TOP 15 CONSENSUS FEATURES ---\n');
fprintf('(Selected across all methods on average)\n\n');
for i = 1:min(15, n_features)
    feat_idx = idx(i);
    feat_name = foot_feature_labels_short{feat_idx};
    avg_pct = avg_sorted(i);

    % Count method agreement
    n_methods_high = sum(selection_pct(:, feat_idx) > 50);
    n_methods_any = sum(selection_pct(:, feat_idx) > 0);

    fprintf('[%2d] %35s: %.1f%% avg  |  ', i, feat_name, avg_pct);
    fprintf('%d methods (>50%%), %d methods (any)\n', n_methods_high, n_methods_any);
end

fprintf('\n--- TOP 5 FEATURES BY METHOD ---\n\n');
for m = 1:n_methods
    method = methods{m};
    pct = selection_pct(m, :);

    % Top 5 features for this method
    [~, method_idx] = sort(pct, 'descend');
    fprintf('%10s: ', upper(method));
    for k = 1:min(5, n_features)
        f = method_idx(k);
        fprintf('%s (%.0f%%) | ', foot_feature_labels_short{f}, pct(f));
    end
    fprintf('\n');
end

fprintf('\n========================================\n');
fprintf('  Output folder: %s\n', fig_dir);
fprintf('========================================\n');
