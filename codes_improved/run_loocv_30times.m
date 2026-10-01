%% Run LOOCV 30 Times with SVR Included + Consistency Analysis
% Performs 30 independent LOOCV runs with all 7 methods (including SVR)
% Saves results to dedicated folder
% Then runs consistency analysis across all 30 runs

clear; close all; clc;

%% Configuration
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
addpath('functions\for_GA');

config = configure_parameters();
out_dir = config.out_dir;

% Create dedicated folders for 30-run results
save_dir = fullfile(out_dir, 'ga_feature_30runs');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
    fprintf('Created: %s\n', save_dir);
end

pic_out_dir = fullfile(out_dir, 'all_in_one_pic', 'consistency_30runs');
if ~exist(pic_out_dir, 'dir')
    mkdir(pic_out_dir);
    fprintf('Created: %s\n', pic_out_dir);
end

% Generate a unique batch ID for this run session
% All .mat files and output plots from this execution will share this ID
batch_id = sprintf('%04x', randi([0 65535]));

fprintf('\n=== LOOCV 30-Run Experiment with SVR ===\n');
fprintf('Batch ID: %s\n', batch_id);
fprintf('Save directory: %s\n', save_dir);
fprintf('Output directory: %s\n\n', pic_out_dir);

%% Load Prerequisites (same as main.m Step 6)
fprintf('Loading prerequisite data...\n');

% Load tile data
load(fullfile(out_dir, 'ga_feature', 'sub_tile.mat'), 'sub_tile');
load(fullfile(out_dir, 'ga_feature', 'sub_tile_foot.mat'), 'sub_tile_foot');

% Construct all_sub_diffs (copy from main.m lines 674-693)
fprintf('Constructing all_sub_diffs with EEG and foot contrasts...\n');
types = fieldnames(sub_tile.sub_1);
type_name_order = [2 1 3 2 5 4];  % Mapping for 3 contrasts: [A1 B1 A2 B2 A3 B3]
sub_list_fields = fieldnames(sub_tile);
all_sub_diffs = struct();

for idx_sub = 1:length(sub_list_fields)
    sub_name = sub_list_fields{idx_sub};
    for j = 1:3  % 3 contrasts
        eeg_type_A = types{type_name_order(2*j - 1)};
        eeg_type_B = types{type_name_order(2*j)};

        % Store tile metadata on first iteration
        if j == 1
            all_sub_diffs.(sub_name).t = sub_tile.(sub_name).(eeg_type_A).t;
            all_sub_diffs.(sub_name).f = sub_tile.(sub_name).(eeg_type_A).f;
        end

        % Compute and store EEG contrast differences
        contrast_field = [eeg_type_A '_vs_' eeg_type_B];
        all_sub_diffs.(sub_name).eeg.(contrast_field) = ...
            sub_tile.(sub_name).(eeg_type_A).vector - sub_tile.(sub_name).(eeg_type_B).vector;

        % Compute and store foot contrast differences
        foot_type_A = name_finder(eeg_type_A);
        foot_type_B = name_finder(eeg_type_B);
        all_sub_diffs.(sub_name).foot.(contrast_field) = ...
            sub_tile_foot.(sub_name).(foot_type_A).vector - sub_tile_foot.(sub_name).(foot_type_B).vector;
    end
end

fprintf('all_sub_diffs constructed: %d subjects, with EEG/foot contrasts\n', length(fieldnames(all_sub_diffs)) - 2);

% --- DIAGNOSTICS: verify all_sub_diffs contents before LOOCV ---
fprintf('\n[DIAG] === all_sub_diffs verification ===\n');
fprintf('[DIAG] Subject list (first 3): %s\n', strjoin(sub_list_fields(1:min(3,end))', ', '));
eeg_contrast_fields = fieldnames(all_sub_diffs.sub_1.eeg);
fprintf('[DIAG] EEG contrast fields: %s\n', strjoin(eeg_contrast_fields', ', '));
foot_contrast_fields = fieldnames(all_sub_diffs.sub_1.foot);
fprintf('[DIAG] Foot contrast fields: %s\n', strjoin(foot_contrast_fields', ', '));
for cf = 1:length(eeg_contrast_fields)
    fn = eeg_contrast_fields{cf};
    v = all_sub_diffs.sub_1.eeg.(fn);
    fprintf('[DIAG] eeg.%s: size=[%s] anyNaN=%d allZero=%d\n', fn, num2str(size(v)), any(isnan(v(:))), all(v(:)==0));
end

% Load foot data ONCE (before loop) to avoid redundant disk reads
fprintf('\nLoading foot data...\n');
load(fullfile(out_dir, 'foot', 'foot_cleaned_session_combined_filtered.mat'), 'foot_filtered');
[foot_features, foot_labels] = generate_foot_features(foot_filtered, config);
fprintf('Foot features: %d subjects × %d features\n', size(foot_features, 1), size(foot_features, 2));
fprintf('[DIAG] foot_features: anyNaN=%d allZero=%d std=%.4f range=[%.4f %.4f]\n', ...
    any(isnan(foot_features(:))), all(foot_features(:)==0), std(foot_features(:)), min(foot_features(:)), max(foot_features(:)));

%% Run LOOCV (5 times for quick verification, exclude SVR)
% For verification: Run 5 loops with 6 methods (exclude slow SVR)
% After verification: Increase to 30 loops with all 7 methods

N_RUNS = 5;  % Quick verification (change to 30 for production)
EXCLUDE_SVR = true;  % Skip SVR for speed (set to false for full 7-method run)

methods_to_use = {'lasso', 'ridge'};
if ~EXCLUDE_SVR
    methods_to_use = {'spearman', 'pca', 'robust', 'pls', 'lasso', 'ridge', 'svr'};
end

fprintf('\n=== STARTING %d LOOCV RUNS ===\n', N_RUNS);
fprintf('Methods: %s\n', strjoin(methods_to_use, ', '));
fprintf('(SVR excluded for speed. To include SVR, set EXCLUDE_SVR = false)\n\n');

for run_idx = 1:N_RUNS
    deleteFig();
    fprintf('\n');
    fprintf('╔════════════════════════════════════╗\n');
    fprintf('║  Run %2d / %2d                       ║\n', run_idx, N_RUNS);
    fprintf('╚════════════════════════════════════╝\n');

    try
        % Run LOOCV with custom save directory, methods parameter, and batch ID
        [~, save_path] = LOOCV_Corrected(all_sub_diffs, foot_features, foot_labels, config, save_dir, methods_to_use, batch_id);

        fprintf('✓ Run %d complete: %s\n', run_idx, save_path);

        % Small pause to avoid file conflicts
        pause(0.5);

    catch ME
        fprintf('✗ Run %d FAILED: %s\n', run_idx, ME.message);
        continue;
    end
end

fprintf('\n╔════════════════════════════════════╗\n');
fprintf('║  All %d LOOCV runs complete!       ║\n', N_RUNS);
fprintf('╚════════════════════════════════════╝\n');

%% Post-processing: Aggregate figures from ALL saved .mat files
% Analyzes every results_LOOCV_*.mat in save_dir — no batch filter.
% Outputs: (1) LOOCV accuracy box plot, (2) tile consistency maps.
%% Configuration
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
addpath('functions\for_GA');
config = configure_parameters();
out_dir = config.out_dir;
save_dir     = fullfile(out_dir, 'ga_feature_30runs');
results_dir  = fullfile(save_dir, 'toResult');  % 100-run archive lives here
pic_out_dir  = fullfile(out_dir, 'all_in_one_pic', 'consistency_30runs');
if ~exist(pic_out_dir, 'dir'), mkdir(pic_out_dir); end

all_mat_files = dir(fullfile(results_dir, 'results_LOOCV_*.mat'));
fprintf('\nFound %d result files in %s\n', length(all_mat_files), results_dir);
if isempty(all_mat_files)
    error('No result files found in %s', results_dir);
end

%% Figure 1: LOOCV Prediction Accuracy (Spearman r) across all runs
method_names = {'spearman', 'pca', 'robust', 'pls', 'lasso', 'ridge', 'svr'};
r_all = struct();
for m = 1:length(method_names)
    r_all.(method_names{m}) = [];
end

for f_idx = 1:length(all_mat_files)
    fpath = fullfile(results_dir, all_mat_files(f_idx).name);
    try
        d = load(fpath, 'results');
        for m = 1:length(method_names)
            method = method_names{m};
            if isfield(d.results, method)
                pred   = d.results.(method).predicted;
                actual = d.results.(method).actual;
                if std(pred) > 1e-5 && std(actual) > 1e-5
                    r_val = corr(pred, actual, 'Type', 'Spearman');
                    r_all.(method) = [r_all.(method); r_val];
                end
            end
        end
    catch ME
        fprintf('  Skipped %s: %s\n', all_mat_files(f_idx).name, ME.message);
    end
end

% Print summary table
fprintf('\n%-12s  %5s  %7s  %7s  %13s\n', 'Method', 'N', 'Mean r', 'Std r', 'Range');
fprintf('%s\n', repmat('-', 1, 52));
for m = 1:length(method_names)
    method = method_names{m};
    vals = r_all.(method);
    if ~isempty(vals)
        fprintf('%-12s  %5d  %7.3f  %7.3f  [%5.3f %5.3f]\n', ...
            method, length(vals), mean(vals), std(vals), min(vals), max(vals));
    else
        fprintf('%-12s  %5s  %7s  %7s  %13s\n', method, '-', '-', '-', 'no data');
    end
end

% Build box plot
methods_with_data = {};
data_matrix_cols  = {};
for m = 1:length(method_names)
    vals = r_all.(method_names{m});
    if ~isempty(vals)
        methods_with_data{end+1} = upper(method_names{m});
        data_matrix_cols{end+1}  = vals;
    end
end

if ~isempty(data_matrix_cols)
    max_n = max(cellfun(@length, data_matrix_cols));
    bx_data = NaN(max_n, length(data_matrix_cols));
    for i = 1:length(data_matrix_cols)
        n = length(data_matrix_cols{i});
        bx_data(1:n, i) = data_matrix_cols{i};
    end

    fig_acc = figure('Color', 'w', 'Position', [100 100 950 480]);
    boxplot(bx_data, 'Labels', methods_with_data, 'Widths', 0.5);
    hold on;
    yline(0, 'k--', 'LineWidth', 0.8);
    hold off;
    ylabel('Spearman r (held-out participants)', 'FontSize', 12);
    title(sprintf('LOOCV prediction accuracy across %d runs', length(all_mat_files)), ...
          'FontSize', 13, 'FontWeight', 'bold');
    ylim([-0.3 1.1]);
    grid on;
    set(gca, 'FontSize', 11);

    acc_path = fullfile(pic_out_dir, 'loocv_accuracy_summary.png');
    exportgraphics(fig_acc, acc_path, 'Resolution', 300);
    fprintf('\nSaved: %s\n', acc_path);
    close(fig_acc);
end

%% Figure 2: Tile Consistency Maps (all runs, no batch filter)
fprintf('\n=== RUNNING CONSISTENCY ANALYSIS (all runs) ===\n');
analyze_loocv_consistency(results_dir, pic_out_dir, ...
    fullfile(out_dir, 'all_in_one_pic'), fullfile(out_dir, 'ga_feature'));

fprintf('\n✓ Post-processing complete.\n');
fprintf('  Files analysed: %d\n', length(all_mat_files));
fprintf('  Plots saved to: %s\n', pic_out_dir);
