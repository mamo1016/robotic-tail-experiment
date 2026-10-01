%% compute_average_loocv_performance.m
% Compute average correlation (r) and p-values across 30 LOOCV runs
% for all 7 body prediction methods (Spearman, PCA, Robust, PLS, LASSO, Ridge, SVR)
%
% Input: All results_LOOCV_*.mat files in output\ga_feature_30runs\toResult\
% Output: Table with mean r, std r, mean p, std p for each method

clear; clc;

% Set paths
results_dir = 'output/ga_feature_30runs/toResult/';
method_names = {'spearman', 'pca', 'robust', 'pls', 'lasso', 'ridge', 'svr'};
n_methods = length(method_names);

% Find all LOOCV result files
loocv_files = dir(fullfile(results_dir, 'results_LOOCV_*.mat'));
fprintf('Found %d LOOCV result files\n', length(loocv_files));

if isempty(loocv_files)
    error('No results_LOOCV_*.mat files found in %s', results_dir);
end

% Initialize storage
results_r = zeros(length(loocv_files), n_methods);  % [n_runs × 7]
results_p = zeros(length(loocv_files), n_methods);

%% Load all files and compute r, p for each method
for f_idx = 1:length(loocv_files)
    file_path = fullfile(results_dir, loocv_files(f_idx).name);

    try
        data = load(file_path, 'results');
        results = data.results;
    catch ME
        fprintf('ERROR loading %s: %s\n', loocv_files(f_idx).name, ME.message);
        continue;
    end

    % For each method
    for m_idx = 1:n_methods
        method = method_names{m_idx};

        if isfield(results, method)
            % Extract predicted and actual values
            predicted = results.(method).predicted;  % [N × 25] or similar structure
            actual = results.(method).actual;        % [N × 25] or similar structure

            % Flatten to vectors (concatenate all subjects/samples)
            predicted_vec = predicted(:);
            actual_vec = actual(:);

            % Compute Spearman correlation (r) and p-value
            [r, p] = corr(predicted_vec, actual_vec, 'Type', 'Spearman');

            results_r(f_idx, m_idx) = r;
            results_p(f_idx, m_idx) = p;
        else
            % Method not available in this run
            results_r(f_idx, m_idx) = NaN;
            results_p(f_idx, m_idx) = NaN;
        end
    end

    fprintf('Run %d/%d: %s - processed\n', f_idx, length(loocv_files), loocv_files(f_idx).name);
end

%% Compute statistics
fprintf('\n========== LOOCV BODY PREDICTION PERFORMANCE (30 RUNS) ==========\n\n');
fprintf('%-12s | %8s ± %6s | %8s ± %8s\n', 'Method', 'Mean r', 'Std r', 'Mean p', 'Std p');
fprintf('%s\n', repmat('-', 60, 1));

mean_r = nanmean(results_r, 1);
std_r  = nanstd(results_r, 0, 1);
mean_p = nanmean(results_p, 1);
std_p  = nanstd(results_p, 0, 1);

for m_idx = 1:n_methods
    method = method_names{m_idx};
    fprintf('%-12s | %8.4f ± %6.4f | %8.6f ± %8.6f\n', ...
        method, mean_r(m_idx), std_r(m_idx), mean_p(m_idx), std_p(m_idx));
end

%% Save to file
fprintf('\n========== PER-RUN DETAILS ==========\n\n');
fprintf('%-12s | ', 'Run');
for m_idx = 1:n_methods
    fprintf('%10s | ', method_names{m_idx});
end
fprintf('\n');
fprintf('%s\n', repmat('-', 100, 1));

for f_idx = 1:length(loocv_files)
    [~, fname, ~] = fileparts(loocv_files(f_idx).name);
    fprintf('%-12s | ', fname(15:end));  % Show timestamp
    for m_idx = 1:n_methods
        fprintf('%10.4f | ', results_r(f_idx, m_idx));
    end
    fprintf('\n');
end

%% Create summary table
summary_table = table(method_names', mean_r', std_r', mean_p', std_p', ...
    'VariableNames', {'Method', 'Mean_r', 'Std_r', 'Mean_p', 'Std_p'});

fprintf('\n=== SUMMARY TABLE ===\n');
disp(summary_table);

%% Save summary to file
output_file = fullfile(results_dir, '..', 'average_loocv_performance_30runs.mat');
save(output_file, 'summary_table', 'results_r', 'results_p', 'method_names');
fprintf('\n✓ Summary saved to: %s\n', output_file);

% Also save as CSV
csv_file = fullfile(results_dir, '..', 'average_loocv_performance_30runs.csv');
writetable(summary_table, csv_file);
fprintf('✓ CSV saved to: %s\n', csv_file);
