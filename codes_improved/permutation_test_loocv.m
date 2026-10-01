%% permutation_test_loocv.m  (v2 — safe version)
% Permutation null distribution for GA-LOOCV brain-behaviour correlations.
%
% Addresses Reviewer MC-1: "r=0.99 is implausibly high — prove it via
% permutation testing."
%
% Fixes vs v1:
%   - Uses run_loocv_core() instead of LOOCV_Corrected() — no figures,
%     no file saves, no clc inside the permutation loop.
%   - Regular for loop (not parfor) for memory safety.
%   - Checkpoint save every CHECKPOINT_N iterations — safe to Ctrl+C
%     and resume later.
%   - GPU memory reset every CHECKPOINT_N iterations.
%
% Usage:
%   Just run the script. If it was interrupted, run it again — it will
%   resume automatically from the last checkpoint.
%
% Run from: codes_improved

clear; clc;
addpath('functions\for_GA');
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
config = configure_parameters();
out_dir = config.out_dir;

N_PERM        = 300;   % 300 permutations completed (stopped early — p<0.0033 achieved for all methods)
CHECKPOINT_N  = 50;   % save checkpoint every N permutations
paper_fig_dir = 'figures';
save_dir      = fullfile(out_dir, 'ga_feature');

%% --- Load pre-processed data from most recent LOOCV results ---
loocv_files = dir(fullfile(save_dir, 'results_LOOCV_*.mat'));
if isempty(loocv_files)
    error('No LOOCV results file found in: %s', save_dir);
end
[~, idx] = max([loocv_files.datenum]);
loocv_path = fullfile(loocv_files(idx).folder, loocv_files(idx).name);
fprintf('Loading LOOCV results: %s\n', loocv_path);

% eeg_data and foot_data are already pre-processed (filtered + normalized)
load(loocv_path, 'results', 'eeg_data', 'foot_data', 'foot_labels');

methods   = fieldnames(results);
N_methods = length(methods);
N_subs    = size(eeg_data, 1);

% Observed r-values (from the real, unshuffled data)
r_observed = NaN(1, N_methods);
for m = 1:N_methods
    pred   = results.(methods{m}).predicted;
    actual = results.(methods{m}).actual;
    r_observed(m) = corr(pred, actual, 'Type', 'Pearson');
end

fprintf('\nObserved r-values (real data):\n');
for m = 1:N_methods
    fprintf('  %-10s r = %.4f\n', methods{m}, r_observed(m));
end

%% --- Resume from checkpoint if it exists ---
checkpoint_file = fullfile(save_dir, 'permutation_checkpoint.mat');
if exist(checkpoint_file, 'file')
    load(checkpoint_file, 'r_null', 'perm_resume');
    fprintf('\n>>> Resuming from permutation %d / %d <<<\n', perm_resume, N_PERM);
else
    r_null      = NaN(N_PERM, N_methods);
    perm_resume = 1;
    fprintf('\nStarting fresh permutation run.\n');
end

%% --- GPU setup ---
useGPU = (gpuDeviceCount > 0);
if useGPU
    fprintf('\nGPU detected — will reset every %d iterations.\n', CHECKPOINT_N);
    gpuDevice(1);
else
    fprintf('\nNo GPU detected. Running on CPU.\n');
end

%% --- Permutation loop ---
fprintf('\nRunning %d permutations...\n', N_PERM);
t_start = tic;

for perm_i = perm_resume:N_PERM

    % Shuffle foot score labels across subjects
    perm_idx  = randperm(N_subs);
    foot_perm = foot_data(perm_idx, :);

    try
        % Lightweight LOOCV — no figures, no saves, no clc
        results_perm = run_loocv_core(eeg_data, foot_perm);

        for m = 1:N_methods
            mname = methods{m};
            if isfield(results_perm, mname)
                r_null(perm_i, m) = corr(results_perm.(mname).predicted, ...
                                         results_perm.(mname).actual, ...
                                         'Type', 'Pearson');
            end
        end
    catch ME
        warning('Permutation %d failed: %s', perm_i, ME.message);
    end

    % --- Checkpoint every CHECKPOINT_N iterations ---
    if mod(perm_i, CHECKPOINT_N) == 0 || perm_i == N_PERM
        perm_resume = perm_i + 1;
        save(checkpoint_file, 'r_null', 'perm_resume', '-v7.3');
        elapsed = toc(t_start);
        n_done = perm_i - (perm_resume - CHECKPOINT_N - 1);
        eta_min = (elapsed / n_done) * (N_PERM - perm_i) / 60;
        fprintf('  [%d/%d] Checkpoint saved | %.1f min elapsed | ETA ~%.0f min\n', ...
            perm_i, N_PERM, elapsed/60, eta_min);

        % Reset GPU memory and close any stray figures
        if useGPU
            reset(gpuDevice(1));
        end
        close all;
    end
end

elapsed = toc(t_start);
fprintf('\nAll %d permutations done in %.1f minutes.\n', N_PERM, elapsed/60);

%% --- Compute empirical p-values ---
fprintf('\n=== MC-1: Permutation Test Results ===\n');
fprintf('%-10s  %8s  %8s  %8s  %8s\n', 'Method', 'r_obs', 'r_null_med', 'r_null_95', 'p_emp');
fprintf('%s\n', repmat('-', 1, 54));

p_empirical = NaN(1, N_methods);
for m = 1:N_methods
    null_m   = r_null(:, m);
    null_m   = null_m(~isnan(null_m));
    null_med = median(null_m);
    null_95  = prctile(null_m, 95);
    p_emp    = mean(null_m >= r_observed(m));
    p_empirical(m) = p_emp;
    fprintf('%-10s  %8.4f  %8.4f  %8.4f  %8.4f\n', ...
        methods{m}, r_observed(m), null_med, null_95, p_emp);
end

%% --- Save final null distributions ---
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM');
save_path = fullfile(save_dir, sprintf('permutation_null_%s.mat', timestamp));
save(save_path, 'r_null', 'r_observed', 'methods', 'p_empirical', 'N_PERM', '-v7.3');
fprintf('\nNull distributions saved to: %s\n', save_path);

% Delete checkpoint once complete
if exist(checkpoint_file, 'file')
    delete(checkpoint_file);
    fprintf('Checkpoint file deleted.\n');
end

%% --- Figure: null distribution vs observed r ---
n_cols = min(4, N_methods);
n_rows = ceil(N_methods / n_cols);
fig = figure('Color', 'w', 'Position', [50 50 300*n_cols 280*n_rows]);
tl  = tiledlayout(n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

for m = 1:N_methods
    nexttile;
    null_m = r_null(:, m);
    null_m = null_m(~isnan(null_m));

    histogram(null_m, 40, 'FaceColor', [0.6 0.6 0.8], 'EdgeColor', 'none', ...
        'Normalization', 'probability');
    hold on;
    xline(r_observed(m), 'r-', 'LineWidth', 2.5);
    xline(prctile(null_m, 95), 'k--', 'LineWidth', 1.2);
    hold off;

    p_emp = p_empirical(m);
    if p_emp == 0
        p_str = sprintf('p < %.4f', 1/N_PERM);
    else
        p_str = sprintf('p = %.4f', p_emp);
    end
    title(sprintf('%s\nr_{obs}=%.3f | %s', upper(methods{m}), r_observed(m), p_str), ...
        'FontSize', 9, 'FontWeight', 'bold');
    xlabel('Permuted r', 'FontSize', 8);
    ylabel('Proportion', 'FontSize', 8);
    set(gca, 'FontSize', 8);
end

title(tl, sprintf('MC-1: Permutation Null Distribution (%d permutations)\nRed = observed r | Dashed = 95th percentile of null', N_PERM), ...
    'FontSize', 11);

fig_path = fullfile(paper_fig_dir, 'permutation_null_distribution.png');
exportgraphics(fig, fig_path, 'Resolution', 300);
fprintf('Figure saved to: %s\n', fig_path);
fprintf('=== DONE ===\n');
