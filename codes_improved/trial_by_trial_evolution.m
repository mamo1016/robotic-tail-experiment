% trial_by_trial_evolution.m
% Calculates trial-by-trial evolution of the foot CoP response to unexpected trials

clear; close all; clc;
addpath('functions\for_obtainERD');
config = configure_parameters();
out_dir = config.out_dir;

% Load filtered foot data
load(fullfile(out_dir, 'foot', 'foot_cleaned_session_combined_filtered.mat'), 'foot_filtered');
% Exclude participant 26 so this analysis uses the SAME N=25 cohort as the EEG
% analyses (participant 26 was a late addition recorded after the EEG cohort was
% locked; see TASKS.md). This guard mirrors trial_by_trial_all_features.m and was
% missing here, so temporal_evolution_data.mat -- and therefore the Bayesian
% equivalence analysis built on it -- ran on 26 participants while the paper
% reported N=25. Added 2026-07-25.
if isfield(foot_filtered, 'subject26')
    foot_filtered = rmfield(foot_filtered, 'subject26');
end
sub_names = fieldnames(foot_filtered);
num_subs = length(sub_names);
fprintf('Analysing %d participants (participant 26 excluded for cohort consistency).\n', num_subs);

% Time window for RMS: 0 to 1.5s post-trigger
% Assuming 2400 Hz and trigger at 2.0s (-2s to +4s epoch)
fs = config.eeg_fs;
idx_start = round(abs(config.epoch_start) * fs); % t=0
idx_end   = idx_start + round(1.5 * fs);         % t=1.5

% Find max number of unexpected trials
max_trials = 0;
for i = 1:num_subs
    n_trials = size(foot_filtered.(sub_names{i}).epochs_unexpected, 1);
    if n_trials > max_trials
        max_trials = n_trials;
    end
end

% Preallocate with NaNs
all_rms = nan(num_subs, max_trials);

for i = 1:num_subs
    sub_data = foot_filtered.(sub_names{i}).epochs_unexpected;
    n_trials = size(sub_data, 1);
    
    for t = 1:n_trials
        trial_signal = sub_data(t, idx_start:idx_end);
        all_rms(i, t) = rms(trial_signal);
    end
end

% Compute Grand Average (ignoring NaNs for subjects with fewer trials)
grand_avg = mean(all_rms, 1, 'omitnan');
sem = std(all_rms, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_rms), 1));

% Save intermediate data for MC-3 Bayesian equivalence testing
cop_rms_per_trial = all_rms;
save(fullfile(out_dir, 'temporal_evolution_data.mat'), 'cop_rms_per_trial');


% --- Plotting ---
figure('Position', [100 100 700 500], 'Color', 'w');
hold on;

% 1. Plot individual subjects in light gray
for i = 1:num_subs
    plot(1:max_trials, all_rms(i, :), 'Color', [0.8 0.8 0.8 0.5], 'LineWidth', 1);
end

% 2. Plot Grand Average
errorbar(1:max_trials, grand_avg, sem, 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 8);

% 3. Fit a linear regression to the grand average to get stats
valid_idx = find(~isnan(grand_avg));
mdl = fitlm(valid_idx', grand_avg(valid_idx)');
m = mdl.Coefficients.Estimate(2); % slope
p = mdl.Coefficients.pValue(2);   % p-value

% Overlay trendline
x_val = 1:max_trials;
y_fit = mdl.Coefficients.Estimate(1) + m * x_val;
plot(x_val, y_fit, 'r--', 'LineWidth', 2);

xlabel('Unexpected Trial Sequence Number', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Foot CoP RMS (0 to 1.5s)', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Trial-by-Trial Evolution of Postural Response\nSlope = %.4f, p = %.3f', m, p), 'FontSize', 16);
box off;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
xlim([0.5 max_trials+0.5]);
ylim([0 50]);  % Clip outlier individual traces; grand average sits ~10-20

% Save Figure to output and to paper figures folder
fig_dir = fullfile(out_dir, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
exportgraphics(gcf, fullfile(fig_dir, 'temporal_evolution.png'), 'Resolution', 300);

paper_fig_path = fullfile(fileparts(pwd), 'paper', 'figures', 'temporal_evolution.png');
exportgraphics(gcf, paper_fig_path, 'Resolution', 300);

fprintf('--- Temporal Evolution Analysis ---\n');
fprintf('Slope: %.5f\n', m);
fprintf('p-value: %.5f\n', p);
if p < 0.05
    fprintf('Result: SIGNIFICANT learning trend over time.\n');
else
    fprintf('Result: No significant linear trend.\n');
end
fprintf('Saved figure to: %s\n', fullfile(fig_dir, 'temporal_evolution.png'));
