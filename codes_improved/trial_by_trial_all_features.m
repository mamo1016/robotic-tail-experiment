% trial_by_trial_all_features.m
% Trial-by-trial trend analysis for ALL 40 foot CoP features.
%
% Problem: each participant has a different number of trials per condition.
% Solution: percentile binning — divide each subject's trial sequence into
% K equal-width bins (deciles by default), compute the metric within each
% bin, then compare across subjects on the common bin axis.
%
% For each of the 40 features (4 conditions x 5 time windows x 2 metrics):
%   1. Extract per-trial values for every subject
%   2. Bin into K percentile bins (each subject -> K values)
%   3. Grand mean +/- bootstrap 95% CI across subjects
%   4. Linear regression on binned grand mean; report slope & p
%   5. Save figure to output folder
%
% Output folder: <out_dir>/trial_evolution_40features/

clear; close all; clc;
addpath('functions\for_obtainERD');
config = configure_parameters();
out_dir = config.out_dir;

%% ---- Settings ----
K = 10;             % number of percentile bins (deciles)
n_boot = 2000;      % bootstrap iterations for CI
ci_alpha = 0.05;    % 95% confidence interval

%% ---- Load foot data ----
load(fullfile(out_dir, 'foot', 'foot_cleaned_session_combined_filtered.mat'), ...
     'foot_filtered');
% Exclude participant 26 so this analysis uses the SAME N=25 cohort as the EEG
% analyses (participant 26 was a late addition recorded after the EEG cohort was
% locked; see TASKS.md). Keeps every analysis on one consistent 25-participant set.
if isfield(foot_filtered, 'subject26')
    foot_filtered = rmfield(foot_filtered, 'subject26');
end
sub_names = fieldnames(foot_filtered);
num_subs  = length(sub_names);
fprintf('Analysing %d participants (participant 26 excluded for cohort consistency).\n', num_subs);

%% ---- Define the 40 features (matches generate_foot_features.m) ----
conditions   = {'epochs_expected', 'epochs_unexpected', ...
                'epochs_right_before', 'epochs_right_after'};
cond_labels  = {'Expected', 'Unexpected', 'RightBefore', 'RightAfter'};

time_windows = {[-0.5, 0], [0, 0.5], [0.5, 1.0], [1.0, 1.5], [0, 1.5]};
win_labels   = {'Baseline(-0.5to0s)', 'Reflex(0to0.5s)', ...
                'EarlyComp(0.5to1s)', 'LateComp(1to1.5s)', ...
                'TotalPost(0to1.5s)'};

metric_names = {'RMS', 'Area'};

% Display labels used ONLY in figure titles. Deliberately kept separate from
% cond_labels/win_labels above, because those two build feat_name -> safe_name
% -> the output FILENAME, and main.tex / supplementary.tex reference those
% filenames directly (e.g. 07_Expected_LateComp1to15s_RMS.png). Renaming them
% would silently break every \includegraphics path.
cond_display = {'Expected', 'Unexpected', 'Right Before', 'Right After'};
win_display  = {'baseline window (-0.5 to 0 s)', 'reflex window (0-0.5 s)', ...
                'early compensation (0.5-1 s)', 'late compensation (1-1.5 s)', ...
                'total post-event (0-1.5 s)'};

% Time vector (foot data resampled to EEG sampling rate)
fs    = config.eeg_fs;   % 1200 Hz
times = config.epoch_start : (1/fs) : config.epoch_end;

%% ---- Create output folder ----
fig_dir = fullfile(out_dir, 'trial_evolution_40features');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% ---- Preallocate summary table ----
n_features = length(conditions) * length(time_windows) * length(metric_names);
summary = struct('feature', {}, 'slope', {}, 'p_value', {}, ...
                 'r_squared', {}, 'direction', {});

feat_idx = 0;

%% ---- Main loop: iterate over all 40 features ----
for c = 1:length(conditions)
    cond_name  = conditions{c};
    cond_label = cond_labels{c};
    cond_disp  = cond_display{c};

    for w = 1:length(time_windows)
        win       = time_windows{w};
        win_label = win_labels{w};
        win_disp  = win_display{w};

        % Find time indices for this window
        t_idx = (times >= win(1)) & (times <= win(2));
        if sum(t_idx) == 0
            warning('No samples in window %s — skipping.', win_label);
            continue;
        end

        for m = 1:length(metric_names)
            metric = metric_names{m};
            feat_idx = feat_idx + 1;
            feat_name = sprintf('%s_%s_%s', cond_label, win_label, metric);

            fprintf('[%2d/40] %s ...', feat_idx, feat_name);

            % ---- Step 1: extract per-trial values for every subject ----
            trial_vals = cell(num_subs, 1);   % each cell: [n_trials x 1]
            trial_counts = zeros(num_subs, 1);

            for s = 1:num_subs
                sub = foot_filtered.(sub_names{s});

                if ~isfield(sub, cond_name)
                    trial_vals{s} = [];
                    continue;
                end

                data = sub.(cond_name);  % [n_trials x n_timepoints]
                n_trials = size(data, 1);
                win_data = data(:, t_idx);  % [n_trials x win_samples]

                vals = zeros(n_trials, 1);
                for tr = 1:n_trials
                    seg = win_data(tr, :);
                    switch metric
                        case 'RMS'
                            vals(tr) = sqrt(mean(seg .^ 2));
                        case 'Area'
                            vals(tr) = trapz(abs(seg));
                    end
                end

                trial_vals{s}   = vals;
                trial_counts(s) = n_trials;
            end

            % Skip if no subject has data for this condition
            valid_subs = find(trial_counts > 0);
            if isempty(valid_subs)
                fprintf(' no data — skipped.\n');
                continue;
            end

            % ---- Step 2: percentile binning ----
            % Each subject's trials -> K bins. Within each bin, take the
            % mean. Result: [num_valid_subs x K] matrix.
            binned = nan(num_subs, K);

            for s = valid_subs'
                v = trial_vals{s};
                n = length(v);
                % Assign each trial to a bin (1..K)
                bin_edges = linspace(0, n, K + 1);
                for b = 1:K
                    lo = floor(bin_edges(b)) + 1;
                    hi = floor(bin_edges(b + 1));
                    if hi < lo, hi = lo; end  % single-trial bin
                    binned(s, b) = mean(v(lo:hi));
                end
            end

            % Keep only valid subjects
            binned_valid = binned(valid_subs, :);
            n_valid = size(binned_valid, 1);

            % ---- Step 3: grand mean + bootstrap 95% CI ----
            grand_mean = mean(binned_valid, 1);

            boot_means = zeros(n_boot, K);
            for iter = 1:n_boot
                idx = randi(n_valid, n_valid, 1);  % resample subjects
                boot_means(iter, :) = mean(binned_valid(idx, :), 1);
            end
            ci_lo = prctile(boot_means, 100 * ci_alpha / 2, 1);
            ci_hi = prctile(boot_means, 100 * (1 - ci_alpha / 2), 1);

            % ---- Step 4: linear regression on grand mean ----
            x = (1:K)';
            mdl = fitlm(x, grand_mean');
            slope   = mdl.Coefficients.Estimate(2);
            p_val   = mdl.Coefficients.pValue(2);
            r2      = mdl.Rsquared.Ordinary;

            if slope > 0
                direction = 'increasing';
            else
                direction = 'decreasing';
            end

            % Also bootstrap the slope for a CI
            boot_slopes = zeros(n_boot, 1);
            for iter = 1:n_boot
                idx = randi(n_valid, n_valid, 1);
                bm  = mean(binned_valid(idx, :), 1);
                b_mdl = fitlm(x, bm');
                boot_slopes(iter) = b_mdl.Coefficients.Estimate(2);
            end
            slope_ci_lo = prctile(boot_slopes, 100 * ci_alpha / 2);
            slope_ci_hi = prctile(boot_slopes, 100 * (1 - ci_alpha / 2));

            % Store summary
            summary(feat_idx).feature   = feat_name;
            summary(feat_idx).slope     = slope;
            summary(feat_idx).p_value   = p_val;
            summary(feat_idx).r_squared = r2;
            summary(feat_idx).direction = direction;
            summary(feat_idx).slope_ci  = [slope_ci_lo, slope_ci_hi];
            summary(feat_idx).n_subs    = n_valid;
            summary(feat_idx).trial_counts = trial_counts(valid_subs)';

            % ---- Step 5: plot and save (two versions) ----
            safe_name = strrep(feat_name, '(', '');
            safe_name = strrep(safe_name, ')', '');
            safe_name = strrep(safe_name, '.', '');

            % ---- VERSION 1: WITH individual subject traces (original) ----
            fig = figure('Position', [100 100 800 500], 'Color', 'w', ...
                         'Visible', 'off');

            % Shaded bootstrap CI
            fill([1:K, K:-1:1], [ci_hi, fliplr(ci_lo)], ...
                 [0.7 0.85 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
            hold on;

            % Individual subject traces (light grey)
            for s = 1:n_valid
                plot(1:K, binned_valid(s, :), '-', ...
                     'Color', [0.8 0.8 0.8 0.4], 'LineWidth', 0.8);
            end

            % Grand mean
            plot(1:K, grand_mean, 'ko-', 'LineWidth', 2, ...
                 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerSize', 7);

            % Regression line
            y_fit = mdl.Coefficients.Estimate(1) + slope * x;
            plot(x, y_fit, 'r--', 'LineWidth', 2);

            % Significance marker
            if p_val < 0.05
                sig_str = sprintf('slope=%.4f, p=%.4f *', slope, p_val);
                sig_col = [0.8 0 0];
            else
                sig_str = sprintf('slope=%.4f, p=%.3f (n.s.)', slope, p_val);
                sig_col = [0.3 0.3 0.3];
            end

            xlabel('Trial Bin (decile)', 'FontSize', 13, 'FontWeight', 'bold');
            ylabel(sprintf('CoP %s', metric), 'FontSize', 13, 'FontWeight', 'bold');
            title(sprintf('%s trials: %s (%s)', cond_disp, win_disp, metric), ...
                  'FontSize', 14, 'FontWeight', 'bold');
            subtitle(sig_str, 'FontSize', 11, 'Color', sig_col);

            % Axis labels: decile percentages
            xticks(1:K);
            pct_labels = cell(1, K);
            for b = 1:K
                pct_labels{b} = sprintf('%d%%', b * (100 / K));
            end
            xticklabels(pct_labels);
            xlim([0.5, K + 0.5]);

            set(gca, 'FontSize', 11, 'LineWidth', 1.2, 'Box', 'off');

            % Trial count annotation
            min_t = min(trial_counts(valid_subs));
            max_t = max(trial_counts(valid_subs));
            annotation_str = sprintf('N=%d participants, trials/participant: %d–%d', ...
                                     n_valid, min_t, max_t);
            text(0.98, 0.02, annotation_str, 'Units', 'normalized', ...
                 'HorizontalAlignment', 'right', 'FontSize', 9, ...
                 'Color', [0.4 0.4 0.4]);

            hold off;

            % Save version 1
            fig_path_v1 = fullfile(fig_dir, sprintf('%02d_%s_with_subjects.png', feat_idx, safe_name));
            exportgraphics(fig, fig_path_v1, 'Resolution', 300);
            close(fig);

            % ---- VERSION 2: CLEAN (without individual subject traces) ----
            fig = figure('Position', [100 100 800 500], 'Color', 'w', ...
                         'Visible', 'off');

            % Shaded bootstrap CI
            fill([1:K, K:-1:1], [ci_hi, fliplr(ci_lo)], ...
                 [0.7 0.85 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
                 'DisplayName', 'Bootstrap 95% CI');
            hold on;

            % Grand mean
            h_mean = plot(1:K, grand_mean, 'ko-', 'LineWidth', 2, ...
                 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerSize', 7, ...
                 'DisplayName', 'Grand Mean');

            % Regression line
            y_fit = mdl.Coefficients.Estimate(1) + slope * x;
            h_reg = plot(x, y_fit, 'r--', 'LineWidth', 2, ...
                 'DisplayName', 'Linear Fit');

            xlabel('Trial Bin (decile)', 'FontSize', 13, 'FontWeight', 'bold');
            ylabel(sprintf('CoP %s', metric), 'FontSize', 13, 'FontWeight', 'bold');
            title(sprintf('%s trials: %s (%s)', cond_disp, win_disp, metric), ...
                  'FontSize', 14, 'FontWeight', 'bold');
            subtitle(sig_str, 'FontSize', 11, 'Color', sig_col);

            % Axis labels: decile percentages
            xticks(1:K);
            xticklabels(pct_labels);
            xlim([0.5, K + 0.5]);

            set(gca, 'FontSize', 11, 'LineWidth', 1.2, 'Box', 'off');

            % Legend. Pinned to northwest rather than 'best': MATLAB's 'best'
            % only avoids DATA, not the trial-count text() below, and it was
            % landing bottom-right and obscuring it (e.g. 07_Expected_LateComp).
            legend('Location', 'northwest', 'FontSize', 11, 'Box', 'on', ...
                   'EdgeColor', 'k', 'LineWidth', 0.8);

            % Trial count annotation
            text(0.98, 0.02, annotation_str, 'Units', 'normalized', ...
                 'HorizontalAlignment', 'right', 'FontSize', 9, ...
                 'Color', [0.4 0.4 0.4]);

            hold off;

            % Save version 2
            fig_path_v2 = fullfile(fig_dir, sprintf('%02d_%s_clean.png', feat_idx, safe_name));
            exportgraphics(fig, fig_path_v2, 'Resolution', 300);
            close(fig);

            fprintf(' done (slope=%.4f, p=%.4f)\n', slope, p_val);
        end
    end
end

%% ---- Save summary ----
% Summary table as CSV
summary_table = struct2table(summary);
csv_path = fullfile(fig_dir, 'trend_summary_40features.csv');

% Build CSV manually (struct2table may have nested fields)
fid = fopen(csv_path, 'w');
fprintf(fid, 'feature,slope,p_value,r_squared,direction,slope_ci_lo,slope_ci_hi,n_subs\n');
for i = 1:length(summary)
    if isempty(summary(i).feature), continue; end
    fprintf(fid, '%s,%.6f,%.6f,%.4f,%s,%.6f,%.6f,%d\n', ...
        summary(i).feature, summary(i).slope, summary(i).p_value, ...
        summary(i).r_squared, summary(i).direction, ...
        summary(i).slope_ci(1), summary(i).slope_ci(2), ...
        summary(i).n_subs);
end
fclose(fid);

% Also save as .mat
mat_path = fullfile(fig_dir, 'trend_summary_40features.mat');
save(mat_path, 'summary');

fprintf('\n========================================\n');
fprintf('  TREND ANALYSIS COMPLETE\n');
fprintf('  Figures saved to: %s\n', fig_dir);
fprintf('  Summary CSV: %s\n', csv_path);
fprintf('========================================\n');

% Print significant features
fprintf('\n--- Features with significant trend (p < 0.05) ---\n');
any_sig = false;
for i = 1:length(summary)
    if isempty(summary(i).feature), continue; end
    if summary(i).p_value < 0.05
        fprintf('  [%02d] %s: slope=%.4f, p=%.4f, R2=%.3f (%s)\n', ...
            i, summary(i).feature, summary(i).slope, ...
            summary(i).p_value, summary(i).r_squared, ...
            summary(i).direction);
        any_sig = true;
    end
end
if ~any_sig
    fprintf('  None.\n');
end
