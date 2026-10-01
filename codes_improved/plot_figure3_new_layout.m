%% plot_figure3_new_layout.m
% Regenerates Figure 3 using Option A layout (2-row x 3-column):
%
%   Row 1: ERSP difference heatmaps
%     (a) Control  (b) Surprise  (c) Learning
%   Row 2: Statistical mask + GA consensus (yellow = selected >=50% folds)
%     (d) Control  (e) Surprise  (f) Learning
%
% Data sources (read-only):
%   - avergae_bootstrap_all_sub.mat  (grand-average bootstrapped ERSP)
%   - results_LOOCV_*.mat            (GA tile selection counts)
%
% NEW STANDALONE SCRIPT — does not modify any existing files.
% Run from: codes_improved
%
% Output: figures\figure3_new_layout.png

clear; clc;
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_GA');
config = configure_parameters();

out_dir   = config.out_dir;
save_path = fullfile(fileparts(pwd), 'paper', 'figures', 'figure3_new_layout.png');

%% --- Load grand-average bootstrapped ERSP ---
boot_file = fullfile(out_dir, 'all_in_one_pic', 'avergae_bootstrap_all_sub.mat');
if ~exist(boot_file, 'file')
    error('Bootstrap average file not found:\n  %s', boot_file);
end
load(boot_file, 'ave_boot_average');

[time_range, freq_range] = range_compute(config);

% Compute the three difference maps (grand average)
ersp_control  = ave_boot_average.expected_duplicated.ersp - ave_boot_average.expected.ersp;
ersp_surprise = ave_boot_average.unexpected.ersp          - ave_boot_average.expected.ersp;
ersp_learning = ave_boot_average.right_after.ersp         - ave_boot_average.right_before.ersp;

%% --- Load LOOCV results for GA consensus ---
ga_dir = fullfile(out_dir, 'ga_feature');

% Load GA selection counts from Mar 5 19:56 (has mask_x_counts; GA ran correctly)
loocv_counts_path = fullfile(ga_dir, 'results_LOOCV_2026-03-05_19-56.mat');
if ~exist(loocv_counts_path, 'file')
    error('Counts LOOCV file not found:\n  %s', loocv_counts_path);
end
fprintf('Loading GA counts: results_LOOCV_2026-03-05_19-56.mat\n');
counts_data = load(loocv_counts_path, 'results', 'run_info');
ga_counts = counts_data.results.spearman.mask_x_counts;  % 1x1200, values 0-25
ga_mask   = counts_data.results.spearman.mask_x;         % binary consensus >=50%
n_subs_loocv = 25;

% tile_t and tile_f come from sub_tile (not run_info)
tmp_tile = load(fullfile(out_dir, 'ga_feature', 'sub_tile.mat'), 'sub_tile');
sub_tile_tmp = tmp_tile.sub_tile;
first_sub = fieldnames(sub_tile_tmp); first_sub = first_sub{1};
first_type = fieldnames(sub_tile_tmp.(first_sub)); first_type = first_type{1};
tile_t   = sub_tile_tmp.(first_sub).(first_type).t;
tile_f   = sub_tile_tmp.(first_sub).(first_type).f;
n_tile_f = length(tile_f);
n_tile_t = length(tile_t);
clear sub_tile_tmp tmp_tile

% Load per-subject raw ERSP data for correct per-phase t-test
% (eeg_data from LOOCV is Z-normalised across all 1200 tiles, which distorts
%  the per-contrast signal structure and makes the three bottom panels look
%  identical.  We instead reconstruct raw tiled difference maps per subject.)
each_sub_file = fullfile(out_dir, 'all_in_one_pic', 'avergae_bootstrap_each_sub.mat');
tmp           = load(each_sub_file, 'boot_average');
boot_each     = tmp.boot_average;
sub_names_raw = fieldnames(boot_each);
n_subs_raw    = length(sub_names_raw);

raw_control  = zeros(n_subs_raw, 400);
raw_surprise = zeros(n_subs_raw, 400);
raw_learning = zeros(n_subs_raw, 400);

for s = 1:n_subs_raw
    sn = sub_names_raw{s};
    tc = ersp_tile(boot_each.(sn).expected_duplicated.ersp - boot_each.(sn).expected.ersp,   time_range, freq_range);
    ts = ersp_tile(boot_each.(sn).unexpected.ersp          - boot_each.(sn).expected.ersp,   time_range, freq_range);
    tl = ersp_tile(boot_each.(sn).right_after.ersp         - boot_each.(sn).right_before.ersp, time_range, freq_range);
    raw_control(s, :)  = tc.vector(:)';
    raw_surprise(s, :) = ts.vector(:)';
    raw_learning(s, :) = tl.vector(:)';
end

[h_ctrl,  ~] = ttest(raw_control,  0, 'Alpha', 0.05);
[h_surp,  ~] = ttest(raw_surprise, 0, 'Alpha', 0.05);
[h_learn, ~] = ttest(raw_learning, 0, 'Alpha', 0.05);

h_sig_phases  = {h_ctrl,                  h_surp,                  h_learn};
mean_phases   = {mean(raw_control,  1),   mean(raw_surprise, 1),   mean(raw_learning, 1)};

% Phase index mapping (as in LOOCV_Plot):
%   1:400   = Expected vs Expected  (control)
%   401:800 = Unexpected vs Expected (surprise)
%   801:1200= Right After vs Right Before (learning)
phase_idx = {1:400, 401:800, 801:1200};

% GA counts and consensus mask from Mar 5 19:56 (GA selection counts are valid)
selection_counts = ga_counts;  % 1x1200, values 0-25
mask_x           = ga_mask;    % binary consensus >=50%

%% --- Build colourmaps ---
cmap_ersp = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ...
             ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];
% Black-to-yellow heatmap for selection counts (matching LOOCV_Plot style)
cmap_count = [zeros(64,1), zeros(64,1), zeros(64,1); ...        % black
              linspace(0,1,64)', zeros(64,1), zeros(64,1); ...   % black->red
              ones(64,1), linspace(0,1,64)', zeros(64,1); ...    % red->yellow
              ones(64,1), ones(64,1), linspace(0,1,64)'];        % yellow->white

%% --- Plot: 2 rows x 3 columns ---
fig = figure('Color', 'w', 'Position', [50 50 1500 700]);
tl  = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ersp_maps   = {ersp_control, ersp_surprise, ersp_learning};
row1_titles = {'(a)  Expected vs Expected  [Control baseline]', ...
               '(b)  Unexpected vs Expected  [Surprise signal]', ...
               '(c)  Right After vs Right Before  [Learning signal]'};
row2_titles = {'(d)  GA Consensus  [Control]', ...
               '(e)  GA Consensus  [Surprise]', ...
               '(f)  GA Consensus  [Learning]'};

clim_ersp = 0.1;   % dB symmetric limits for ERSP maps
clim_sig  = 0.1;   % same limits for masked GA panels

%% -- Row 1: ERSP heatmaps --
for col = 1:3
    nexttile(col);
    imagesc(time_range, freq_range, ersp_maps{col}, [-clim_ersp clim_ersp]);
    axis xy;
    colormap(gca, cmap_ersp);
    cb = colorbar; cb.Label.String = '\DeltaERSP (dB)'; cb.FontSize = 8;
    xlabel('Time (s)', 'FontSize', 10);
    ylabel('Frequency (Hz)', 'FontSize', 10);
    ylim([5 35]);
    hold on; xline(0, 'k--', 'LineWidth', 0.8); hold off;
    title(row1_titles{col}, 'FontSize', 10, 'FontWeight', 'bold');
    set(gca, 'FontSize', 9);
end

%% -- Row 2: Selection-count heatmap + GA consensus circles (>=50% LOOCV folds) --
for col = 1:3
    nexttile(col + 3);

    p_idx = phase_idx{col};

    % Per-phase selection counts reshaped to tile grid
    counts_phase = selection_counts(p_idx);
    counts_map   = reshape(counts_phase, [n_tile_f, n_tile_t]);

    % GA consensus circles from mask_x (>=50% of 25 folds)
    vis_mask = mask_x(p_idx);
    [sel_t, sel_f] = mask_visualiser(vis_mask, tile_f, tile_t);
    n_selected = sum(vis_mask);

    imagesc(tile_t, tile_f, counts_map, [0 n_subs_loocv]);
    axis xy;
    colormap(gca, cmap_count);
    cb = colorbar; cb.Label.String = 'Selection count (/ 25 folds)'; cb.FontSize = 8;
    xlabel('Time (s)', 'FontSize', 10);
    ylabel('Frequency (Hz)', 'FontSize', 10);
    ylim([5 35]);
    hold on;
    xline(0, 'k--', 'LineWidth', 0.8, 'Color', [0.7 0.7 0.7]);
    if ~isempty(sel_t)
        plot(sel_t, sel_f, 'o', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'none', ...
             'MarkerSize', 9, 'LineWidth', 1.5);
    end
    hold off;
    title(sprintf('%s | \\geq50%%: %d time-frequency points | Max: %d/25', ...
          row2_titles{col}, n_selected, max(counts_phase)), ...
          'FontSize', 9, 'FontWeight', 'bold');
    set(gca, 'FontSize', 9);
end

%% -- Overall title --
title(tl, {'\bfEEG ERSP Difference Maps and GA Consensus Features', ...
    'Row 1: Grand-average \DeltaERSP (dB) | Row 2: Time-frequency point selection frequency across 25 LOOCV folds (cyan \circ = \geq50% consensus)'}, ...
    'FontSize', 11);

%% --- Save ---
exportgraphics(fig, save_path, 'Resolution', 300);
fprintf('\nFigure saved to: %s\n', save_path);
