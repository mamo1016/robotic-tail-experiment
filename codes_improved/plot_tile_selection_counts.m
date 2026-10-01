%% plot_tile_selection_counts.m
% Supplementary figure: GA selection frequency across the 25 LOOCV folds for a
% single method (SVR), split by the three ERSP contrasts.
%
% Output: <repo>/paper/figures/tile_counts_svr.png
%
% PROVENANCE (resolved 2026-07-28)
%   Now pinned to the SAME run as Supplementary Fig. S3
%   (tile_selection_counts.png, all seven methods), so S2 is the SVR close-up
%   of S3 and their numbers agree:
%       ga_feature_30runs/toResult/results_LOOCV_2026-03-29_03-38.mat
%   SVR control-phase counts in this run: control 4, surprise 5, learning 3.
%
%   HISTORY: until 2026-07-28 this figure was built from
%   ga_feature/results_LOOCV_2026-03-06_14-47.mat, whose stored config carries
%   save_filename_override = loocv_results_CONTROL.mat, i.e. the CONTROL-contrast
%   validation run behind control_loocv.png. Its GA optimised against a null
%   target, so its selection pattern did not represent the main analysis, and its
%   counts (2/2/1) contradicted the caption's claim that the control contrast
%   shows no consistent convergence. Author decision: regenerate from S3's run.
%
%   Earlier still, the script loaded whichever file was newest in ga_feature/,
%   which is why the provenance had to be reconstructed after the fact. Do not
%   revert to "newest".

clear; close all; clc;

addpath('functions\for_obtainERD');
addpath('functions\for_GA');
config = configure_parameters();

SOURCE_FILE = 'results_LOOCV_2026-03-29_03-38.mat';
load_path   = fullfile(config.out_dir, 'ga_feature_30runs', 'toResult', SOURCE_FILE);
if ~isfile(load_path)
    error('plot_tile_selection_counts:sourceMissing', ...
          'Pinned source run not found:\n  %s', load_path);
end
fprintf('Loading: %s\n', load_path);
load(load_path, 'results', 'run_info');

method_to_plot = 'svr';
n_folds        = 25;

% tile_t / tile_f: use run_info if available, else use ersp_tile.m definition
if isfield(run_info, 'tile_t') && isfield(run_info, 'tile_f')
    tile_t = run_info.tile_t;
    tile_f = run_info.tile_f;
else
    tile_t = -2 : 0.2 : (3.0 - 0.2);   % 25 bins, -2s to 2.8s
    tile_f = 0  : 3   : (50  - 3);      % 16 bins, 0Hz to 47Hz
end
n_freq               = length(tile_f);
n_time               = length(tile_t);
n_tiles_per_contrast = n_freq * n_time;

% Phase names follow the manuscript standard: Control / Surprise / Learning.
c1 = 'Control (Expected vs. Expected)';
c2 = 'Surprise (Unexpected vs. Expected)';
c3 = 'Learning (Right After vs. Right Before)';
contrast_names = {c1, c2, c3};

if ~isfield(results, method_to_plot)
    avail = strjoin(fieldnames(results), ', ');
    error('Method "%s" not found. Available: %s', method_to_plot, avail);
end

counts = results.(method_to_plot).mask_x_counts;

fig = figure('Color', 'w', 'Position', [50 50 1400 450]);
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for c = 1:3
    i_start  = (c-1)*n_tiles_per_contrast + 1;
    i_end    = c*n_tiles_per_contrast;
    tile_map = reshape(counts(i_start:i_end), n_freq, n_time);

    nexttile;
    imagesc(tile_t, tile_f, tile_map, [0 n_folds]);
    axis xy;
    colorbar;
    colormap(hot);
    hold on;

    threshold = ceil(n_folds / 2);
    contour(tile_t, tile_f, tile_map, [threshold threshold], 'w-', 'LineWidth', 2);

    [max_val, max_idx] = max(tile_map(:));
    [max_f_bin, max_t_bin] = ind2sub(size(tile_map), max_idx);
    plot(tile_t(max_t_bin), tile_f(max_f_bin), 'co', 'MarkerSize', 12, 'LineWidth', 2);
    n_above = sum(tile_map(:) >= threshold);
    if n_above == 1
        unit_str = 'time-frequency point';
    else
        unit_str = 'time-frequency points';
    end
    % Two lines: the full contrast name plus the counts will not fit on one
    % line at this panel width and runs into the neighbouring panel.
    title({contrast_names{c}, ...
           sprintf('\\geq50%%: %d %s | Max: %d/%d', n_above, unit_str, max_val, n_folds)}, ...
          'FontSize', 9);
    xlabel('Time (s)', 'FontSize', 10);
    ylabel('Frequency (Hz)', 'FontSize', 10);
    ylim([5 35]);
    xline(0, 'k--', 'LineWidth', 0.8);
end

sgtitle(sprintf('Time-frequency point selection frequency [%s], n=%d folds', ...
        upper(method_to_plot), n_folds), 'FontSize', 13, 'FontWeight', 'bold');

fprintf('\n--- Selection summary [%s] ---\n', upper(method_to_plot));
fprintf('Threshold: >= %d of %d folds\n\n', ceil(n_folds/2), n_folds);
for c = 1:3
    i_start = (c-1)*n_tiles_per_contrast + 1;
    i_end   = c*n_tiles_per_contrast;
    blk     = counts(i_start:i_end);
    fprintf('%-42s %d above threshold, max=%d\n', ...
            contrast_names{c}, sum(blk >= ceil(n_folds/2)), max(blk));
end

paper_fig_dir = fullfile(fileparts(pwd), 'paper', 'figures');
out_path = fullfile(paper_fig_dir, 'tile_counts_svr.png');
exportgraphics(fig, out_path, 'Resolution', 300);
fprintf('Saved: %s\n', out_path);
