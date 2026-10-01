%% plot_loocv_brain_body_single_run.m
% Diagnostic: cross-validation brain–body scatter + tile selection counts
% from ONE results_LOOCV_*.mat file.
%
% Purpose:
%   Answers "why is the consistency map blank even though r is high?"
%   Shows (A) per-method predicted vs actual scatter and (B) raw tile
%   selection counts across the 25 LOOCV folds — so you can see whether
%   tiles are being selected at all and whether any reach the 50% threshold.
%
% Usage:
%   1. Set MAT_FILE below to the .mat you want to inspect (or leave empty
%      to pick the most-recent file automatically).
%   2. Run the script (F5 or Run).
%   3. Two figures are saved to OUT_DIR.

clear; close all; clc;

%% ── Configuration ─────────────────────────────────────────────────────────
% Which run to plot:
%   'canonical'  -> results_LOOCV_2026-02-24_15-32.mat (ga_feature/). This is
%                   the paper's primary run and drives panel A. It has NO
%                   mask_x_counts field, so panel B comes out BLANK.
%   'tilecounts' -> results_LOOCV_2026-03-29_03-38.mat (ga_feature_30runs/
%                   toResult/). Carries per-fold mask_x_counts, so both panels
%                   render. This is the run behind the supplementary figure
%                   tile_selection_counts.png.
%   'newest'/''  -> auto-pick the most recent file (diagnostic use only)
%   or an explicit path to inspect one particular run.
MAT_FILE = 'tilecounts';

% Where to save the output PNGs (defaults to all_in_one_pic)
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
addpath('functions\for_GA');
config    = configure_parameters();
OUT_DIR   = fullfile(config.out_dir, 'all_in_one_pic', 'brain_body_diagnostic');
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

% Tile metadata (needed to reshape counts into 2-D maps)
TILE_DATA_DIR = fullfile(config.out_dir, 'ga_feature');

%% ── Resolve MAT_FILE ───────────────────────────────────────────────────────
% The canonical run is the paper's primary analysis. Figures built from any
% other run will not match the numbers quoted in the captions -- this is what
% caused tile_selection_counts.png to disagree with its supplementary caption.
CANONICAL_FILE = 'results_LOOCV_2026-02-24_15-32.mat';
TILECOUNT_FILE = 'results_LOOCV_2026-03-29_03-38.mat';

if strcmp(MAT_FILE, 'canonical')
    MAT_FILE = fullfile(config.out_dir, 'ga_feature', CANONICAL_FILE);
    if ~isfile(MAT_FILE)
        error('plot_loocv_brain_body:canonicalMissing', ...
              'Canonical results file not found:\n  %s', MAT_FILE);
    end
elseif strcmp(MAT_FILE, 'tilecounts')
    MAT_FILE = fullfile(config.out_dir, 'ga_feature_30runs', 'toResult', TILECOUNT_FILE);
    if ~isfile(MAT_FILE)
        error('plot_loocv_brain_body:tilecountMissing', ...
              'Selection-count results file not found:\n  %s', MAT_FILE);
    end
elseif isempty(MAT_FILE) || strcmp(MAT_FILE, 'newest')
    search_dirs = {
        fullfile(config.out_dir, 'ga_feature_30runs', 'toResult'),
        fullfile(config.out_dir, 'ga_feature')
    };
    all_files = [];
    for d = 1:length(search_dirs)
        hits = dir(fullfile(search_dirs{d}, 'results_LOOCV_*.mat'));
        if ~isempty(hits)
            all_files = [all_files; hits]; %#ok<AGROW>
        end
    end
    if isempty(all_files)
        error('No results_LOOCV_*.mat files found. Set MAT_FILE manually.');
    end
    [~, newest_idx] = max([all_files.datenum]);
    MAT_FILE = fullfile(all_files(newest_idx).folder, all_files(newest_idx).name);
end

fprintf('Loading: %s\n', MAT_FILE);
d = load(MAT_FILE, 'results');
results = d.results;
method_names = fieldnames(results);
fprintf('Methods found: %s\n', strjoin(method_names', ', '));

%% ── Load tile grid ─────────────────────────────────────────────────────────
load(fullfile(TILE_DATA_DIR, 'sub_tile.mat'), 'sub_tile');
first_sub  = fieldnames(sub_tile); first_sub  = first_sub{1};
first_type = fieldnames(sub_tile.(first_sub)); first_type = first_type{1};
tile_t = sub_tile.(first_sub).(first_type).t;
tile_f = sub_tile.(first_sub).(first_type).f;
n_tile_f = length(tile_f);
n_tile_t = length(tile_t);
phase_names  = {'Control', 'Surprise', 'Learning'};
phase_ranges = {1:400, 401:800, 801:1200};

n_methods = length(method_names);
colors = lines(n_methods);

%% ══════════════════════════════════════════════════════════════════════════
%% FIGURE A: Brain–Body Scatter (predicted vs actual, each subject = 1 dot)
%% ══════════════════════════════════════════════════════════════════════════
n_cols = ceil(n_methods / 2);
n_rows = 2;
fig_a = figure('Color', 'w', 'Position', [50 50 300*n_cols 560]);
tiledlayout(n_rows, n_cols, 'Padding', 'compact', 'TileSpacing', 'compact');

fprintf('\n%-12s  %6s  %7s\n', 'Method', 'r', 'p');
fprintf('%s\n', repmat('-', 1, 30));

for m = 1:n_methods
    mname = method_names{m};
    pred  = results.(mname).predicted;
    act   = results.(mname).actual;

    % Guard: skip if both are flat
    if std(pred) < 1e-6 && std(act) < 1e-6
        fprintf('%-12s  SKIPPED (all zeros)\n', mname);
        continue;
    end

    [r_val, p_val] = corr(pred, act, 'Type', 'Spearman');
    fprintf('%-12s  %6.3f  %7.4f\n', mname, r_val, p_val);

    nexttile;
    scatter(pred, act, 70, colors(m,:), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    hold on;
    % Least-squares line through the scatter
    p_fit = polyfit(pred, act, 1);
    x_range = linspace(min(pred), max(pred), 100);
    plot(x_range, polyval(p_fit, x_range), 'k-', 'LineWidth', 1.5);
    hold off;
    title(sprintf('%s\nr=%.3f, p=%.4f', upper(mname), r_val, p_val), ...
          'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Predicted (brain)', 'FontSize', 9);
    ylabel('Actual (body)', 'FontSize', 9);
    grid on; box on;
    axis tight;
end

sgtitle('LOOCV Brain–Body Predictions', ...
        'FontSize', 10, 'Interpreter', 'none');

fig_a_path = fullfile(OUT_DIR, 'brain_body_scatter.png');
exportgraphics(fig_a, fig_a_path, 'Resolution', 300);
fprintf('\nSaved: %s\n', fig_a_path);
close(fig_a);

%% ══════════════════════════════════════════════════════════════════════════
%% FIGURE B: Tile Selection Counts (raw fold counts, before 50% threshold)
%% Purpose:  If counts are all low → GA is selecting very few tiles per
%%           fold → sparsity penalty is dominating → blank consistency map
%% ══════════════════════════════════════════════════════════════════════════
n_subs = length(results.(method_names{1}).predicted);  % = 25
threshold = n_subs / 2;  % 50% of folds = 12.5

fig_b = figure('Color', 'w', 'Position', [50 50 420*3 300*n_methods]);
tiledlayout(n_methods, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

fprintf('\n%-12s  %-10s  %8s  %8s  %8s\n', 'Method', 'Phase', 'MaxCount', '>50%', 'TotalSel');
fprintf('%s\n', repmat('-', 1, 55));

for m = 1:n_methods
    mname = method_names{m};

    if ~isfield(results.(mname), 'mask_x_counts')
        fprintf('%-12s  no mask_x_counts field — skipping\n', mname);
        continue;
    end

    counts = results.(mname).mask_x_counts;  % [1 × 1200]

    for p = 1:3
        phase_range = phase_ranges{p};
        phase_counts = counts(phase_range);          % [1 × 400]
        count_map    = reshape(phase_counts, n_tile_f, n_tile_t);

        above_thresh = sum(phase_counts >= threshold);
        fprintf('%-12s  %-10s  %8.0f  %8d  %8d\n', ...
            mname, phase_names{p}, max(phase_counts), above_thresh, sum(phase_counts > 0));

        nexttile;
        imagesc(tile_t, tile_f, count_map);
        colormap(gca, 'hot');
        colorbar;
        clim([0 n_subs]);  % Full range: 0 = never selected, 25 = always selected

        % Overlay 50% contour in cyan
        hold on;
        contour(tile_t, tile_f, count_map, [threshold threshold], 'c-', 'LineWidth', 1.5);
        hold off;

        title(sprintf('%s — %s\nmax=%d, above 50%%=%d time-frequency points', ...
              upper(mname), phase_names{p}, max(phase_counts), above_thresh), ...
              'FontSize', 8);
        xlabel('Time (s)', 'FontSize', 7);
        ylabel('Frequency (Hz)', 'FontSize', 7);
    end
end

sgtitle(sprintf('Time-frequency point selection counts per fold (0=never, %d=always)\nCyan contour = 50%% threshold (need ≥%.0f/%d folds)', ...
        n_subs, threshold, n_subs), ...
        'FontSize', 9, 'Interpreter', 'none');

fig_b_path = fullfile(OUT_DIR, 'tile_selection_counts.png');
exportgraphics(fig_b, fig_b_path, 'Resolution', 300);
fprintf('\nSaved: %s\n', fig_b_path);
close(fig_b);

fprintf('\n── Summary ────────────────────────────────────────────────────\n');
fprintf('If "above 50%%" is 0 for all phases → consistency map will be BLANK\n');
fprintf('because no tile was selected in ≥%.0f of the %d LOOCV folds.\n', threshold, n_subs);
fprintf('This is expected behaviour when sparsity penalty is ON (USE_SPARSITY_PENALTY=1)\n');
fprintf('and the GA makes different sparse selections each fold.\n');
fprintf('Output folder: %s\n', OUT_DIR);
