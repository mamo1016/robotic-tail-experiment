function verify_svr_consensus_tiles()
% verify_svr_consensus_tiles
% ----------------------------------------------------------------------
% Purpose
%   Verifies the SVR consensus-tile claims that were marked with
%   <span class="verify"> in paper/notes/plain_language_summaries.html.
%   Reads the canonical 100-run LOOCV archive, reproduces the consensus
%   detection logic of analyze_loocv_consistency.m, and prints a report
%   of every consensus tile's frequency, time, ΔERSP, and selection count.
%
% Inputs (hard-coded paths — edit if your output tree differs)
%   loocv_dir  : output/ga_feature_30runs/toResult
%   ersp_path  : output/all_in_one_pic/avergae_bootstrap_all_sub.mat
%   tile_path  : output/ga_feature/sub_tile.mat
%
% Outputs (saved to output/svr_verification/)
%   svr_verification_report.txt    — text report (one block per phase)
%   svr_consensus_tiles_visual.png — 1×3 figure with masked ΔERSP overlay
%
% Notes
%   - The "toResult" folder currently holds 102 LOOCV files but the source
%     code (run_loocv_30times.m line 154) labels it the "100-run archive".
%     This script enforces that intent: it sorts the files chronologically
%     and keeps the first 100, dropping the 2 latest extras.
%   - Consensus threshold matches analyze_loocv_consistency.m: a tile is a
%     consensus tile if it was selected in ≥50 % of runs (≥50/100 here).
%   - ΔERSP per tile is the mean of the pixel-resolution ERSP difference
%     within each tile's [time × frequency] cell, exactly as the right-hand
%     panel of loocv_consistency_combined_*_svr.png is produced.
% ----------------------------------------------------------------------

clear; close all; clc;

%% --- Paths --------------------------------------------------------------
loocv_dir   = 'output/ga_feature_30runs/toResult';
ersp_path   = 'output/all_in_one_pic/avergae_bootstrap_all_sub.mat';
tile_path   = 'output/ga_feature/sub_tile.mat';
output_dir  = 'output/svr_verification';

N_RUNS_CANONICAL = 100;

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

%% --- Add MATLAB paths so config + range_compute resolve -----------------
addpath('functions/for_obtainERD');
addpath('functions/for_ersp');

config = configure_parameters();

%% --- Locate and select the canonical 100 LOOCV files -------------------
loocv_files = dir(fullfile(loocv_dir, 'results_LOOCV_*.mat'));
fprintf('Found %d LOOCV result files in:\n  %s\n', length(loocv_files), loocv_dir);

if length(loocv_files) < N_RUNS_CANONICAL
    error('Need at least %d files, found %d.', N_RUNS_CANONICAL, length(loocv_files));
end

% Sort by filename (timestamp ordering), keep first N_RUNS_CANONICAL
[~, sort_idx] = sort({loocv_files.name});
loocv_files   = loocv_files(sort_idx);
loocv_files   = loocv_files(1:N_RUNS_CANONICAL);

fprintf('Using first %d files chronologically:\n', N_RUNS_CANONICAL);
fprintf('  First: %s\n', loocv_files(1).name);
fprintf('  Last : %s\n\n', loocv_files(end).name);

%% --- Load tile grid -----------------------------------------------------
load(tile_path, 'sub_tile');
first_sub  = fieldnames(sub_tile); first_sub  = first_sub{1};
first_type = fieldnames(sub_tile.(first_sub)); first_type = first_type{1};
tile_t = sub_tile.(first_sub).(first_type).t;
tile_f = sub_tile.(first_sub).(first_type).f;
n_tile_f = length(tile_f);
n_tile_t = length(tile_t);
fprintf('Tile grid: %d freq × %d time = %d tiles per phase\n', ...
    n_tile_f, n_tile_t, n_tile_f * n_tile_t);

%% --- Stack SVR masks across the 100 runs -------------------------------
n_tiles_total  = 1200;
svr_mask_stack = NaN(N_RUNS_CANONICAL, n_tiles_total);

for f_idx = 1:N_RUNS_CANONICAL
    fpath = fullfile(loocv_dir, loocv_files(f_idx).name);
    try
        d = load(fpath, 'results');
        if isfield(d.results, 'svr') && isfield(d.results.svr, 'mask_x')
            svr_mask_stack(f_idx, :) = d.results.svr.mask_x;
        else
            fprintf('  Warning: %s has no results.svr.mask_x — skipping\n', loocv_files(f_idx).name);
        end
    catch ME
        fprintf('  Warning: failed to load %s (%s)\n', loocv_files(f_idx).name, ME.message);
    end
end

valid_rows     = ~all(isnan(svr_mask_stack), 2);
n_valid        = sum(valid_rows);
svr_mask_stack = svr_mask_stack(valid_rows, :);
svr_mask_stack(isnan(svr_mask_stack)) = 0;

if n_valid == 0
    error('No SVR masks loaded — cannot proceed.');
end

fprintf('SVR masks loaded from %d / %d files\n\n', n_valid, N_RUNS_CANONICAL);

selection_count = sum(svr_mask_stack, 1);   % [1 × 1200]
threshold       = ceil(n_valid / 2);        % consensus = ≥50 %

%% --- Load ERSP phase differences ---------------------------------------
load(ersp_path, 'ave_boot_average');
[time_range, freq_range] = range_compute(config);

phase_names      = {'Control', 'Surprise', 'Learning'};
phase_idx_ranges = {1:400, 401:800, 801:1200};
phase_diffs      = { ...
    ave_boot_average.expected_duplicated.ersp - ave_boot_average.expected.ersp, ...
    ave_boot_average.unexpected.ersp          - ave_boot_average.expected.ersp, ...
    ave_boot_average.right_after.ersp         - ave_boot_average.right_before.ersp};

tile_dt = tile_t(2) - tile_t(1);
tile_df = tile_f(2) - tile_f(1);

%% --- Open report file --------------------------------------------------
report_path = fullfile(output_dir, 'svr_verification_report.txt');
fid         = fopen(report_path, 'w');
emit        = @(s) (fprintf('%s', s) + fprintf(fid, '%s', s));

emit(sprintf('=== SVR consensus tile verification ===\n'));
emit(sprintf('LOOCV directory : %s\n', loocv_dir));
emit(sprintf('Files used      : %d (first %d chronologically)\n', n_valid, N_RUNS_CANONICAL));
emit(sprintf('First file      : %s\n', loocv_files(1).name));
emit(sprintf('Last file       : %s\n', loocv_files(end).name));
emit(sprintf('Consensus rule  : selected in ≥%d / %d files (≥50%%)\n\n', threshold, n_valid));

%% --- Set up figure -----------------------------------------------------
fig = figure('Color', 'w', 'Position', [50 50 1500 500]);
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

cmap_bwr = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ...
            ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];

summary = struct();

%% --- Per-phase analysis ------------------------------------------------
for p = 1:length(phase_names)
    p_name     = phase_names{p};
    p_range    = phase_idx_ranges{p};
    p_count    = selection_count(p_range);
    consensus  = (p_count >= threshold);
    n_consens  = sum(consensus);

    emit(sprintf('--- Phase: %s ---\n', p_name));
    emit(sprintf('Consensus tiles (≥%d/%d): %d\n', threshold, n_valid, n_consens));

    % ΔERSP per tile (downsample full pixel grid to tile cells)
    dif_phase       = phase_diffs{p};
    delta_ersp_tile = NaN(n_tile_f, n_tile_t);
    for fi = 1:n_tile_f
        for ti = 1:n_tile_t
            row_px = freq_range >= (tile_f(fi) - tile_df/2) & freq_range <= (tile_f(fi) + tile_df/2);
            col_px = time_range >= (tile_t(ti) - tile_dt/2) & time_range <= (tile_t(ti) + tile_dt/2);
            patch  = dif_phase(row_px, col_px);
            if ~isempty(patch)
                delta_ersp_tile(fi, ti) = mean(patch(:), 'omitnan');
            end
        end
    end

    % Build per-tile records for the consensus set
    consens_idx = find(consensus);
    n_consens   = length(consens_idx);
    tile_idx_v  = zeros(n_consens, 1);
    fi_v        = zeros(n_consens, 1);
    ti_v        = zeros(n_consens, 1);
    freq_v      = zeros(n_consens, 1);
    time_v      = zeros(n_consens, 1);
    delta_v     = zeros(n_consens, 1);
    count_v     = zeros(n_consens, 1);

    for k = 1:n_consens
        tidx        = consens_idx(k);
        fi          = mod(tidx - 1, n_tile_f) + 1;
        ti          = floor((tidx - 1) / n_tile_f) + 1;
        tile_idx_v(k) = tidx;
        fi_v(k)       = fi;
        ti_v(k)       = ti;
        freq_v(k)     = tile_f(fi);
        time_v(k)     = tile_t(ti);
        delta_v(k)    = delta_ersp_tile(fi, ti);
        count_v(k)    = p_count(tidx);
    end

    % Sort by time then frequency
    if n_consens > 0
        [~, ord]   = sortrows([time_v, freq_v]);
        tile_idx_v = tile_idx_v(ord);
        freq_v     = freq_v(ord);
        time_v     = time_v(ord);
        delta_v    = delta_v(ord);
        count_v    = count_v(ord);
    end

    % Print table
    emit(sprintf('  %-6s  %-7s  %-8s  %-10s  %s\n', 'Idx', 'Freq', 'Time', 'dERSP', 'Count'));
    emit(sprintf('  %s\n', repmat('-', 1, 50)));
    for k = 1:n_consens
        emit(sprintf('  %-6d  %5.1fHz  %+6.2fs  %+8.4fdB  %d/%d\n', ...
            tile_idx_v(k), freq_v(k), time_v(k), delta_v(k), count_v(k), n_valid));
    end
    emit(sprintf('\n'));

    % Per-phase summary line
    if n_consens > 0
        emit(sprintf('  Summary:\n'));
        emit(sprintf('    N tiles    : %d\n', n_consens));
        emit(sprintf('    Freq range : %.1f – %.1f Hz\n', min(freq_v), max(freq_v)));
        time_str = strjoin(arrayfun(@(x) sprintf('%+.2f', x), time_v(:)', 'UniformOutput', false), ', ');
        emit(sprintf('    Time pts   : {%s} s\n', time_str));
        emit(sprintf('    dERSP rng  : %+.4f to %+.4f dB\n', min(delta_v), max(delta_v)));
        emit(sprintf('    Max count  : %d/%d\n\n', max(count_v), n_valid));

        summary.(p_name).n_tiles    = n_consens;
        summary.(p_name).freq_range = [min(freq_v), max(freq_v)];
        summary.(p_name).times      = time_v(:)';
        summary.(p_name).delta_rng  = [min(delta_v), max(delta_v)];
        summary.(p_name).max_count  = max(count_v);
    else
        emit(sprintf('  No consensus tiles.\n\n'));
        summary.(p_name).n_tiles = 0;
    end

    %% --- Visual subplot ------------------------------------------------
    nexttile;

    % Pixel mask for >=50% consensus tiles
    pixel_mask   = zeros(size(dif_phase));
    consens_2d   = reshape(consensus, [n_tile_f, n_tile_t]);
    for fi = 1:n_tile_f
        for ti = 1:n_tile_t
            if consens_2d(fi, ti) == 1
                row_keep = freq_range >= (tile_f(fi) - tile_df/2) & freq_range <= (tile_f(fi) + tile_df/2);
                col_keep = time_range >= (tile_t(ti) - tile_dt/2) & time_range <= (tile_t(ti) + tile_dt/2);
                pixel_mask(row_keep, col_keep) = 1;
            end
        end
    end
    masked_ersp = dif_phase .* pixel_mask;

    imagesc(time_range, freq_range, masked_ersp, [-0.1 0.1]);
    axis xy;
    colormap(gca, cmap_bwr);
    cbar = colorbar;
    cbar.Label.String = 'ΔERSP (dB)';
    hold on;
    xline(0, 'k--', 'LineWidth', 0.8);

    % Label each consensus tile
    for k = 1:n_consens
        text(time_v(k), freq_v(k), ...
            sprintf('%.0fHz\n%+.2fs\n%+.3fdB\n%d/%d', ...
                freq_v(k), time_v(k), delta_v(k), count_v(k), n_valid), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 7, 'FontWeight', 'bold', 'Color', 'k', ...
            'BackgroundColor', [1 1 1 0.75], 'Margin', 1);
    end
    hold off;

    ylim([5 35]);
    xlabel('Time (s)', 'FontSize', 11);
    ylabel('Frequency (Hz)', 'FontSize', 11);
    title(sprintf('%s — SVR consensus (%d tiles, >=%d/%d)', ...
        p_name, n_consens, threshold, n_valid), 'FontWeight', 'bold');
    set(gca, 'FontSize', 10);
end

%% --- Compact summary block (HTML cross-check) --------------------------
emit(sprintf('=== COMPACT SUMMARY for plain_language_summaries.html ===\n'));
emit(sprintf('Use these to confirm/replace the red <span class="verify"> values.\n\n'));

html_label = struct('Surprise', 'Surprise (HTML §2b)', ...
                    'Learning', 'Adaptation (HTML §2c)', ...
                    'Control',  'Control/Baseline (HTML §2d)');

for p_cell = {'Surprise', 'Learning', 'Control'}
    pn = p_cell{1};
    if summary.(pn).n_tiles > 0
        s        = summary.(pn);
        time_str = strjoin(arrayfun(@(x) sprintf('%+.2f', x), s.times, 'UniformOutput', false), ', ');
        emit(sprintf('%-26s: %d tiles | %.1f-%.1f Hz | times {%s} s | dERSP %+.4f to %+.4f dB | max %d/%d\n', ...
            html_label.(pn), s.n_tiles, s.freq_range(1), s.freq_range(2), ...
            time_str, s.delta_rng(1), s.delta_rng(2), s.max_count, n_valid));
    else
        emit(sprintf('%-26s: 0 consensus tiles\n', html_label.(pn)));
    end
end

emit(sprintf('\n=== END OF REPORT ===\n'));
fclose(fid);

%% --- Save figure -------------------------------------------------------
fig_path = fullfile(output_dir, 'svr_consensus_tiles_visual.png');
exportgraphics(fig, fig_path, 'Resolution', 300);

fprintf('\nReport : %s\n', report_path);
fprintf('Figure : %s\n', fig_path);

end
