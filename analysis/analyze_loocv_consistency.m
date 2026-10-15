function analyze_loocv_consistency(loocv_dir, pic_out_dir, ersp_data_dir, tile_data_dir, batch_id)
% analyze_loocv_consistency(loocv_dir, pic_out_dir, ersp_data_dir, tile_data_dir)
% Purpose: Find which tiles are consistently selected across multiple LOOCV result files
% Groups results by: Phase (Control/Surprise/Learning) × Method (spearman/pca/robust/pls/lasso/ridge/svr)
% Handles missing methods (not all runs have all 7 methods)
% Output: Heatmaps showing tile selection consistency across runs
%
% Inputs:
%   loocv_dir      — folder containing results_LOOCV_*.mat files
%   pic_out_dir    — folder where PNG outputs are saved
%   ersp_data_dir  — folder containing avergae_bootstrap_all_sub.mat (optional, defaults to config)
%   tile_data_dir  — folder containing sub_tile.mat (optional, defaults to config.out_dir/ga_feature)

config = configure_parameters();

if nargin < 3 || isempty(ersp_data_dir)
    ersp_data_dir = fullfile(config.out_dir, 'all_in_one_pic');
end

if nargin < 4 || isempty(tile_data_dir)
    tile_data_dir = fullfile(config.out_dir, 'ga_feature');
end

if nargin < 5 || isempty(batch_id)
    batch_id = '';
end

ga_dir = loocv_dir;

%% Find LOOCV result files (filtered by batch_id if provided)
if isempty(batch_id)
    loocv_files = dir(fullfile(ga_dir, 'results_LOOCV_*.mat'));
    fprintf('Found %d LOOCV result files (all batches).\n', length(loocv_files));
else
    loocv_files = dir(fullfile(ga_dir, sprintf('results_LOOCV_*_%s.mat', batch_id)));
    fprintf('Found %d LOOCV result files for batch_id=%s.\n', length(loocv_files), batch_id);
end

if isempty(loocv_files)
    error('No LOOCV results found in %s', ga_dir);
end

%% Load tile coordinates
load(fullfile(tile_data_dir, 'sub_tile.mat'), 'sub_tile');
first_sub  = fieldnames(sub_tile); first_sub  = first_sub{1};
first_type = fieldnames(sub_tile.(first_sub)); first_type = first_type{1};
tile_t = sub_tile.(first_sub).(first_type).t;
tile_f = sub_tile.(first_sub).(first_type).f;
n_tile_f = length(tile_f);
n_tile_t = length(tile_t);

fprintf('Tile grid: %d freq × %d time = %d tiles per phase\n', n_tile_f, n_tile_t, n_tile_f*n_tile_t);

%% Define phases and methods
phase_names      = {'Control', 'Surprise', 'Learning'};
phase_idx_ranges = {1:400, 401:800, 801:1200};
method_names     = {'spearman', 'pca', 'robust', 'pls', 'lasso', 'ridge', 'svr'};
n_phases = length(phase_names);
n_methods = length(method_names);
n_tiles_per_phase = 400;

%% Initialize consistency matrix
% consistency{phase_idx}{method_idx} = [n_files × 400] binary matrix
% Each row = one LOOCV file, each column = one tile
consistency = cell(n_phases, n_methods);
for p = 1:n_phases
    for m = 1:n_methods
        consistency{p, m} = [];
    end
end

file_dates = {}; % Track which file each row corresponds to

%% Load all LOOCV files and extract masks
for f_idx = 1:length(loocv_files)
    loocv_file = fullfile(ga_dir, loocv_files(f_idx).name);
    [~, file_name, ~] = fileparts(loocv_files(f_idx).name);

    try
        loocv_data = load(loocv_file, 'results');
    catch ME
        fprintf('ERROR loading %s: %s\n', file_name, ME.message);
        continue;
    end

    results = loocv_data.results;
    available_methods = fieldnames(results);

    fprintf('File %d/%d: %s\n  Available methods: %s\n', f_idx, length(loocv_files), file_name, strjoin(available_methods, ', '));

    file_dates{f_idx} = file_name;

    % For each phase
    for p = 1:n_phases
        p_range = phase_idx_ranges{p};

        % For each method (including missing ones)
        for m = 1:n_methods
            method = method_names{m};

            if isfield(results, method)
                % Extract mask for this phase
                mask_full = results.(method).mask_x;
                mask_phase = mask_full(p_range);  % 1×400 binary

                % Append to consistency matrix
                consistency{p, m} = [consistency{p, m}; mask_phase];
            else
                % Method not available in this file: append NaN row
                consistency{p, m} = [consistency{p, m}; NaN(1, n_tiles_per_phase)];
            end
        end
    end
end

fprintf('\n=== LOOCV Files Processed ===\n');
fprintf('Total files: %d\n', length(file_dates));

%% Compute consistency statistics for each phase × method
fprintf('\n=== Tile Selection Consistency ===\n\n');

if ~exist(pic_out_dir, 'dir')
    mkdir(pic_out_dir);
end

% Pre-load ERSP data for combined plotting
load(fullfile(ersp_data_dir, 'avergae_bootstrap_all_sub.mat'), 'ave_boot_average');
[time_range, freq_range] = range_compute(config);

% Define phase ERSP differences
phase_diffs = { ...
    ave_boot_average.expected_duplicated.ersp - ave_boot_average.expected.ersp, ...
    ave_boot_average.unexpected.ersp          - ave_boot_average.expected.ersp, ...
    ave_boot_average.right_after.ersp         - ave_boot_average.right_before.ersp };

% Colormap for ERSP
cmap_bwr = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ...
            ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];

tile_dt = tile_t(2) - tile_t(1);
tile_df = tile_f(2) - tile_f(1);

for p = 1:n_phases
    p_name = phase_names{p};
    fprintf('Phase: %s\n', p_name);
    dif_phase = phase_diffs{p};

    for m = 1:n_methods
        method = method_names{m};
        mask_matrix = consistency{p, m};  % [n_files × 400]

        % Count how many files have data for this method
        n_valid_files = sum(~all(isnan(mask_matrix), 2));

        if n_valid_files == 0
            fprintf('  %s: No data\n', method);
            continue;
        end

        % Replace NaN with 0 (missing methods = not selected)
        mask_matrix_clean = mask_matrix;
        mask_matrix_clean(isnan(mask_matrix)) = 0;

        % Count selection frequency: how many files selected each tile
        selection_count = sum(mask_matrix_clean, 1);  % [1 × 400]

        % Tiles selected in ALL valid files
        all_selected = (selection_count == n_valid_files);
        n_all = sum(all_selected);

        % Tiles selected in MOST files (>= 50%)
        most_selected = (selection_count >= ceil(n_valid_files / 2));
        n_most = sum(most_selected);

        % Tiles selected in SOME files (>= 1)
        some_selected = (selection_count >= 1);
        n_some = sum(some_selected);

        fprintf('  %s: %d files | %d tiles (all) | %d tiles (≥50%%) | %d tiles (≥1)\n', ...
            method, n_valid_files, n_all, n_most, n_some);

        % --- Create combined figure with consistency (LEFT) and masked ERSP (RIGHT) ---
        selection_2d = reshape(selection_count, [n_tile_f, n_tile_t]);

        % Compute >50% consistency mask for masked ERSP
        consistency_threshold = ceil(n_valid_files / 2);
        consistent_mask = (selection_count >= consistency_threshold);  % [1 × 400]
        n_consistent = sum(consistent_mask);

        % Create figure with 1 row, 2 columns
        fig = figure('Color', 'w', 'Position', [100 100 1200 420]);
        tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

        % LEFT: Consistency heatmap
        nexttile;
        imagesc(tile_t, tile_f, selection_2d, [0 n_valid_files]);
        axis xy;
        colormap(gca, 'hot');
        cbar_left = colorbar;
        cbar_left.Label.String = sprintf('Selection count\n(out of %d runs)', n_valid_files);
        hold on;
        xline(0, 'c--', 'LineWidth', 0.8);
        hold off;
        xlabel('Time (s)', 'FontSize', 12);
        ylabel('Frequency (Hz)', 'FontSize', 12);
        title(sprintf('%s — %s\nConsistency across runs', p_name, upper(method)), 'FontSize', 11, 'FontWeight', 'bold');
        set(gca, 'FontSize', 10);
        ylim([5 35]);

        % RIGHT: Masked ERSP (if there are consistent tiles)
        nexttile;
        if n_consistent > 0
            % Build full-resolution pixel mask for >50% consistent tiles
            pixel_mask = zeros(size(dif_phase));
            consistent_2d = reshape(consistent_mask, [n_tile_f, n_tile_t]);

            for fi = 1:n_tile_f
                for ti = 1:n_tile_t
                    if consistent_2d(fi, ti) == 1
                        row_keep = freq_range >= (tile_f(fi) - tile_df/2) & ...
                                   freq_range <= (tile_f(fi) + tile_df/2);
                        col_keep = time_range >= (tile_t(ti) - tile_dt/2) & ...
                                   time_range <= (tile_t(ti) + tile_dt/2);
                        pixel_mask(row_keep, col_keep) = 1;
                    end
                end
            end

            % Apply mask to ERSP
            masked_ersp = dif_phase .* pixel_mask;

            % Downsample to tile grid
            masked_tile = zeros(n_tile_f, n_tile_t);
            for fi = 1:n_tile_f
                for ti = 1:n_tile_t
                    row_px = freq_range >= (tile_f(fi) - tile_df/2) & freq_range <= (tile_f(fi) + tile_df/2);
                    col_px = time_range >= (tile_t(ti) - tile_dt/2) & time_range <= (tile_t(ti) + tile_dt/2);
                    patch = masked_ersp(row_px, col_px);
                    if isempty(patch)
                        masked_tile(fi, ti) = 0;
                    else
                        masked_tile(fi, ti) = mean(patch(:));
                    end
                end
            end
            masked_tile(isnan(masked_tile)) = 0;

            imagesc(tile_t, tile_f, masked_tile, [-0.1 0.1]);
            colormap(gca, cmap_bwr);
            axis xy;
            cbar_right = colorbar;
            cbar_right.Label.String = 'ΔERSP (dB)';
            hold on;
            xline(0, 'k--', 'LineWidth', 0.8);
            hold off;
            ylim([5 35]);
            xlabel('Time (s)', 'FontSize', 12);
            ylabel('Frequency (Hz)', 'FontSize', 12);
            title(sprintf('ΔERSP masked (≥%d runs)', consistency_threshold), 'FontSize', 11, 'FontWeight', 'bold');
            set(gca, 'FontSize', 10);
        else
            % No consistent tiles
            imagesc([-0.1 0.1], [5 35], zeros(2, 2));
            colormap(gca, cmap_bwr);
            axis xy;
            cbar_right = colorbar;
            cbar_right.Label.String = 'ΔERSP (dB)';
            xlabel('Time (s)', 'FontSize', 12);
            ylabel('Frequency (Hz)', 'FontSize', 12);
            title('No consistent time-frequency points', 'FontSize', 11, 'FontWeight', 'bold');
            set(gca, 'FontSize', 10);
        end

        if isempty(batch_id)
            out_path = fullfile(pic_out_dir, sprintf('loocv_consistency_combined_%s_%s.png', lower(p_name), lower(method)));
        else
            out_path = fullfile(pic_out_dir, sprintf('loocv_consistency_combined_%s_%s_%s.png', lower(p_name), lower(method), batch_id));
        end
        exportgraphics(fig, out_path, 'Resolution', 300);
        fprintf('    → Saved: %s\n', out_path);
        close(fig);
    end
    fprintf('\n');
end

%% Summary: Highly consistent tiles (selected in all or most runs)
fprintf('\n=== HIGHLY CONSISTENT TILES (ACROSS ALL RUNS) ===\n\n');

for p = 1:n_phases
    p_name = phase_names{p};
    fprintf('Phase: %s\n', p_name);

    for m = 1:n_methods
        method = method_names{m};
        mask_matrix = consistency{p, m};

        n_valid_files = sum(~all(isnan(mask_matrix), 2));
        if n_valid_files < 2
            continue;
        end

        mask_matrix_clean = mask_matrix;
        mask_matrix_clean(isnan(mask_matrix)) = 0;
        selection_count = sum(mask_matrix_clean, 1);

        % Tiles selected in ALL valid files
        all_selected_idx = find(selection_count == n_valid_files);

        if ~isempty(all_selected_idx)
            fprintf('  %s (all %d files):\n', method, n_valid_files);
            for i = 1:length(all_selected_idx)
                tile_idx = all_selected_idx(i);
                fi = mod(tile_idx - 1, n_tile_f) + 1;
                ti = floor((tile_idx - 1) / n_tile_f) + 1;
                fprintf('    Tile %d: t=%+.1fs f=%2dHz\n', tile_idx, tile_t(ti), tile_f(fi));
            end
        end
    end
end

fprintf('\n=== ANALYSIS COMPLETE ===\n');
end
