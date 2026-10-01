%% Improved Pipeline: EEG and Foot CoP Analysis aa
% This script replaces obtainSurpriseERDFaster_main.m
% 
% Improvements:
% 1. Removed redundant clear/clc between sections.
% 2. Fixed Data Leakage in GA validation (Relies exclusively on LOOCV for final reporting).
% 3. Fixed Statistical Flaw in Control Split (Ensures disjoint subsets for expected vs expected_duplicated).
% 4. Consolidated path variables and removed double dipping.
% 5. (Optimized performance by removing massive redundant loads in Step 3 and 4)

%% ------------------------------------------------------------------------- %
% 1. EEG Offline Signal Processing Pipeline(Epoching) %
% ------------------------------------------------------------------------- %
clear; close all; clc;
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
addpath('functions\for_GA');
config = configure_parameters();
filter_para = design_filter(config);
parallelPool_gpu_starter(config);
out_dir = config.out_dir;
dirs = ["all_in_one_pic", "epochs", "ersp", "ersp_ave", "ersp_each_sub_ave", "foot", "foot_pic", "foot2", "ga_feature", "interpolate", "rightafter_rightbefore_pic", "unexpected_expected_pic"];
for i = 1:length(dirs), chk_dir(fullfile(out_dir, dirs(i))); end

fprintf('\n=== STEP 1: Epoch Extraction ===\n');
subjects = dir(fullfile(config.baseDir, 'Sub*'));
fig_epoch = figure('Color', 'w');

for sub_num = 1:length(subjects) 
    subject_folder = fullfile(config.baseDir, subjects(sub_num).name);
    fprintf('Processing %s...\n', subjects(sub_num).name);
    files = all_files(subject_folder);

    for session = 1:length(files.eeg_files)
        [matchedfiles, timeDiff] = eeg_emg_timeoffset(files, session, config);
        adjusted_times.eeg_to_emg_offset = seconds(timeDiff);

        single_session_dataset = load_data_optimized(subject_folder, matchedfiles, config);
        single_session_dataset = footdist(single_session_dataset);
        single_session_dataset = footCompute_optimized(single_session_dataset);

        single_session_dataset.eeg_filtered = filter_eeg_vectorized(single_session_dataset.eeg_raw, filter_para, config);

        epoch_data = struct();
        epoch_data.times_eeg_st_expected = epoch_index(single_session_dataset, adjusted_times.eeg_to_emg_offset, config, 1);
        epoch_data.times_eeg_st_unexpected = epoch_index(single_session_dataset, adjusted_times.eeg_to_emg_offset, config, 2);
        epoch_data.times_eeg_st_right_before = epoch_index(single_session_dataset, adjusted_times.eeg_to_emg_offset, config, 3);
        epoch_data.times_eeg_st_right_after = epoch_index(single_session_dataset, adjusted_times.eeg_to_emg_offset, config, 4);

        [single_session_dataset] = foot_load_average(single_session_dataset, config);

        if length(epoch_data.times_eeg_st_expected) <= 10 || length(epoch_data.times_eeg_st_unexpected) <= 3 
            disp(sub_num)
            continue;
        end
        epoch_data = epoch_organise(config, single_session_dataset, adjusted_times, epoch_data, 0);

        plot_switch = 1;
        epoch_data.epochs_expected = epochs_compute(epoch_data, epoch_data.idxs_expected, plot_switch, fig_epoch, 'blue');
        epoch_data.epochs_unexpected = epochs_compute(epoch_data, epoch_data.idxs_unexpected, plot_switch, fig_epoch, 'red');
        epoch_data.epochs_right_before = epochs_compute(epoch_data, epoch_data.idxs_right_before, plot_switch, fig_epoch, 'green');
        epoch_data.epochs_right_after = epochs_compute(epoch_data, epoch_data.idxs_right_after, plot_switch, fig_epoch, 'magenta');

        filename = fullfile(out_dir, "epochs", "subject_" + sub_num + "___session" + session + "___epoch_data.mat");
        save(filename, 'epoch_data');
    end 
end

%% ------------------------------------------------------------------------- %
% 2. Epoch Interpolation & ERSP Computation %
% ------------------------------------------------------------------------- %
clear; close all; clc;
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
addpath('functions\for_GA');
config = configure_parameters();
out_dir = config.out_dir;

fprintf('\n=== STEP 2: Interpolation & ERSP Computation ===\n');
all_epoch = all_epoch_organise(config);
emg_ave_plot(all_epoch, config);

interpolate_all_epoch = struct();
interpolate_all_epoch = cleaning_epoch(all_epoch, interpolate_all_epoch);
save(fullfile(out_dir, "interpolate", "interpolate_all_epoch.mat"), 'interpolate_all_epoch', '-v7.3');

fields_interp = fieldnames(interpolate_all_epoch);
ersp_compute(fields_interp, interpolate_all_epoch, config, out_dir);

boot_strap_pick_num = min_epoch(interpolate_all_epoch);
save(fullfile(out_dir, "boot_strap_pick_num.mat"), 'boot_strap_pick_num', '-v7.3');


%% ------------------------------------------------------------------------- %
% 3. EEG Bootstrapping (Fixed Split Logic) %
% ------------------------------------------------------------------------- %
clear; close all; clc;
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
addpath('functions\for_GA');
config = configure_parameters();
out_dir = config.out_dir;
% Load boot_strap_pick_num from disk
load(fullfile(out_dir, "boot_strap_pick_num.mat"), 'boot_strap_pick_num');

fprintf('\n=== STEP 3: EEG Bootstrapping (Control Fix & Optimized Memory Processing) ===\n');
subjects_ersp = dir(fullfile(out_dir, "ersp_ave", '*_sub_*'));
subject_number_list = zeros(length(subjects_ersp), 1);

for list_num = 1:length(subject_number_list) 
    data = load(fullfile(subjects_ersp(list_num).folder, "ersp_epoch_ave_sub_" + list_num + ".mat")).ersp_epoch;
    current_name = fieldnames(data);
    tokens = regexp(current_name{1}, 'subject(\d+)', 'tokens');
    subject_number_list(list_num) = str2double(tokens{1}{1});
end 
save(fullfile(out_dir, "subject_number_list.mat"), 'subject_number_list', '-v7.3');

trial_type = {'expected', 'unexpected', 'right_before', 'right_after'};
unique_num = unique(subject_number_list);

fprintf('Pre-loading EEG ERSP data for all subjects into RAM...\n');
all_sub_ersp_data = struct();
for struct_num = 1:length(unique_num) 
    current_sub_num = find(subject_number_list == unique_num(struct_num));
    for trial_type_num = 1:length(trial_type) 
        trial_type_name = trial_type{trial_type_num};

        erspA = [];
        erspB = [];
        erspC = [];
        for i = 1:length(current_sub_num) 
            data = load(fullfile(subjects_ersp(1).folder, "ersp_epoch_ave_sub_" + current_sub_num(i) + ".mat")).ersp_epoch;
            current_name = fieldnames(data);
            ersp_data_boot = data.(current_name{1}).(trial_type_name).data_ersp.all_trial_ersp_db_average;
            if i == 1 
                erspA = ersp_data_boot;
            elseif i == 2
                erspB = ersp_data_boot;
            elseif i == 3 
                erspC = ersp_data_boot;
            end 
        end 
        all_sub_ersp_data(struct_num).(trial_type_name) = cat(1, erspA, erspB, erspC);
    end 
end

num_loops = 200;
boot_average = struct();

fprintf('Running %d fast in-memory bootstrap iterations...\n', num_loops);
for struct_num = 1:length(unique_num) 
    sum_ersp_epoch_bootstraped = struct();
    for loop = 1:num_loops 
        ersp_epoch_bootstraped = struct();

        for trial_type_num = 1:length(trial_type) 
            trial_type_name = trial_type{trial_type_num};
            ersp_data_boot = all_sub_ersp_data(struct_num).(trial_type_name);
            [n_trials, ~, ~] = size(ersp_data_boot);
            rng('shuffle');

            if strcmp(trial_type_name, 'expected')
                if n_trials >= 2 * boot_strap_pick_num 
                    idx_both = randperm(n_trials, 2 * boot_strap_pick_num);
                    select_idx_dup = idx_both(1:boot_strap_pick_num);
                    select_idx_exp = idx_both(boot_strap_pick_num + 1:end);
                else 
                    half = floor(n_trials / 2);
                    idx_both = randperm(n_trials);
                    pool_dup = idx_both(1:half);
                    pool_exp = idx_both(half + 1:end);
                    select_idx_dup = pool_dup(randi(length(pool_dup), 1, boot_strap_pick_num));
                    select_idx_exp = pool_exp(randi(length(pool_exp), 1, boot_strap_pick_num));
                end

                [final_ersp, final_std] = boot_strap(config, ersp_data_boot, boot_strap_pick_num, select_idx_dup);
                ersp_epoch_bootstraped.expected_duplicated.final_ersp = final_ersp;
                ersp_epoch_bootstraped.expected_duplicated.final_std = final_std;

                [final_ersp, final_std] = boot_strap(config, ersp_data_boot, boot_strap_pick_num, select_idx_exp);
                ersp_epoch_bootstraped.expected.final_ersp = final_ersp;
                ersp_epoch_bootstraped.expected.final_std = final_std;

            else 
                select_idx = randperm(n_trials, boot_strap_pick_num);
                [final_ersp, final_std] = boot_strap(config, ersp_data_boot, boot_strap_pick_num, select_idx);
                ersp_epoch_bootstraped.(trial_type_name).final_ersp = final_ersp;
                ersp_epoch_bootstraped.(trial_type_name).final_std = final_std;
            end 
        end

        types = fieldnames(ersp_epoch_bootstraped);
        for type_num = 1:length(types) 
            t_name = types{type_num};
            if loop == 1 
                sum_ersp_epoch_bootstraped.(t_name).final_ersp = ersp_epoch_bootstraped.(t_name).final_ersp;
                sum_ersp_epoch_bootstraped.(t_name).final_std = ersp_epoch_bootstraped.(t_name).final_std;
            else 
                sum_ersp_epoch_bootstraped.(t_name).final_ersp = sum_ersp_epoch_bootstraped.(t_name).final_ersp + ersp_epoch_bootstraped.(t_name).final_ersp;
                sum_ersp_epoch_bootstraped.(t_name).final_std = sum_ersp_epoch_bootstraped.(t_name).final_std + ersp_epoch_bootstraped.(t_name).final_std;
            end 
        end 
    end

    ave_single_sub = struct();
    for type_num = 1:length(types) 
        t_name = types{type_num};
        ave_single_sub.(t_name).ersp = sum_ersp_epoch_bootstraped.(t_name).final_ersp / num_loops;
        ave_single_sub.(t_name).std = sum_ersp_epoch_bootstraped.(t_name).final_std / num_loops;
    end

    filename = fullfile(out_dir, "ersp", sprintf("ersp_epoch_bootstraped_subject%d_final.mat", struct_num));
    save(filename, 'ave_single_sub', '-v7.3');
    boot_average.("sub_" + struct_num) = ave_single_sub;
end

[time_range, freq_range] = range_compute(config);
fig_ersp = figure('Color', 'w');

for sub_num = 1:length(unique_num) 
    ave_ersp_single_sub_optimized(boot_average.("sub_" + sub_num), 1, time_range, freq_range, sub_num, fig_ersp);
end 
save(fullfile(out_dir, "all_in_one_pic", "avergae_bootstrap_each_sub.mat"), 'boot_average', '-v7.3');

% Compute Grand Average over all subjects 
ave_boot_average = struct();
num_valid_subs = length(unique_num);
for sub_num = 1:num_valid_subs 
    types = fieldnames(boot_average.("sub_" + sub_num));
    for type_num = 1:length(types) 
        if sub_num == 1 
            ave_boot_average.(types{type_num}).ersp = boot_average.("sub_" + sub_num).(types{type_num}).ersp;
        else 
            ave_boot_average.(types{type_num}).ersp = ave_boot_average.(types{type_num}).ersp + boot_average.("sub_" + sub_num).(types{type_num}).ersp;
        end 
        if sub_num == num_valid_subs 
            ave_boot_average.(types{type_num}).ersp = ave_boot_average.(types{type_num}).ersp / num_valid_subs;
        end 
    end 
end 
save(fullfile(out_dir, "all_in_one_pic", "avergae_bootstrap_all_sub.mat"), 'ave_boot_average', '-v7.3');

%% ------------------------------------------------------------------------- %
% 3b. Generate Figure 3 Clean (Top Row Only, a-c)                             %
% Run this section standalone after Step 3 has been run at least once.        %
% Loads saved grand average from disk — no recomputation needed.              %
% ------------------------------------------------------------------------- %
clear; close all; clc;
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
addpath('functions\for_GA');
config = configure_parameters();
out_dir = config.out_dir;

fprintf('\n=== STEP 3b: Generate Figure 3 Clean (Top Row Only) ===\n');

load(fullfile(out_dir, "all_in_one_pic", "avergae_bootstrap_all_sub.mat"), 'ave_boot_average');
[time_range, freq_range] = range_compute(config);
types = fieldnames(ave_boot_average);  % expected, expected_duplicated, unexpected, right_before, right_after

fig_figure3_clean = figure('Color', 'w', 'Position', [100 100 1400 400]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel (a): Expected vs Expected (Control baseline)
dif_ctrl = ave_boot_average.expected_duplicated.ersp - ave_boot_average.expected.ersp;
visualisation_ersp(dif_ctrl, fig_figure3_clean, time_range, freq_range, [-0.5 0.5]);
title("(a) Expected vs Expected" + newline + "[Control baseline]", 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Frequency (Hz)'); xlabel('Time (s)');
ylim([5 40]); xlim([-2 3]);

% Panel (b): Unexpected vs Expected (Surprise signal)
dif_surprise = ave_boot_average.unexpected.ersp - ave_boot_average.expected.ersp;
visualisation_ersp(dif_surprise, fig_figure3_clean, time_range, freq_range, [-0.5 0.5]);
title("(b) Unexpected vs Expected" + newline + "[Surprise signal]", 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Frequency (Hz)'); xlabel('Time (s)');
ylim([5 40]); xlim([-2 3]);

% Panel (c): Right After vs Right Before (Learning signal)
dif_learning = ave_boot_average.right_after.ersp - ave_boot_average.right_before.ersp;
visualisation_ersp(dif_learning, fig_figure3_clean, time_range, freq_range, [-0.5 0.5]);
title("(c) Right After vs Right Before" + newline + "[Learning signal]", 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Frequency (Hz)'); xlabel('Time (s)');
ylim([5 40]); xlim([-2 3]);

out_path_fig3 = fullfile(out_dir, "all_in_one_pic", "figure3_top_row_only.png");
exportgraphics(fig_figure3_clean, out_path_fig3, 'Resolution', 300);
fprintf('Saved: Figure 3 Top Row (a-c) Clean Version\n  → %s\n', out_path_fig3);

%% ------------------------------------------------------------------------- %
% 3c. Consensus Tiles Only — GA-selected ERSP regions (Surprise contrast)     %
% Standalone section. Loads all data from disk.                               %
% Shows grand-average ERSP with only the >=50% GA consensus tiles visible.    %
% ------------------------------------------------------------------------- %
clear; close all; clc;
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
addpath('functions\for_GA');
config = configure_parameters();
out_dir = config.out_dir;

fprintf('\n=== STEP 3c: Consensus Tiles Only (All 3 Phases) ===\n');

% --- Load data ---
load(fullfile(out_dir, 'all_in_one_pic', 'avergae_bootstrap_all_sub.mat'), 'ave_boot_average');
load(fullfile(out_dir, 'ga_feature', 'sub_tile.mat'), 'sub_tile');
ga_dir      = fullfile(out_dir, 'ga_feature');
today_str   = datestr(now, 'yyyy-mm-dd');
today_files = dir(fullfile(ga_dir, ['results_LOOCV_' today_str '_*.mat']));
if isempty(today_files)
    error('No LOOCV files found for today (%s). Run Step 6 first.', today_str);
end
fprintf('Found %d LOOCV file(s) for today (%s).\n', length(today_files), today_str);
[time_range, freq_range] = range_compute(config);

% Tile grid coordinates
first_sub  = fieldnames(sub_tile); first_sub  = first_sub{1};
first_type = fieldnames(sub_tile.(first_sub)); first_type = first_type{1};
tile_t = sub_tile.(first_sub).(first_type).t;
tile_f = sub_tile.(first_sub).(first_type).f;
n_tile_f = length(tile_f);
n_tile_t = length(tile_t);
tile_dt = tile_t(2) - tile_t(1);
tile_df = tile_f(2) - tile_f(1);

% Grand-average ERSP differences for all 3 phases (full resolution)
phase_names      = {'Control', 'Surprise', 'Learning'};
phase_idx_ranges = {1:400, 401:800, 801:1200};
phase_diffs      = { ...
    ave_boot_average.expected_duplicated.ersp - ave_boot_average.expected.ersp, ...
    ave_boot_average.unexpected.ersp          - ave_boot_average.expected.ersp, ...
    ave_boot_average.right_after.ersp         - ave_boot_average.right_before.ersp };
str_cell = {"spearman","pca","robust","pls","lasso","ridge","svr"};
for f_idx = 1:length(today_files)
    loocv_file = fullfile(ga_dir, today_files(f_idx).name);
    [~, loocv_name, ~] = fileparts(today_files(f_idx).name);
    fprintf('\n--- File %d/%d: %s ---\n', f_idx, length(today_files), loocv_name);
    loocv_d = load(loocv_file, 'results');
    for str_i = 1:6
    for p_idx = 1:3
        dif_phase = phase_diffs{p_idx};
        p_name    = phase_names{p_idx};
        p_range   = phase_idx_ranges{p_idx};
        deleteFig();
        % GA consensus mask for current phase
        mask_phase   = loocv_d.results.(str_cell{str_i}).mask_x(p_range);        % 1x400 binary
        counts_phase = loocv_d.results.(str_cell{str_i}).mask_x_counts(p_range);
        mask_2d   = reshape(mask_phase,   [n_tile_f, n_tile_t]);
        counts_2d = reshape(counts_phase, [n_tile_f, n_tile_t]);
        
        % Build full-resolution binary mask: keep only pixels inside consensus tiles
        pixel_mask = zeros(size(dif_phase));
        for fi = 1:n_tile_f
            for ti = 1:n_tile_t
                if mask_2d(fi, ti) == 1
                    row_keep = freq_range >= (tile_f(fi) - tile_df/2) & ...
                               freq_range <= (tile_f(fi) + tile_df/2);
                    col_keep = time_range >= (tile_t(ti) - tile_dt/2) & ...
                               time_range <= (tile_t(ti) + tile_dt/2);
                    pixel_mask(row_keep, col_keep) = 1;
                end
            end
        end
        
        % Apply mask — zero out non-consensus pixels
        masked_ersp = dif_phase .* pixel_mask;
        
        % --- Downsample masked_ersp to 400-tile grid (same as LOOCV_Plot style) ---
        % Average pixels within each tile to get [n_tile_f x n_tile_t] matrix
        masked_tile = zeros(n_tile_f, n_tile_t);
        for fi = 1:n_tile_f
            for ti = 1:n_tile_t
                row_px = freq_range >= (tile_f(fi) - tile_df/2) & ...
                         freq_range <= (tile_f(fi) + tile_df/2);
                col_px = time_range >= (tile_t(ti) - tile_dt/2) & ...
                         time_range <= (tile_t(ti) + tile_dt/2);
                patch = masked_ersp(row_px, col_px);
                if isempty(patch)
                    masked_tile(fi, ti) = 0;
                else
                    masked_tile(fi, ti) = mean(patch(:));
                end
            end
        end
        masked_tile(isnan(masked_tile)) = 0;  % replace any NaN with 0 (shows as white)
        
        % --- Downsample full ERSP to 400-tile grid (unmasked) ---
        % Average pixels within each tile to get [n_tile_f x n_tile_t] matrix
        full_tile = zeros(n_tile_f, n_tile_t);
        for fi = 1:n_tile_f
            for ti = 1:n_tile_t
                row_px = freq_range >= (tile_f(fi) - tile_df/2) & ...
                         freq_range <= (tile_f(fi) + tile_df/2);
                col_px = time_range >= (tile_t(ti) - tile_dt/2) & ...
                         time_range <= (tile_t(ti) + tile_dt/2);
                patch = dif_phase(row_px, col_px);
                if isempty(patch)
                    full_tile(fi, ti) = 0;
                else
                    full_tile(fi, ti) = mean(patch(:));
                end
            end
        end
        full_tile(isnan(full_tile)) = 0;  % replace any NaN with 0
        
        % --- Plot tiled version ---
        cmap_bwr = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ...
                    ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];
        
        fig_masked = figure('Color', 'w', 'Position', [100 100 600 420]);
        tiledlayout(1, 1, 'Padding', 'compact');
        nexttile;
        imagesc(tile_t, tile_f, masked_tile, [-0.5 0.5]);
        axis xy;
        colormap(cmap_bwr); colorbar;
        hold on;
        xline(0, 'k--', 'LineWidth', 0.8);
        % --- Diagnostic: print ΔERSP value at each consensus tile ---
        sel_tiles = find(mask_phase == 1);
        fprintf('\n=== Diagnostic: ΔERSP at consensus tiles (phase: %s, method: %s, file: %s) ===\n', p_name, str_cell{str_i}, loocv_name);
        for dbg_i = 1:length(sel_tiles)
            ti_idx  = sel_tiles(dbg_i);
            fi_tile = mod(ti_idx-1, n_tile_f) + 1;
            tt_tile = floor((ti_idx-1)/n_tile_f) + 1;
            val     = masked_tile(fi_tile, tt_tile);
            cnt     = counts_2d(fi_tile, tt_tile);
            fprintf('  tile %3d: t=%+.1fs  f=%2dHz  ΔERSP=%+.4f  counts=%d/25\n', ...
                ti_idx, tile_t(tt_tile), tile_f(fi_tile), val, cnt);
        end

        % Add count labels on each consensus tile
        for i = 1:length(sel_tiles)
            ti_idx  = sel_tiles(i);
            fi_tile = mod(ti_idx-1, n_tile_f) + 1;
            tt_tile = floor((ti_idx-1)/n_tile_f) + 1;
            cnt     = counts_2d(fi_tile, tt_tile);
            text(tile_t(tt_tile), tile_f(fi_tile+1), sprintf('%d/25', cnt), ...
                 'Color', 'k', 'FontSize', 9, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        hold off;
        ylim([5 35]);
        xlabel('Time (s)', 'FontSize', 12);
        ylabel('Frequency (Hz)', 'FontSize', 12);
        title({['\bfGA Consensus Tiles — ' p_name], ...
               sprintf('Grand-average \\DeltaERSP (dB, tiled) | %d tiles \\geq50%% of 25 LOOCV folds | non-selected = 0', sum(mask_phase)), ...
               sprintf('method: %s', str_cell{str_i})}, ...
              'FontSize', 10);
        set(gca, 'FontSize', 11);

        out_path_masked = fullfile(out_dir, 'all_in_one_pic', string("figure3_consensus_masked_"+ loocv_name +"_"+ p_name +"_"+ str_cell{str_i} +".png"));
        exportgraphics(fig_masked, out_path_masked, 'Resolution', 300);
        fprintf('Saved: Consensus Masked ERSP (tiled)\n  → %s\n', out_path_masked);
    end  % p_idx phase loop
    end  % str_i methods loop
end  % f_idx files loop
% --- Plot full (unmasked) ERSP map for all 3 phases ---
cmap_bwr_full = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ...
                 ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];
for p_idx = 1:3
    dif_phase   = phase_diffs{p_idx};
    p_name      = phase_names{p_idx};
    full_tile_p = zeros(n_tile_f, n_tile_t);
    for fi = 1:n_tile_f
        for ti = 1:n_tile_t
            row_px = freq_range >= (tile_f(fi) - tile_df/2) & freq_range <= (tile_f(fi) + tile_df/2);
            col_px = time_range >= (tile_t(ti) - tile_dt/2) & time_range <= (tile_t(ti) + tile_dt/2);
            patch  = dif_phase(row_px, col_px);
            if isempty(patch); full_tile_p(fi,ti) = 0;
            else;              full_tile_p(fi,ti) = mean(patch(:)); end
        end
    end
    full_tile_p(isnan(full_tile_p)) = 0;
    fig_full = figure('Color', 'w', 'Position', [100 100 600 420]);
    tiledlayout(1, 1, 'Padding', 'compact');
    nexttile;
    imagesc(tile_t, tile_f, full_tile_p, [-0.5 0.5]);
    axis xy;
    colormap(cmap_bwr_full); colorbar;
    hold on; xline(0, 'k--', 'LineWidth', 0.8); hold off;
    ylim([5 35]);
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('Frequency (Hz)', 'FontSize', 12);
    title({['\bfFull ERSP Map — ' p_name], ...
           'Grand-average \DeltaERSP (dB, tiled) | All 400 tiles shown'}, ...
          'FontSize', 10);
    set(gca, 'FontSize', 11);
    out_path_full = fullfile(out_dir, 'all_in_one_pic', ['figure3_full_ersp_' p_name '.png']);
    exportgraphics(fig_full, out_path_full, 'Resolution', 300);
    fprintf('Saved: Full ERSP Map (%s)\n  → %s\n', p_name, out_path_full);
end

%% ------------------------------------------------------------------------- %
% 4. Foot Data Bottstrapping (Fixed Split Logic) %
% ------------------------------------------------------------------------- %
clear; close all; clc;
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
addpath('functions\for_GA');
config = configure_parameters();
out_dir = config.out_dir;

fprintf('\n=== STEP 4: Foot Bootstrapping (Control Fix) ===\n');

foot_cleaned = foot_cleaning(config);
foot_data_boot_combined = session_combine_foot(foot_cleaned);
clear foot_cleaned

fig_foot = figure('Color', 'w');
x = config.epoch_start:1/config.eeg_fs:config.epoch_end;
name_list = fields(foot_data_boot_combined);
foot_filtered = struct();
min_size = inf;

for sub_num = 1:length(name_list) 
    type_list = fields(foot_data_boot_combined.(name_list{sub_num}));
    for type_loop = 1:length(type_list) 
        data = foot_data_boot_combined.(name_list{sub_num}).(type_list{type_loop});
        temp_data = zeros(size(data));
        for i = 1:size(data, 1) 
            y = lowpass(data(i, :), config.foot_fpass(2), config.eeg_fs);
            y = y - mean(y(1, 1:1200));
            temp_data(i, :) = y;
        end

        if size(temp_data, 1) < min_size 
            min_size = size(temp_data, 1);
        end 
        foot_filtered.(name_list{sub_num}).(type_list{type_loop}) = temp_data;
    end 
end 
save(fullfile(out_dir, "foot", "foot_cleaned_session_combined_filtered.mat"), 'foot_filtered', '-v7.3');
boot_strap_pick_num_foot = min_size;

num_loops_foot = 10;
fprintf('Running %d foot bootstrap iterations entirely in memory...\n', num_loops_foot);
for sub_num = 1:length(name_list) 
    type_list = fields(foot_filtered.(name_list{sub_num}));

    sum_foot_epoch = struct();
    cached_ersp = struct();
    for type_loop = 1:length(type_list) 
        f_data = foot_filtered.(name_list{sub_num}).(type_list{type_loop});
        cached_ersp.(type_list{type_loop}) = ersp_foot(config, f_data);
    end

    for loop = 1:num_loops_foot
        foot_epoch_bootstraped = struct();
        for type_loop = 1:length(type_list) 
            foot_data_boot = foot_filtered.(name_list{sub_num}).(type_list{type_loop});
            ersp_data = cached_ersp.(type_list{type_loop});
            [n_trials, ~] = size(foot_data_boot);
            rng('shuffle');

            % --- FIXED LOGIC: Disjoint Split for Expected Condition ---
            if strcmp(type_list{type_loop}, 'epochs_expected')
                if n_trials >= 2 * boot_strap_pick_num_foot
                    idx_both = randperm(n_trials, 2 * boot_strap_pick_num_foot);
                    select_idx_dup = idx_both(1:boot_strap_pick_num_foot);
                    select_idx_exp = idx_both(boot_strap_pick_num_foot + 1:end);
                else 
                    half = floor(n_trials / 2);
                    idx_both = randperm(n_trials);
                    pool_dup = idx_both(1:half);
                    pool_exp = idx_both(half + 1:end);
                    select_idx_dup = pool_dup(randi(length(pool_dup), 1, boot_strap_pick_num_foot));
                    select_idx_exp = pool_exp(randi(length(pool_exp), 1, boot_strap_pick_num_foot));
                end

                [final_boot, final_std] = boot_strap_2d(config, foot_data_boot, boot_strap_pick_num_foot, select_idx_dup);
                [final_ersp, final_ersp_std] = boot_strap(config, ersp_data.all_trial_ersp_db, boot_strap_pick_num_foot, select_idx_dup);

                foot_epoch_bootstraped.expected_duplicated.final_boot = final_boot;
                foot_epoch_bootstraped.expected_duplicated.final_std = final_std;
                foot_epoch_bootstraped.expected_duplicated.final_ersp = final_ersp;
                foot_epoch_bootstraped.expected_duplicated.final_ersp_std = final_ersp_std;

                [final_boot, final_std] = boot_strap_2d(config, foot_data_boot, boot_strap_pick_num_foot, select_idx_exp);
                [final_ersp, final_ersp_std] = boot_strap(config, ersp_data.all_trial_ersp_db, boot_strap_pick_num_foot, select_idx_exp);

                foot_epoch_bootstraped.epochs_expected.final_boot = final_boot;
                foot_epoch_bootstraped.epochs_expected.final_std = final_std;
                foot_epoch_bootstraped.epochs_expected.final_ersp = final_ersp;
                foot_epoch_bootstraped.epochs_expected.final_ersp_std = final_ersp_std;
            else 
                select_idx = randperm(n_trials, boot_strap_pick_num_foot);
                [final_boot, final_std] = boot_strap_2d(config, foot_data_boot, boot_strap_pick_num_foot, select_idx);
                [final_ersp, final_ersp_std] = boot_strap(config, ersp_data.all_trial_ersp_db, boot_strap_pick_num_foot, select_idx);

                foot_epoch_bootstraped.(type_list{type_loop}).final_boot = final_boot;
                foot_epoch_bootstraped.(type_list{type_loop}).final_std = final_std;
                foot_epoch_bootstraped.(type_list{type_loop}).final_ersp = final_ersp;
                foot_epoch_bootstraped.(type_list{type_loop}).final_ersp_std = final_ersp_std;
            end 
        end

        b_types = fields(foot_epoch_bootstraped);
        for t = 1:length(b_types) 
            if loop == 1 
                sum_foot_epoch.(b_types{t}) = foot_epoch_bootstraped.(b_types{t});
            else 
                sum_foot_epoch.(b_types{t}).final_boot = sum_foot_epoch.(b_types{t}).final_boot + foot_epoch_bootstraped.(b_types{t}).final_boot;
                sum_foot_epoch.(b_types{t}).final_std = sum_foot_epoch.(b_types{t}).final_std + foot_epoch_bootstraped.(b_types{t}).final_std;
                sum_foot_epoch.(b_types{t}).final_ersp = sum_foot_epoch.(b_types{t}).final_ersp + foot_epoch_bootstraped.(b_types{t}).final_ersp;
                sum_foot_epoch.(b_types{t}).final_ersp_std = sum_foot_epoch.(b_types{t}).final_ersp_std + foot_epoch_bootstraped.(b_types{t}).final_ersp_std;
            end 
        end 
    end

    ave_foot_data = struct();
    for t = 1:length(b_types) 
        ave_foot_data.(b_types{t}).final_boot = sum_foot_epoch.(b_types{t}).final_boot / num_loops_foot;
        ave_foot_data.(b_types{t}).final_std = sum_foot_epoch.(b_types{t}).final_std / num_loops_foot;
        ave_foot_data.(b_types{t}).final_ersp = sum_foot_epoch.(b_types{t}).final_ersp / num_loops_foot;
        ave_foot_data.(b_types{t}).final_ersp_std = sum_foot_epoch.(b_types{t}).final_ersp_std / num_loops_foot;
    end 
    save(fullfile(out_dir, "foot2", "average_foot_data_bootstrapped_sub" + sub_num + "_.mat"), 'ave_foot_data', '-v7.3');
end

%% ------------------------------------------------------------------------- %
% 5. Format ERSP Tiles and Compute Subject - Level Differences %
% ------------------------------------------------------------------------- %
clear; close all; clc;
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
addpath('functions\for_GA');
config = configure_parameters();
out_dir = config.out_dir;

fprintf('\n=== STEP 5: Tiling and Subject Differences ===\n');

[time_range, freq_range] = range_compute(config);
ersp_boot_average = load(fullfile(out_dir, "all_in_one_pic", "avergae_bootstrap_each_sub.mat")).boot_average;
fields_boot = fieldnames(ersp_boot_average);
sub_tile = struct();
sub_tile_foot = struct();

for struct_num = 1:length(fields_boot) 
    current_name = fields_boot{struct_num};
    types = fieldnames(ersp_boot_average.(current_name));
    ave_foot_data = load(fullfile(out_dir, "foot2", "average_foot_data_bootstrapped_sub" + struct_num + "_.mat")).ave_foot_data;

    for type_num = 1:length(types) 
        current_type_name = types{type_num};
        sub_tile.(current_name).(current_type_name) = ersp_tile(ersp_boot_average.(current_name).(current_type_name).ersp, time_range, freq_range);

        current_type_name_foot = name_finder(current_type_name);
        sub_tile_foot.(current_name).(current_type_name_foot) = ersp_tile(ave_foot_data.(current_type_name_foot).final_ersp, time_range, freq_range);
    end 
end 
save(fullfile(out_dir, "ga_feature", "sub_tile.mat"), 'sub_tile', '-v7.3');
save(fullfile(out_dir, "ga_feature", "sub_tile_foot.mat"), 'sub_tile_foot', '-v7.3');


%% ------------------------------------------------------------------------- %
% 6. Model Training & Strict Validation (LOOCV Corrected) %
% ------------------------------------------------------------------------- %
clear; close all; clc;
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
addpath('functions\for_GA');
config = configure_parameters();
out_dir = config.out_dir;
% Make sure to load the previous dependencies
load(fullfile(out_dir, "ga_feature", "sub_tile.mat"), 'sub_tile');
load(fullfile(out_dir, "ga_feature", "sub_tile_foot.mat"), 'sub_tile_foot');
types = fieldnames(sub_tile.sub_1);
type_name_order = [2 1 3 2 5 4];
sub_list = fieldnames(sub_tile);
num_subs = length(sub_list);
all_sub_diffs = struct();
for i = 1:num_subs 
    sub_name = sub_list{i};
    for j = 1:3 
        eeg_type_A = types{type_name_order(2 * j - 1)};
        eeg_type_B = types{type_name_order(2 * j)};
        if j == 1 
            all_sub_diffs.(sub_name).t = sub_tile.(sub_name).(eeg_type_A).t;
            all_sub_diffs.(sub_name).f = sub_tile.(sub_name).(eeg_type_A).f;
        end
        foot_type_A = name_finder(eeg_type_A);
        foot_type_B = name_finder(eeg_type_B);
        all_sub_diffs.(sub_name).eeg.(string(eeg_type_A + "_vs_" + eeg_type_B)) = sub_tile.(sub_name).(eeg_type_A).vector - sub_tile.(sub_name).(eeg_type_B).vector;
        all_sub_diffs.(sub_name).foot.(string(eeg_type_A + "_vs_" + eeg_type_B)) = sub_tile_foot.(sub_name).(foot_type_A).vector - sub_tile_foot.(sub_name).(foot_type_B).vector;
    end 
end

fprintf('\n=== STEP 6a: LOOCV Computation (All 7 Methods) ===\n');
fprintf('Only re-run when changing computation logic. Results saved to disk.\n');

for i = 1:10
    deleteFig();
    load(fullfile(out_dir, "foot", "foot_cleaned_session_combined_filtered.mat"), "foot_filtered");
    [foot_features, foot_labels] = generate_foot_features(foot_filtered, config);
    loocv_save_path = LOOCV_Corrected(all_sub_diffs, foot_features, foot_labels, config);
    pause(1);
end

clear; close all; clc;
fprintf('\n=== STEP 6b: LOOCV Plotting (Loads Saved Results) ===\n');
fprintf('Self-contained: loads most recent results_LOOCV_*.mat automatically.\n');
LOOCV_Plot([]);

fprintf('\n=== PIPELINE DEMONSTRATION COMPLETE ===\n');

%% ------------------------------------------------------------------------- %
% 6c. Generate Supplementary Figure S2 (Methodology Breakdown)                %
% Run this section standalone after Step 6 has been run at least once.        %
% Loads all_sub_diffs + LOOCV results from disk — no recomputation needed.    %
% ------------------------------------------------------------------------- %
clear; close all; clc;
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
addpath('functions\for_GA');
config = configure_parameters();
out_dir = config.out_dir;

fprintf('\n=== STEP 6c: Generate Supplementary Figure S2 (Methodology Breakdown) ===\n');

% Load required data
load(fullfile(out_dir, "all_in_one_pic", "avergae_bootstrap_all_sub.mat"), 'ave_boot_average');
load(fullfile(out_dir, "ga_feature", "sub_tile.mat"), 'sub_tile');
load(fullfile(out_dir, "ga_feature", "sub_tile_foot.mat"), 'sub_tile_foot');
[time_range, freq_range] = range_compute(config);

% Reconstruct all_sub_diffs (needed for t-test significant map)
types_tile = fieldnames(sub_tile.sub_1);
type_name_order = [2 1 3 2 5 4];
sub_list = fieldnames(sub_tile);
num_subs = length(sub_list);
all_sub_diffs = struct();
for i = 1:num_subs
    sub_name = sub_list{i};
    for j = 1:3
        eeg_type_A = types_tile{type_name_order(2 * j - 1)};
        eeg_type_B = types_tile{type_name_order(2 * j)};
        if j == 1
            all_sub_diffs.(sub_name).t = sub_tile.(sub_name).(eeg_type_A).t;
            all_sub_diffs.(sub_name).f = sub_tile.(sub_name).(eeg_type_A).f;
        end
        all_sub_diffs.(sub_name).eeg.(string(eeg_type_A + "_vs_" + eeg_type_B)) = ...
            sub_tile.(sub_name).(eeg_type_A).vector - sub_tile.(sub_name).(eeg_type_B).vector;
    end
end

% Load the correct LOOCV results (Feb 24 15:32 — validated run)
loocv_file = fullfile(out_dir, 'ga_feature', 'results_LOOCV_2026-02-24_15-32.mat');
if ~exist(loocv_file, 'file')
    fprintf('ERROR: Correct LOOCV file not found: %s\n', loocv_file);
    return;
end
loocv_data = load(loocv_file);
fprintf('Loaded LOOCV results: results_LOOCV_2026-02-24_15-32.mat\n');

% GA consensus mask is stored in results.spearman.mask_x (1x1200)
best_mask_x = loocv_data.results.spearman.mask_x;

% Compute grand-average surprise difference
dif_surprise = ave_boot_average.unexpected.ersp - ave_boot_average.expected.ersp;

% Compute group-level t-test significance across all subjects (surprise contrast)
n_tiles_total = numel(sub_tile.sub_1.unexpected.vector);
eeg_diffs_matrix = zeros(num_subs, n_tiles_total);
for i = 1:num_subs
    sub_name = sub_list{i};
    eeg_diffs_matrix(i, :) = sub_tile.(sub_name).unexpected.vector(:)' - sub_tile.(sub_name).expected.vector(:)';
end
[~, ~, ~, ~] = ttest(eeg_diffs_matrix);  % dummy call to confirm ttest available
[h_tile, ~] = ttest(eeg_diffs_matrix);   % h=1 where p<0.05
t_vec = sub_tile.sub_1.unexpected.t;
f_vec = sub_tile.sub_1.unexpected.f;
n_freqs = length(f_vec);
n_times = length(t_vec);
sig_tile_mask = reshape(h_tile, n_freqs, n_times);

% Reconstruct full-resolution significance map by upsampling tiles to ERSP grid
[T_grid, F_grid] = meshgrid(time_range, freq_range);
sig_map = zeros(size(dif_surprise));
for fi = 1:n_freqs
    for ti = 1:n_times
        % Find ERSP pixels belonging to this tile
        t_lo = t_vec(ti) - (t_vec(2)-t_vec(1))/2;
        t_hi = t_vec(ti) + (t_vec(2)-t_vec(1))/2;
        f_lo = f_vec(fi) - (f_vec(2)-f_vec(1))/2;
        f_hi = f_vec(fi) + (f_vec(2)-f_vec(1))/2;
        row_mask = freq_range >= f_lo & freq_range <= f_hi;
        col_mask = time_range >= t_lo & time_range <= t_hi;
        if sig_tile_mask(fi, ti)
            sig_map(row_mask, col_mask) = dif_surprise(row_mask, col_mask);
        end
    end
end

% Generate the figure
fprintf('Generating Supplementary Figure S2: Full Methodology Breakdown...\n');
fig_supp_s2 = figure('Color', 'w', 'Position', [100 100 1600 500]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel (a): Raw grand-average ERSP (no masking)
visualisation_ersp(dif_surprise, fig_supp_s2, time_range, freq_range);
title("(a) Raw Grand-Average ERSP" + newline + "Unexpected - Expected", 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Frequency (Hz)'); xlabel('Time (s)');
ylim([5 40]); xlim([-2 3]);

% Panel (b): Group-level t-test significant regions (blue = significant ERD)
visualisation_ersp(sig_map, fig_supp_s2, time_range, freq_range);
title("(b) Group-Level Significance" + newline + "One-sample t-test, \alpha=0.05, N=25", 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Frequency (Hz)'); xlabel('Time (s)');
ylim([5 40]); xlim([-2 3]);

% Panel (c): GA consensus features (yellow circles overlay)
visualisation_ersp(dif_surprise, fig_supp_s2, time_range, freq_range);
hold on;
if ~isempty(best_mask_x)
    % EEG tiles are in positions 401-800 of best_mask_x (surprise contrast)
    eeg_mask_section = best_mask_x(401:800);
    [selected_t, selected_f] = mask_visualiser(eeg_mask_section, f_vec, t_vec);
    plot(selected_t, selected_f, 'yo', 'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
end
hold off;
title("(c) GA Consensus Features" + newline + "Yellow = \geq50% LOOCV folds", 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Frequency (Hz)'); xlabel('Time (s)');
ylim([5 40]); xlim([-2 3]);

out_path_s2 = fullfile(out_dir, "all_in_one_pic", "supplementary_figure_S2_methodology.png");
exportgraphics(fig_supp_s2, out_path_s2, 'Resolution', 300);
fprintf('Saved: Supplementary Figure S2 - Methodology Breakdown\n  → %s\n', out_path_s2);



%% ------------------------------------------------------------------------- %
%% ------------------------------------------------------------------------- %
% 7. Control Condition Analysis (Delta_control + GA)                          %
% Run this AFTER Step 6. Runs full GA + LOOCV on the null baseline contrast   %
% ------------------------------------------------------------------------- %
fprintf('\n=== STEP 7: Control Condition Analysis (GA + LOOCV) ===\n');

load(fullfile(out_dir, 'ga_feature', 'sub_tile.mat'), 'sub_tile');
load(fullfile(out_dir, 'ga_feature', 'sub_tile_foot.mat'), 'sub_tile_foot');
load(fullfile(out_dir, 'foot', 'foot_cleaned_session_combined_filtered.mat'), 'foot_filtered');
[foot_features, ~] = generate_foot_features(foot_filtered, config);

sub_list = fieldnames(sub_tile);
n_subs = length(sub_list);
n_tiles = 400;
tile_f = sub_tile.(sub_list{1}).expected_duplicated.f;

eeg_control = zeros(n_subs, n_tiles);
foot_control = zeros(n_subs, n_tiles);
for i = 1:n_subs
    sub_name = sub_list{i};
    eeg_control(i, :) = reshape(sub_tile.(sub_name).expected_duplicated.vector, 1, []) - reshape(sub_tile.(sub_name).expected.vector, 1, []);
    foot_control(i, :) = reshape(sub_tile_foot.(sub_name).expected_duplicated.vector, 1, []) - reshape(sub_tile_foot.(sub_name).epochs_expected.vector, 1, []);
end

rmv_idx = tile_f < 10 | tile_f > 30;
n_freqs = length(tile_f);
n_times = n_tiles / n_freqs;
for i = 1:n_subs
    block = reshape(eeg_control(i,:), n_freqs, n_times);
    block(rmv_idx, :) = 0;
    eeg_control(i,:) = block(:)';
end

eeg_control = normalize(eeg_control, 2);
foot_control = normalize(foot_features, 2); % Test against real foot features to see if null EEG predicts it

% Modify config to output specific files for control
config.save_filename_override = 'loocv_results_CONTROL.mat';

fprintf('\nRunning full LOOCV Pipeline on Control Contrast... This will take a few minutes.\n');
loocv_save_path_control = LOOCV_Corrected(eeg_control, foot_features, foot_labels, config);

fprintf('\n=== STEP 7a COMPLETE ===\n');

%% ------------------------------------------------------------------------- %
% 7b. Plot Control Condition Results                                          %
% ------------------------------------------------------------------------- %
fprintf('\n=== STEP 7b: Plot Control Condition ===\n');
if exist('loocv_save_path_control', 'var')
    LOOCV_Plot(loocv_save_path_control);
else
    loocv_save_path_control = fullfile(config.out_dir, 'ga_feature', 'loocv_results_CONTROL.mat');
    if exist(loocv_save_path_control, 'file')
        LOOCV_Plot(loocv_save_path_control);
    else
        fprintf('Control results not found. Run Step 7a first.\n');
    end
end
fprintf('\n=== STEP 7b COMPLETE ===\n');

%% -------------------------------------------------------------------------
% 8a. Supervisor Check — Issue 6: Per-Condition ERSP for Each Consensus Tile
% Standalone. Loads saved data from disk.
% For each GA consensus tile (>=50%), plots UNEXPECTED vs BASELINE and
% EXPECTED vs BASELINE separately to determine the source of the contrast:
%   Unexpected↑, Expected≈0 → surprise-driven elevation
%   Unexpected≈0, Expected↓ → expected events suppressed; surprise is "normal"
%   Both↑, Unexpected more  → graded response
% -------------------------------------------------------------------------
clear; close all; clc;
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
addpath('functions\for_GA');
config = configure_parameters();
out_dir = config.out_dir;

fprintf('\n=== STEP 8a: Issue 6 — Per-Condition ERSP per Consensus Tile ===\n');

load(fullfile(out_dir, 'all_in_one_pic', 'avergae_bootstrap_all_sub.mat'), 'ave_boot_average');
load(fullfile(out_dir, 'ga_feature', 'sub_tile.mat'), 'sub_tile');
[time_range, freq_range] = range_compute(config);

% Auto-detect all LOOCV files from today; loop over files and methods
ga_dir_8a   = fullfile(out_dir, 'ga_feature');
today_str   = datestr(now, 'yyyy-mm-dd');
today_8a    = dir(fullfile(ga_dir_8a, ['results_LOOCV_' today_str '_*.mat']));
if isempty(today_8a)
    error('No LOOCV files found for today (%s). Run Step 6 first.', today_str);
end
fprintf('Found %d LOOCV file(s) for today (%s).\n', length(today_8a), today_str);
str_cell_8a = {"spearman","pca","robust","pls","lasso","ridge","svr"};
for f8_idx = 1:length(today_8a)
    loocv_file_8a = fullfile(ga_dir_8a, today_8a(f8_idx).name);
    [~, loocv_name_8a, ~] = fileparts(today_8a(f8_idx).name);
    fprintf('\n--- File %d/%d: %s ---\n', f8_idx, length(today_8a), loocv_name_8a);
    loocv_data = load(loocv_file_8a, 'results');
for m8_idx = 1:length(str_cell_8a)
    method_name = str_cell_8a{m8_idx};
    fprintf('  Method: %s\n', method_name);

% Tile grid parameters
first_sub  = fieldnames(sub_tile); first_sub  = first_sub{1};
first_type = fieldnames(sub_tile.(first_sub)); first_type = first_type{1};
tile_t   = sub_tile.(first_sub).(first_type).t;
tile_f   = sub_tile.(first_sub).(first_type).f;
n_tile_f = length(tile_f);
n_tile_t = length(tile_t);
tile_dt  = tile_t(2) - tile_t(1);   % 0.2 s
tile_df  = tile_f(2) - tile_f(1);   % 3 Hz

% Consensus mask — surprise phase is tiles 401-800 in the 1200-tile vector
% Use ridge mask (same as Section 3c) — this gives the 5-tile result
mask_surprise   = loocv_data.results.(method_name).mask_x(401:800);
counts_surprise = loocv_data.results.(method_name).mask_x_counts(401:800);
n_tile_total    = length(mask_surprise);
assert(n_tile_total == n_tile_f * n_tile_t, ...
    'Mask length (%d) != n_tile_f*n_tile_t (%d)', n_tile_total, n_tile_f*n_tile_t);
mask_2d   = reshape(mask_surprise,   n_tile_f, n_tile_t);
counts_2d = reshape(counts_surprise, n_tile_f, n_tile_t);

sel_idx = find(mask_surprise == 1);
n_sel   = length(sel_idx);
fprintf('Consensus tiles found: %d\n', n_sel);

% Extract per-tile stats — sort by time for readability
t_centers = zeros(n_sel, 1);
f_centers = zeros(n_sel, 1);
unexp_vals    = zeros(n_sel, 1);
exp_vals      = zeros(n_sel, 1);
contrast_vals = zeros(n_sel, 1);
count_vals    = zeros(n_sel, 1);

for i = 1:n_sel
    ti_lin  = sel_idx(i);
    fi_tile = mod(ti_lin - 1, n_tile_f) + 1;
    tt_tile = floor((ti_lin - 1) / n_tile_f) + 1;

    f_ctr = tile_f(fi_tile);
    t_ctr = tile_t(tt_tile);

    row_px = freq_range >= (f_ctr - tile_df/2) & freq_range < (f_ctr + tile_df/2);
    col_px = time_range >= (t_ctr - tile_dt/2) & time_range < (t_ctr + tile_dt/2);

    t_centers(i)    = t_ctr;
    f_centers(i)    = f_ctr;
    unexp_vals(i)   = mean(ave_boot_average.unexpected.ersp(row_px, col_px), 'all');
    exp_vals(i)     = mean(ave_boot_average.expected.ersp(row_px, col_px),   'all');
    contrast_vals(i)= unexp_vals(i) - exp_vals(i);
    count_vals(i)   = counts_2d(fi_tile, tt_tile);
end

[~, sort_order] = sort(t_centers);
t_centers    = t_centers(sort_order);
f_centers    = f_centers(sort_order);
unexp_vals   = unexp_vals(sort_order);
exp_vals     = exp_vals(sort_order);
contrast_vals= contrast_vals(sort_order);
count_vals   = count_vals(sort_order);

% Print table
fprintf('\n%-22s  %10s  %10s  %10s\n', 'Tile', 'Unexpected', 'Expected', 'Contrast');
fprintf('%s\n', repmat('-', 1, 57));
for i = 1:n_sel
    fprintf('t=%+.1fs  f=%2dHz (%2d/25)  %+10.4f  %+10.4f  %+10.4f\n', ...
        t_centers(i), f_centers(i), count_vals(i), unexp_vals(i), exp_vals(i), contrast_vals(i));
end
fprintf('\nBaseline = 0 dB. Positive = elevated above baseline.\n');
fprintf('Interpretation guide:\n');
fprintf('  Unexpected↑ + Expected≈0  → surprise-driven elevation\n');
fprintf('  Unexpected≈0 + Expected↓  → expected suppression; surprise = default state\n');
fprintf('  Both positive, Unexp > Exp → graded response\n');

% Bar chart: grouped bars per tile
tick_labels = cell(n_sel, 1);
for i = 1:n_sel
    tick_labels{i} = sprintf('t=%+.1fs\nf=%dHz\n(%d/25)', t_centers(i), f_centers(i), count_vals(i));
end

fig_8a = figure('Color', 'w', 'Position', [100, 100, max(700, n_sel*140), 450]);
bar_data = [unexp_vals, exp_vals];
b = bar(1:n_sel, bar_data, 'grouped');
b(1).FaceColor = [0.85 0.20 0.20];
b(2).FaceColor = [0.20 0.40 0.85];
hold on;
yline(0, 'k-', 'LineWidth', 1.2);
hold off;
set(gca, 'XTick', 1:n_sel, 'XTickLabel', tick_labels, 'FontSize', 10);
ylabel('\DeltaERSP vs Baseline (dB)', 'FontSize', 12);
title({'Supervisor Check — Issue 6: Per-Condition ERSP for Each Consensus Tile', ...
       sprintf('File: %s | Method: %s', loocv_name_8a, method_name), ...
       'Red = Unexpected | Blue = Expected | Baseline = 0 dB'}, ...
      'FontSize', 11, 'FontWeight', 'bold');
legend({'Unexpected', 'Expected'}, 'Location', 'best');
grid on;

out_path_8a = fullfile(out_dir, 'all_in_one_pic', ...
    ['supervisor_check_issue6_percondition_' loocv_name_8a '_' method_name '.png']);
exportgraphics(fig_8a, out_path_8a, 'Resolution', 300);
fprintf('\nSaved: %s\n', out_path_8a);
end  % m8_idx methods loop
end  % f8_idx files loop
fprintf('=== STEP 8a COMPLETE ===\n');

%% -------------------------------------------------------------------------
% 8b. Supervisor Check — Issue 1: Frequency Spectrum at t≈0 (Doublet Check)
% Standalone. Loads grand-average ERSP from disk.
% Plots ΔERSP (Unexpected - Expected) as a function of frequency at t=0±0.1s.
% If one smooth hump at 21-25 Hz → tiles 2&3 are spectral smearing of one source.
% If two distinct peaks with a dip → doublet is more credible.
% Also reports exact tile frequency boundaries for the ~21 Hz and ~24 Hz tiles.
% -------------------------------------------------------------------------
clear; close all; clc;
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
addpath('functions\for_GA');
config = configure_parameters();
out_dir = config.out_dir;

fprintf('\n=== STEP 8b: Issue 1 — Frequency Spectrum at t≈0 (Doublet Check) ===\n');

load(fullfile(out_dir, 'all_in_one_pic', 'avergae_bootstrap_all_sub.mat'), 'ave_boot_average');
[time_range, freq_range] = range_compute(config);

dif_surprise = ave_boot_average.unexpected.ersp - ave_boot_average.expected.ersp;

tile_df = 3; tile_dt = 0.2;

% Report tile frequency boundaries
fprintf('\nTile boundary report (tile_df = %.0f Hz, tile_dt = %.2f s):\n', tile_df, tile_dt);
fprintf('  f=21 Hz tile: covers 21–24 Hz\n');
fprintf('  f=24 Hz tile: covers 24–27 Hz\n');
fprintf('  → Shared boundary at 24 Hz (no gap between tiles).\n');
fprintf('  → If spectrum shows one hump peaking ~21–25 Hz: tiles are spectral smearing.\n');
fprintf('  → If two peaks with a dip at ~24 Hz: doublet more credible.\n\n');

% Extract frequency slice at t = 0 ± half tile width (±0.1 s)
t_center = 0;
t_half   = tile_dt / 2;
col_idx  = time_range >= (t_center - t_half) & time_range <= (t_center + t_half);
freq_slice = mean(dif_surprise(:, col_idx), 2);

% ROI: 5–40 Hz
f_roi = freq_range >= 5 & freq_range <= 40;

fig_8b = figure('Color', 'w', 'Position', [100, 100, 640, 420]);
plot(freq_range(f_roi), freq_slice(f_roi), 'b-', 'LineWidth', 2);
hold on;
yline(0, 'k--', 'LineWidth', 0.8);
% Mark tile left edges (tile_f stores start of each tile, not centre)
% Tile at f=21 covers 21-24 Hz (centre 22.5 Hz); tile at f=24 covers 24-27 Hz (centre 25.5 Hz)
for f_left = [21, 24]
    xline(f_left, 'r--', 'LineWidth', 1.0, 'Alpha', 0.7);
    text(f_left, max(freq_slice(f_roi))*0.85, sprintf('f=%dHz\ntile left edge', f_left), ...
         'FontSize', 9, 'Color', [0.7 0 0], 'HorizontalAlignment', 'center');
end
% Mark shared tile boundary at 24 Hz (right edge of f=21 tile = left edge of f=24 tile)
xline(24, 'm-', 'LineWidth', 1.5, 'Alpha', 0.5, 'Label', 'shared boundary');
hold off;
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('\DeltaERSP at t\approx0 (dB)', 'FontSize', 12);
title({'Supervisor Check — Issue 1: Frequency Spectrum at t\approx0', ...
       'Dashed red = tile left edges (21 Hz, 24 Hz) | Magenta = shared boundary | One hump vs two peaks?'}, ...
      'FontSize', 11, 'FontWeight', 'bold');
xlim([5 40]);
grid on;

% Also print peak locations near 21-27 Hz
beta_roi = freq_range >= 18 & freq_range <= 30;
[peak_val, peak_idx] = max(freq_slice(beta_roi));
beta_freqs = freq_range(beta_roi);
fprintf('Peak in 18–30 Hz band: %.4f dB at %.1f Hz\n', peak_val, beta_freqs(peak_idx));

out_path_8b = fullfile(out_dir, 'all_in_one_pic', 'supervisor_check_issue1_freq_spectrum.png');
exportgraphics(fig_8b, out_path_8b, 'Resolution', 300);
fprintf('Saved: %s\n', out_path_8b);
fprintf('=== STEP 8b COMPLETE ===\n');

%% -------------------------------------------------------------------------
% 8c. Supervisor Check — Issue 5: Post-Stimulus Beta Time Course
% Standalone. Loads grand-average ERSP from disk.
% Plots ΔERSP averaged over 18–27 Hz from t=0 to +3 s.
% If single broad hump → tiles 4 & 5 are one rebound; drop two-phase claim.
% If trough between t≈1.2 s and t≈2.6 s → two phases are plausible.
% Note: tile 4 (f≈24 Hz) and tile 5 (f≈18 Hz) are at different freq ranges,
% so also plots them separately to show whether they are frequency-dissociable.
% -------------------------------------------------------------------------
clear; close all; clc;
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
addpath('functions\for_GA');
config = configure_parameters();
out_dir = config.out_dir;

fprintf('\n=== STEP 8c: Issue 5 — Post-Stimulus Beta Time Course ===\n');

load(fullfile(out_dir, 'all_in_one_pic', 'avergae_bootstrap_all_sub.mat'), 'ave_boot_average');
[time_range, freq_range] = range_compute(config);

dif_surprise = ave_boot_average.unexpected.ersp - ave_boot_average.expected.ersp;
unexp_ersp   = ave_boot_average.unexpected.ersp;
exp_ersp     = ave_boot_average.expected.ersp;

% Tile 4: f=24 Hz (covers 24-27 Hz), t≈+1.2 s
% Tile 5: f=18 Hz (covers 18-21 Hz), t≈+2.6 s
% Combined beta band: 18-27 Hz
t_lo = 0; t_hi = 3.0;
t_roi = time_range >= t_lo & time_range <= t_hi;

% Broad band: 18-27 Hz
f_broad = freq_range >= 18 & freq_range <= 27;
tc_contrast_broad = mean(dif_surprise(f_broad, :), 1);
tc_unexp_broad    = mean(unexp_ersp(f_broad,   :), 1);
tc_exp_broad      = mean(exp_ersp(f_broad,     :), 1);

% Tile-4 band only: 24-27 Hz
f_tile4 = freq_range >= 24 & freq_range < 27;
tc_contrast_t4 = mean(dif_surprise(f_tile4, :), 1);

% Tile-5 band only: 18-21 Hz
f_tile5 = freq_range >= 18 & freq_range < 21;
tc_contrast_t5 = mean(dif_surprise(f_tile5, :), 1);

fprintf('Tile 4 (upper beta): f=24 Hz, covers 24-27 Hz, expected peak at t≈+1.2 s\n');
fprintf('Tile 5 (mid-beta):   f=18 Hz, covers 18-21 Hz, expected peak at t≈+2.6 s\n');
fprintf('These are different frequency ranges — argument shifts from time-phase to frequency-band.\n\n');

fig_8c = figure('Color', 'w', 'Position', [100, 100, 1000, 420]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel (a): Broad-band contrast time course
nexttile;
plot(time_range(t_roi), tc_contrast_broad(t_roi), 'k-', 'LineWidth', 2);
hold on;
yline(0, 'k--', 'LineWidth', 0.8);
xline(1.2, 'r--', 'LineWidth', 1.2, 'Alpha', 0.8);
xline(2.6, 'b--', 'LineWidth', 1.2, 'Alpha', 0.8);
text(1.2, max(tc_contrast_broad(t_roi))*0.9, 'Tile4\nt=1.2s', ...
     'FontSize', 9, 'Color', [0.7 0 0], 'HorizontalAlignment', 'center');
text(2.6, max(tc_contrast_broad(t_roi))*0.9, 'Tile5\nt=2.6s', ...
     'FontSize', 9, 'Color', [0 0 0.7], 'HorizontalAlignment', 'center');
hold off;
xlabel('Time (s)', 'FontSize', 11);
ylabel('\DeltaERSP (dB)', 'FontSize', 11);
title({'(a) Contrast: 18–27 Hz', 'One hump = one rebound'}, 'FontSize', 10, 'FontWeight', 'bold');
xlim([t_lo t_hi]); grid on;

% Panel (b): Per-condition broad-band
nexttile;
plot(time_range(t_roi), tc_unexp_broad(t_roi), 'r-', 'LineWidth', 2); hold on;
plot(time_range(t_roi), tc_exp_broad(t_roi),   'b-', 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 0.8);
xline(1.2, 'r--', 'LineWidth', 0.8, 'Alpha', 0.5);
xline(2.6, 'b--', 'LineWidth', 0.8, 'Alpha', 0.5);
hold off;
xlabel('Time (s)', 'FontSize', 11);
ylabel('\DeltaERSP vs Baseline (dB)', 'FontSize', 11);
title({'(b) Per-Condition: 18–27 Hz', 'Red=Unexpected, Blue=Expected'}, 'FontSize', 10, 'FontWeight', 'bold');
legend({'Unexpected', 'Expected'}, 'Location', 'northeast', 'FontSize', 9);
xlim([t_lo t_hi]); grid on;

% Panel (c): Tile-specific bands (dissociation check)
nexttile;
plot(time_range(t_roi), tc_contrast_t4(t_roi), 'r-',  'LineWidth', 2); hold on;
plot(time_range(t_roi), tc_contrast_t5(t_roi), 'b-',  'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 0.8);
xline(1.2, 'k:', 'LineWidth', 0.8);
xline(2.6, 'k:', 'LineWidth', 0.8);
hold off;
xlabel('Time (s)', 'FontSize', 11);
ylabel('\DeltaERSP (dB)', 'FontSize', 11);
title({'(c) Tile-Specific Bands', 'Red=24-27Hz (Tile4), Blue=18-21Hz (Tile5)'}, ...
      'FontSize', 10, 'FontWeight', 'bold');
legend({'24-27 Hz (Tile4)', '18-21 Hz (Tile5)'}, 'Location', 'northeast', 'FontSize', 9);
xlim([t_lo t_hi]); grid on;

sgtitle({'Supervisor Check — Issue 5: Post-Stimulus Beta Time Course', ...
         'Single hump = one rebound; Trough between dashed lines = two phases possible'}, ...
        'FontSize', 11, 'FontWeight', 'bold');

out_path_8c = fullfile(out_dir, 'all_in_one_pic', 'supervisor_check_issue5_beta_timecourse.png');
exportgraphics(fig_8c, out_path_8c, 'Resolution', 300);
fprintf('Saved: %s\n', out_path_8c);
fprintf('=== STEP 8c COMPLETE ===\n');

%% -------------------------------------------------------------------------
% 8d. Supervisor Check — Issue 3: Alpha Band Tile Coverage (8–12 Hz)
% Standalone. Loads grand-average ERSP and per-subject tile data from disk.
% Reports: (1) which tiles cover 8–12 Hz, (2) grand mean and std of surprise
% contrast in those tiles.
%   High std >> |mean| → consistent signal present but cancelling across subs
%   Both ≈0            → genuinely absent alpha effect
% Also shows an ERSP heatmap restricted to the 6–15 Hz band.
% -------------------------------------------------------------------------
clear; close all; clc;
addpath('functions\for_obtainERD');
addpath('functions\for_ersp');
addpath('functions\for_foot');
addpath('functions\for_GA');
config = configure_parameters();
out_dir = config.out_dir;

fprintf('\n=== STEP 8d: Issue 3 — Alpha Band Tile Coverage (8-12 Hz) ===\n');

load(fullfile(out_dir, 'all_in_one_pic', 'avergae_bootstrap_all_sub.mat'), 'ave_boot_average');
load(fullfile(out_dir, 'ga_feature', 'sub_tile.mat'), 'sub_tile');
[time_range, freq_range] = range_compute(config);

% Tile grid
first_sub  = fieldnames(sub_tile); first_sub  = first_sub{1};
first_type = fieldnames(sub_tile.(first_sub)); first_type = first_type{1};
tile_t   = sub_tile.(first_sub).(first_type).t;
tile_f   = sub_tile.(first_sub).(first_type).f;
n_tile_t = length(tile_t);
tile_df  = tile_f(2) - tile_f(1);  % 3 Hz
tile_dt  = tile_t(2) - tile_t(1);  % 0.2 s

% Tiles overlapping 8–12 Hz: tile range [f, f+tile_df) must intersect [8, 12)
alpha_mask  = (tile_f + tile_df > 8) & (tile_f < 12);
alpha_f_idx = find(alpha_mask);
fprintf('\nTile frequency resolution: %.0f Hz per tile\n', tile_df);
fprintf('Tiles overlapping 8–12 Hz (alpha band):\n');
for k = alpha_f_idx'
    fprintf('  f=%dHz tile: covers %d-%d Hz\n', tile_f(k), tile_f(k), tile_f(k)+tile_df);
end
if isempty(alpha_f_idx)
    fprintf('  WARNING: No tiles overlap the 8–12 Hz alpha band!\n');
end

% Per-subject surprise contrast in alpha tiles
sub_list = fieldnames(sub_tile);
n_subs   = length(sub_list);
n_alpha  = length(alpha_f_idx);
alpha_sub_diff = zeros(n_subs, n_alpha * n_tile_t);

for s = 1:n_subs
    sub_name = sub_list{s};
    diff_2d  = sub_tile.(sub_name).unexpected.vector - sub_tile.(sub_name).expected.vector;
    % diff_2d is [n_tile_f × n_tile_t]; extract alpha-band rows
    alpha_sub_diff(s, :) = reshape(diff_2d(alpha_f_idx, :), 1, []);
end

grand_mean_alpha = mean(alpha_sub_diff, 1);
grand_std_alpha  = std(alpha_sub_diff, 0, 1);
snr_ratio        = mean(abs(grand_mean_alpha)) / (mean(grand_std_alpha) + eps);

fprintf('\nAlpha tile statistics (Unexpected - Expected, N=%d subjects):\n', n_subs);
fprintf('  Mean of |grand mean| across all alpha tiles: %.4f dB\n', mean(abs(grand_mean_alpha)));
fprintf('  Mean of std            across all alpha tiles: %.4f dB\n', mean(grand_std_alpha));
fprintf('  Signal-to-variability ratio: %.3f\n', snr_ratio);
fprintf('\n  Interpretation:\n');
fprintf('  - Ratio >> 1: consistent signal present (mean >> noise)\n');
fprintf('  - Ratio << 1: either genuinely absent, or high inter-subject variability\n');
fprintf('  - High std + low mean: consistent signal present but cancelling across subjects\n');

% Grand-average ERSP heatmap: 6–15 Hz band
dif_surprise = ave_boot_average.unexpected.ersp - ave_boot_average.expected.ersp;
f_vis_roi    = freq_range >= 6 & freq_range <= 15;
alpha_ersp   = dif_surprise(f_vis_roi, :);
alpha_freqs  = freq_range(f_vis_roi);

cmap_bwr = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ...
            ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];

fig_8d = figure('Color', 'w', 'Position', [100, 100, 900, 480]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(time_range, alpha_freqs, alpha_ersp, [-0.5 0.5]);
axis xy; colormap(cmap_bwr); colorbar;
hold on;
xline(0,  'k--', 'LineWidth', 1.2);
yline(8,  'w:',  'LineWidth', 1.5);
yline(12, 'w:',  'LineWidth', 1.5);
text(-1.8, 10.5, 'Alpha',   'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');
text(-1.8, 9.0, '8-12Hz', 'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');
hold off;
xlabel('Time (s)', 'FontSize', 11);
ylabel('Frequency (Hz)', 'FontSize', 11);
title({'(a) ERSP: 6–15 Hz', 'Unexpected − Expected'}, 'FontSize', 10, 'FontWeight', 'bold');
xlim([-2 3]);

nexttile;
x_idx = 1:length(grand_mean_alpha);
plot(x_idx, grand_mean_alpha, 'k-',  'LineWidth', 1.8); hold on;
plot(x_idx, grand_std_alpha,  'r--', 'LineWidth', 1.2);
yline(0, 'k:', 'LineWidth', 0.8);
hold off;
xlabel(sprintf('Alpha tile index (%d freq bands × %d time steps)', n_alpha, n_tile_t), 'FontSize', 10);
ylabel('\DeltaERSP (dB)', 'FontSize', 11);
title({'(b) Alpha Tile: Grand Mean vs Std', ...
       'Black=mean | Red=std | std>>mean → variability; both≈0 → absent'}, ...
      'FontSize', 10, 'FontWeight', 'bold');
legend({'Grand mean', 'Std dev'}, 'Location', 'best', 'FontSize', 9);
grid on;

sgtitle({'Supervisor Check — Issue 3: Alpha Band Coverage (8–12 Hz)', ...
         sprintf('SNR ratio = %.3f  |  Dotted white lines = 8-12 Hz alpha window', snr_ratio)}, ...
        'FontSize', 11, 'FontWeight', 'bold');

out_path_8d = fullfile(out_dir, 'all_in_one_pic', 'supervisor_check_issue3_alpha_coverage.png');
exportgraphics(fig_8d, out_path_8d, 'Resolution', 300);
fprintf('\nSaved: %s\n', out_path_8d);
fprintf('=== STEP 8d COMPLETE ===\n');
