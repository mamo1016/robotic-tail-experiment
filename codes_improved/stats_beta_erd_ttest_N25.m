%% stats_beta_erd_ttest_N25.m
% Re-computes the group-level one-sample t-test on the Surprise Beta ERD
% (Unexpected - Expected, 13-30 Hz, 0-1 s post-trigger) per channel,
% restricted to the N = 25 EEG cohort (poor-EEG participant excluded), so
% that N matches the GA-LOOCV coupling analysis.
%
% WHY THIS EXISTS
%   stats_beta_erd_ttest.m reads ersp_ave/ (pre-exclusion) and reports
%   N = 26 (df = 25). The poor-EEG participant was excluded downstream,
%   when ga_feature/sub_tile.mat was built (N = 25). This script restricts
%   the t-test to the same participant set sub_tile represents -> N = 25.
%
% SAFETY / AUDIT
%   - Writes ONLY new CSVs. Never overwrites the original N=26 CSV or data.
%   - Step A reproduces the N=26 result so it can be checked against the
%     existing CSV (validates the ersp_ave reading is faithful).
%   - Per-participant CPz values are printed with participant IDs so the
%     excluded participant is visible before the N=25 t-test is computed.
%
% Run from: codes_improved   (MATLAB R2026a)

clear; clc;
addpath('functions\for_obtainERD');
config = configure_parameters();

%% --- Parameters (identical to stats_beta_erd_ttest.m) ---
out_dir         = config.out_dir;
ersp_dir        = fullfile(out_dir, 'ersp_ave');

beta_range      = [13 30];   % Hz
analysis_window = [0.0 1.0]; % seconds post-trigger
n_ch            = 9;

ch_labels = {'FC1','FCz','FC2','C1','Cz','C2','CP1','CPz','CP2'};
cpz_ch    = find(strcmp(ch_labels, 'CPz'));   % column index for reporting

%% --- Load session -> participant mapping ---
list_file = fullfile(out_dir, 'subject_number_list.mat');
if ~exist(list_file, 'file')
    error('subject_number_list.mat not found in: %s', out_dir);
end
tmp = load(list_file, 'subject_number_list');
subject_number_list = tmp.subject_number_list;      % length = n_sessions
unique_num = unique(subject_number_list);            % sorted participant IDs

fprintf('subject_number_list: %d sessions, %d unique participant IDs.\n', ...
    numel(subject_number_list), numel(unique_num));
fprintf('  participant IDs: %s\n', mat2str(unique_num(:)'));

%% --- Derive the EEG keep-set from sub_tile.mat (positional -> ID map) ---
subtile_file = fullfile(out_dir, 'ga_feature', 'sub_tile.mat');
if ~exist(subtile_file, 'file')
    error('sub_tile.mat not found in: %s', fullfile(out_dir, 'ga_feature'));
end
st = load(subtile_file, 'sub_tile');
st_fields = fieldnames(st.sub_tile);                           % {'sub_1',...}
st_pos    = sort(cellfun(@(x) sscanf(x, 'sub_%d'), st_fields)); % positional keys
fprintf('\nsub_tile: %d entries, positional keys sub_%d..sub_%d.\n', ...
    numel(st_pos), min(st_pos), max(st_pos));

% Map each positional key sub_K -> participant ID unique_num(K).
% Valid only if subject_number_list is unchanged since sub_tile was built.
% Validated below by the N=26 reproduction matching the existing CSV.
if max(st_pos) > numel(unique_num)
    error(['sub_tile positional key (%d) exceeds number of participant IDs ' ...
           '(%d). subject_number_list may have changed since sub_tile was ' ...
           'built; cannot map positional->ID safely.'], ...
           max(st_pos), numel(unique_num));
end
keep_ids = sort(reshape(unique_num(st_pos), 1, []));   % participant IDs in EEG set
excluded_ids = setdiff(reshape(unique_num, 1, []), keep_ids);
fprintf('Mapped EEG keep-set (N=%d): %s\n', numel(keep_ids), mat2str(keep_ids));
fprintf('Excluded participant ID(s): %s\n', mat2str(excluded_ids));

%% --- Load per-session beta ERD, group by participant ID (ALL participants) ---
files = dir(fullfile(ersp_dir, 'ersp_epoch_ave_sub_*.mat'));
if isempty(files)
    error('No ERSP data found in: %s', ersp_dir);
end
fnames = {files.name};
fnums  = cellfun(@(x) sscanf(x, 'ersp_epoch_ave_sub_%d.mat'), fnames);
[~, sort_idx] = sort(fnums);
files = files(sort_idx);

if length(files) ~= length(subject_number_list)
    warning('File count (%d) does not match subject_number_list length (%d).', ...
        length(files), length(subject_number_list));
end

sub_data = containers.Map('KeyType','char','ValueType','any');  % pid -> [sessions x ch]
for f = 1:length(files)
    if f > length(subject_number_list)
        fprintf('No participant ID for session %d -- skipping.\n', f);
        continue;
    end
    sub_id = num2str(subject_number_list(f));   % participant ID as string

    d = load(fullfile(ersp_dir, files(f).name));
    if ~isfield(d, 'ersp_epoch'), continue; end

    sessions = fieldnames(d.ersp_epoch);
    for s = 1:length(sessions)
        sname = sessions{s};
        if ~isfield(d.ersp_epoch.(sname), 'unexpected') || ...
           ~isfield(d.ersp_epoch.(sname), 'expected')
            continue;
        end
        ersp_unexp = d.ersp_epoch.(sname).unexpected.data_ersp;
        ersp_exp   = d.ersp_epoch.(sname).expected.data_ersp;
        if ~isfield(ersp_unexp, 'freq') || ~isfield(ersp_unexp, 'time'), continue; end

        freqs    = ersp_unexp.freq;
        times    = ersp_unexp.time;
        beta_idx = freqs >= beta_range(1) & freqs <= beta_range(2);
        time_idx = times >= analysis_window(1) & times <= analysis_window(2);
        if sum(beta_idx)==0 || sum(time_idx)==0, continue; end

        ch_vals = NaN(1, n_ch);
        for ch = 1:n_ch
            field_name = sprintf('all_trial_ersp_db_single_ch_%d', ch);
            if ~isfield(ersp_unexp, field_name) || ~isfield(ersp_exp, field_name)
                continue;
            end
            map_unexp = squeeze(mean(ersp_unexp.(field_name), 1));
            map_exp   = squeeze(mean(ersp_exp.(field_name),   1));
            diff_map  = map_unexp - map_exp;
            ch_vals(ch) = mean(mean(diff_map(beta_idx, time_idx)));
        end

        if isKey(sub_data, sub_id)
            sub_data(sub_id) = [sub_data(sub_id); ch_vals];
        else
            sub_data(sub_id) = ch_vals;
        end
    end
end

%% --- Average within subject; keep participant IDs aligned to rows ---
sub_ids_str = keys(sub_data);
pid_vec     = cellfun(@str2double, sub_ids_str);
[pid_vec, order] = sort(pid_vec);
sub_ids_str = sub_ids_str(order);
n_all       = numel(pid_vec);

sub_means = NaN(n_all, n_ch);
for i = 1:n_all
    sub_means(i, :) = mean(sub_data(sub_ids_str{i}), 1, 'omitnan');
end

%% --- Per-participant CPz table (so the excluded subject is visible) ---
fprintf('\n=== Per-participant CPz beta-ERD (all %d found) ===\n', n_all);
fprintf('%-13s  %10s  %9s\n', 'ParticipantID', 'CPz(dB)', 'InEEGset');
fprintf('%s\n', repmat('-', 1, 36));
for i = 1:n_all
    inset = ismember(pid_vec(i), keep_ids);
    fprintf('%-13d  %10.4f  %9s\n', pid_vec(i), sub_means(i, cpz_ch), string(inset));
end

%% --- Step A: reproduce N=26 (validation vs existing CSV) ---
mask_all = true(n_all, 1);
ttest_table('STEP A  (reproduce, all participants)', mask_all, ...
    sub_means, ch_labels, n_ch, ...
    fullfile(out_dir, 'stats_beta_erd_ttest_results_REPRO_all.csv'));

%% --- Step B: N=25 (EEG keep-set only) ---
mask_keep = ismember(pid_vec(:), keep_ids(:));
ttest_table('STEP B  (N=25 EEG keep-set)', mask_keep, ...
    sub_means, ch_labels, n_ch, ...
    fullfile(out_dir, 'stats_beta_erd_ttest_results_N25.csv'));

fprintf('\nBonferroni threshold for 9 channels: p < %.4f\n', 0.05/9);
fprintf('Excluded participant ID(s) from N=25: %s\n', mat2str(excluded_ids));


%% ===================== local functions (must be at end) =================
function ttest_table(tag, mask, sub_means, ch_labels, n_ch, out_csv)
    vals_mask = sub_means(mask, :);
    fprintf('\n=== %s : N = %d subjects ===\n', tag, sum(mask));
    fprintf('%-6s  %8s  %8s  %8s  %8s  %8s  %4s\n', ...
        'Chan','Mean(dB)','SD','t','p','Cohen_d','df');
    fprintf('%s\n', repmat('-', 1, 64));
    fid = fopen(out_csv, 'w');
    if fid == -1
        warning('Could not open %s for writing; printing only.', out_csv);
    else
        fprintf(fid, 'Channel,Mean_dB,SD_dB,t,p,Cohen_d,N,df\n');
    end
    for ch = 1:n_ch
        v = vals_mask(:, ch);
        v = v(~isnan(v));
        n = numel(v);  m = mean(v);  sd = std(v);
        se = sd / sqrt(n);  t = m / se;
        p = 2 * tcdf(-abs(t), n-1);  dcoh = m / sd;
        fprintf('%-6s  %8.4f  %8.4f  %8.3f  %8.4f  %8.3f  %4d\n', ...
            ch_labels{ch}, m, sd, t, p, dcoh, n-1);
        if fid ~= -1
            fprintf(fid, '%s,%.4f,%.4f,%.3f,%.4f,%.3f,%d,%d\n', ...
                ch_labels{ch}, m, sd, t, p, dcoh, n, n-1);
        end
    end
    if fid ~= -1
        fclose(fid);
        fprintf('Saved: %s\n', out_csv);
    end
end
