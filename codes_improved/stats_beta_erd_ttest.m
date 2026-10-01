%% stats_beta_erd_ttest.m
% Computes group-level one-sample t-test on the Surprise Beta ERD
% (Unexpected - Expected, 13-30 Hz, 0-1s post-trigger) per channel.
%
% Loads subject_number_list.mat to correctly map session files (1-77) to
% participant IDs (1-26), averages sessions within each participant, then
% runs a one-sample t-test against 0 for each channel.
%
% NOTE: ersp_ave/ contains one file per SESSION (77 total).
%       subject_number_list.mat maps session index -> participant ID.
%       This script groups sessions by participant before averaging.
%
% NEW STANDALONE SCRIPT — does not modify any existing files.
% Run from: codes_improved
%
% Output: printed t-stat, p-value, Cohen's d per channel (for paper text)

clear; clc;
addpath('functions\for_obtainERD');
config = configure_parameters();

%% --- Parameters ---
out_dir         = config.out_dir;
ersp_dir        = fullfile(out_dir, 'ersp_ave');

beta_range      = [13 30];   % Hz
analysis_window = [0.0 1.0]; % seconds post-trigger
n_ch            = 9;

ch_labels = {'FC1','FCz','FC2','C1','Cz','C2','CP1','CPz','CP2'};

%% --- Load session->participant mapping ---
list_file = fullfile(out_dir, 'subject_number_list.mat');
if ~exist(list_file, 'file')
    error('subject_number_list.mat not found in: %s', out_dir);
end
tmp = load(list_file, 'subject_number_list');
subject_number_list = tmp.subject_number_list;  % length = n_sessions (77)

%% --- Load per-session beta ERD, group by participant ID ---
files = dir(fullfile(ersp_dir, 'ersp_epoch_ave_sub_*.mat'));
if isempty(files)
    error('No ERSP data found in: %s', ersp_dir);
end

% Sort files numerically (sub_1, sub_2, ..., sub_77) to match list order
fnames = {files.name};
fnums  = cellfun(@(x) sscanf(x, 'ersp_epoch_ave_sub_%d.mat'), fnames);
[~, sort_idx] = sort(fnums);
files = files(sort_idx);

if length(files) ~= length(subject_number_list)
    warning('File count (%d) does not match subject_number_list length (%d).', ...
        length(files), length(subject_number_list));
end

% Store: participant_id (string) -> [sessions x channels] beta ERD
sub_data = containers.Map('KeyType','char','ValueType','any');

for f = 1:length(files)
    fname = files(f).name;

    % Use participant ID from subject_number_list (not file number)
    if f > length(subject_number_list)
        fprintf('No participant ID for session %d — skipping.\n', f);
        continue;
    end
    sub_id = num2str(subject_number_list(f));  % true participant ID

    d = load(fullfile(ersp_dir, fname));
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

%% --- Average within subject ---
sub_ids  = keys(sub_data);
n_subs   = length(sub_ids);
sub_means = NaN(n_subs, n_ch);

for i = 1:n_subs
    sub_means(i, :) = mean(sub_data(sub_ids{i}), 1, 'omitnan');
end

fprintf('\nN = %d subjects loaded.\n\n', n_subs);
fprintf('%-6s  %8s  %8s  %8s  %8s  %8s\n', ...
    'Chan', 'Mean(dB)', 'SD', 't', 'p', 'Cohen_d');
fprintf('%s\n', repmat('-', 1, 58));

for ch = 1:n_ch
    vals = sub_means(:, ch);
    vals = vals(~isnan(vals));
    n    = length(vals);
    m    = mean(vals);
    sd   = std(vals);
    se   = sd / sqrt(n);
    t    = m / se;
    p    = 2 * tcdf(-abs(t), n-1);   % two-tailed
    d    = m / sd;                   % Cohen's d vs 0
    fprintf('%-6s  %8.4f  %8.4f  %8.3f  %8.4f  %8.3f\n', ...
        ch_labels{ch}, m, sd, t, p, d);
end

fprintf('\nNote: use Bonferroni correction threshold p < %.4f for 9 channels.\n', 0.05/9);

%% --- Save results to CSV ---
out_file = fullfile(out_dir, 'stats_beta_erd_ttest_results.csv');
fid = fopen(out_file, 'w');
fprintf(fid, 'Channel,Mean_dB,SD_dB,t,p,Cohen_d,N\n');
for ch = 1:n_ch
    vals = sub_means(:, ch);
    vals = vals(~isnan(vals));
    n    = length(vals);
    m    = mean(vals);
    sd   = std(vals);
    se   = sd / sqrt(n);
    t    = m / se;
    p    = 2 * tcdf(-abs(t), n-1);
    d    = m / sd;
    fprintf(fid, '%s,%.4f,%.4f,%.3f,%.4f,%.3f,%d\n', ch_labels{ch}, m, sd, t, p, d, n);
end
fclose(fid);
fprintf('\nResults saved to: %s\n', out_file);
