%% plot_topomap_surprise_erd.m
% Generates a scalp topographic map of the Surprise Beta ERD
% (Unexpected - Expected, 13-30 Hz, 0-1s post-trigger), averaged across
% all subjects/sessions.
%
% NEW STANDALONE SCRIPT — does not modify any existing files.
% Run from: codes_improved
%
% Output: figures\topomap_surprise_erd.png

clear; clc;
addpath('functions\for_obtainERD');
config = configure_parameters();

%% --- Parameters ---
out_dir         = config.out_dir;
ersp_dir        = fullfile(out_dir, 'ersp_ave');
script_dir      = fileparts(mfilename('fullpath'));
save_path       = fullfile(script_dir, '..', 'paper', 'figures', 'topomap_surprise_erd.png');

beta_range      = [13 30];   % Hz  (full beta band as reported in paper)
analysis_window = [0.0 1.0]; % seconds post-trigger
n_ch            = 9;         % channels of interest used in the pipeline

% -----------------------------------------------------------------------
% Channel mapping (confirmed from headset layout diagram):
%   config.coi = [13  9  5  15  12  1  14  11  10]
%   Ch1(13)=FC1  Ch2(9)=FCz  Ch3(5)=FC2
%   Ch4(15)=C1   Ch5(12)=Cz  Ch6(1)=C2
%   Ch7(14)=CP1  Ch8(11)=CPz Ch9(10)=CP2
%
% Coordinates are standard 10-20 positions normalised to unit head radius.
% -----------------------------------------------------------------------
ch_labels = {'FC1','FCz','FC2','C1','Cz','C2','CP1','CPz','CP2'};
ch_x = [-0.22,  0.00,  0.22, -0.25,  0.00,  0.25, -0.22,  0.00,  0.22];
ch_y = [ 0.27,  0.32,  0.27,  0.00,  0.00,  0.00, -0.27, -0.32, -0.27];

%% --- Load ERSP data, average WITHIN participant, exclude #26 (N=25 cohort) ---
% Participant 26 is NOT in the N=25 GA cohort: diag_fingerprint.m confirms
% boot_average sub_1..sub_25 map 1:1 to PIDs 1..25 (r >= 0.99, 2nd <= 0.43),
% leaving #26 out. We mirror the diagnostic's file->participant mapping
% (sort files numerically, index via subject_number_list) and average sessions
% within each participant before the grand average, matching the subject-level
% averaging used for the reported CPz/Cz values (-0.94 / -0.38 dB).
EXCLUDE_PID = 26;

S   = load(fullfile(out_dir, 'subject_number_list.mat'), 'subject_number_list');
snl = S.subject_number_list(:);
pid_list = unique(snl);
pid_list = pid_list(pid_list ~= EXCLUDE_PID);   % N=25 cohort

files = dir(fullfile(ersp_dir, 'ersp_epoch_ave_sub_*.mat'));
if isempty(files)
    error('No ERSP data found in: %s\nCheck config.out_dir in configure_parameters.m', ersp_dir);
end
fn = arrayfun(@(x) sscanf(x.name, 'ersp_epoch_ave_sub_%d.mat'), files);
[~, si] = sort(fn);  files = files(si);   % numeric order to align with subject_number_list

% Per-participant accumulation of per-channel beta ERD (sessions averaged within PID)
part_sum = zeros(numel(pid_list), n_ch);
part_cnt = zeros(numel(pid_list), n_ch);

for f = 1:length(files)
    if f > numel(snl), break; end
    pid = snl(f);
    if pid == EXCLUDE_PID
        fprintf('Excluding %s (participant %d, not in N=25 cohort).\n', files(f).name, pid);
        continue;
    end
    r = find(pid_list == pid, 1);

    d = load(fullfile(files(f).folder, files(f).name));  % loads variable 'ersp_epoch'
    if ~isfield(d, 'ersp_epoch')
        fprintf('Skipping %s — no ersp_epoch variable found.\n', files(f).name);
        continue;
    end

    sessions = fieldnames(d.ersp_epoch);
    for s = 1:length(sessions)
        sname = sessions{s};
        if ~isfield(d.ersp_epoch.(sname), 'unexpected') || ...
           ~isfield(d.ersp_epoch.(sname), 'expected')
            continue;
        end

        ersp_unexp = d.ersp_epoch.(sname).unexpected.data_ersp;
        ersp_exp   = d.ersp_epoch.(sname).expected.data_ersp;
        if ~isfield(ersp_unexp, 'freq') || ~isfield(ersp_unexp, 'time')
            continue;
        end

        freqs    = ersp_unexp.freq;
        times    = ersp_unexp.time;
        beta_idx = freqs >= beta_range(1) & freqs <= beta_range(2);
        time_idx = times >= analysis_window(1) & times <= analysis_window(2);
        if sum(beta_idx) == 0 || sum(time_idx) == 0
            continue;
        end

        for ch = 1:n_ch
            field_name = sprintf('all_trial_ersp_db_single_ch_%d', ch);
            if ~isfield(ersp_unexp, field_name) || ~isfield(ersp_exp, field_name)
                continue;
            end
            % Each field is [trials x freq x time]; average over trials
            map_unexp = squeeze(mean(ersp_unexp.(field_name), 1));  % [freq x time]
            map_exp   = squeeze(mean(ersp_exp.(field_name),   1));  % [freq x time]
            diff_map  = map_unexp - map_exp;                        % Surprise signal
            val       = mean(mean(diff_map(beta_idx, time_idx)));
            if ~isnan(val)
                part_sum(r, ch) = part_sum(r, ch) + val;
                part_cnt(r, ch) = part_cnt(r, ch) + 1;
            end
        end
    end
end

% Within-participant mean, then grand mean across the N=25 participants
part_avg = part_sum ./ max(part_cnt, 1);
part_avg(part_cnt == 0) = NaN;
n_used = sum(any(part_cnt > 0, 2));
if n_used == 0
    error('No valid participant data loaded. Check data paths and field naming.');
end
fprintf('\nIncluded %d participants (target N=25).\n', n_used);
grand_avg = mean(part_avg, 1, 'omitnan');  % [1 x 9]
fprintf('Grand-average beta ERD per channel (dB):\n');
for c = 1:n_ch
    fprintf('  %s: %.4f\n', ch_labels{c}, grand_avg(c));
end

%% --- Topographic Plot ---
fig = figure('Color', 'w', 'Position', [100 100 540 540]);

head_r = 0.62;
theta  = linspace(0, 2*pi, 360);

% Create interpolation grid
grid_n = 220;
[xi, yi] = meshgrid(linspace(-head_r, head_r, grid_n), ...
                    linspace(-head_r, head_r, grid_n));
inside = xi.^2 + yi.^2 <= (head_r * 0.97)^2;

% Scattered natural-neighbour interpolation
F  = scatteredInterpolant(ch_x', ch_y', grand_avg', 'natural', 'none');
zi = F(xi, yi);
zi(~inside) = NaN;

% Symmetric colour limits
clim_val = max(abs(grand_avg)) * 1.15;
if clim_val < 0.01, clim_val = 0.5; end

% Filled contour surface
contourf(xi, yi, zi, 60, 'LineStyle', 'none');
hold on;
colormap(redblue_cmap(256));
clim([-clim_val, clim_val]);
cb = colorbar('eastoutside');
cb.Label.String = '\DeltaERSP (dB)  [Unexpected - Expected]';
cb.Label.FontSize = 10;
cb.FontSize = 9;

% Head outline
plot(head_r * cos(theta), head_r * sin(theta), 'k-', 'LineWidth', 2.2);

% Nose (at top, positive-y)
nose_x = [-0.07, 0,  0.07];
nose_y = [head_r * 0.95, head_r + 0.08, head_r * 0.95];
plot(nose_x, nose_y, 'k-', 'LineWidth', 2.2);

% Left ear (negative-x side)
ear_ang_L = linspace(pi - 0.35, pi + 0.35, 30);
plot(head_r * cos(ear_ang_L) - 0.02, head_r * sin(ear_ang_L), 'k-', 'LineWidth', 2.2);

% Right ear (positive-x side)
ear_ang_R = linspace(-0.35, 0.35, 30);
plot(head_r * cos(ear_ang_R) + 0.02, head_r * sin(ear_ang_R), 'k-', 'LineWidth', 2.2);

% Electrode markers
scatter(ch_x, ch_y, 60, 'k', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);

% Labels (offset slightly so they don't overlap the dot)
label_offsets_x = [-0.08, 0.00,  0.06, -0.10,  0.00,  0.06, -0.08,  0.00,  0.06];
label_offsets_y = [ 0.05, 0.05,  0.05,  0.05,  0.05,  0.05, -0.06, -0.06, -0.06];
for c = 1:n_ch
    text(ch_x(c) + label_offsets_x(c), ch_y(c) + label_offsets_y(c), ...
         ch_labels{c}, 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k', ...
         'HorizontalAlignment', 'center');
end

axis equal off;
title({'\bfSurprise Beta ERD — Scalp Topography', ...
       sprintf('Unexpected − Expected  |  %d–%d Hz  |  %.1f–%.1f s post-trigger  |  N=25', ...
               beta_range(1), beta_range(2), analysis_window(1), analysis_window(2))}, ...
      'FontSize', 11, 'FontWeight', 'bold');

% Save high-resolution PNG
exportgraphics(fig, save_path, 'Resolution', 300);
fprintf('\nFigure saved to: %s\n', save_path);

%% --- Helper: diverging red-blue colormap (blue=ERD, red=ERS) ---
function cmap = redblue_cmap(n)
    half = floor(n / 2);
    r = [linspace(0, 1, half), ones(1, n - half)];
    g = [linspace(0, 1, half), linspace(1, 0, n - half)];
    b = [ones(1, half), linspace(1, 0, n - half)];
    cmap = [r(:), g(:), b(:)];
end
