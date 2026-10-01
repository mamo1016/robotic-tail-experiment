%% bayesian_equivalence_learning.m
% Bayesian equivalence test for the temporal evolution of the postural
% response across unexpected trials.
%
% Addresses Reviewer MC-3: "p=0.241 does not prove rapid learning."
% This script formally tests whether the observed linear slope is
% equivalent to zero using:
%   (1) TOST (Two One-Sided Tests) — frequentist equivalence
%   (2) Bayes Factor (BF01) — how much more likely is H0 (flat) vs H1
%
% Output: printed results + figure saved to paper/figures/
%
% NEW STANDALONE SCRIPT — does not modify any existing files.

clear; clc;
addpath('functions\for_obtainERD');
config = configure_parameters();
out_dir = config.out_dir;

paper_fig_dir = fullfile(fileparts(pwd), 'paper', 'figures');

%% --- Load temporal evolution data ---
% Expected file: temporal_evolution_data.mat OR embedded in LOOCV results
% Try common output locations
candidate_files = {
    fullfile(out_dir, 'temporal_evolution_data.mat'),
    fullfile(out_dir, 'ga_feature', 'temporal_evolution_data.mat'),
    fullfile(out_dir, 'all_in_one_pic', 'temporal_evolution_data.mat')
};

data_loaded = false;
for ci = 1:length(candidate_files)
    if exist(candidate_files{ci}, 'file')
        load(candidate_files{ci});
        data_loaded = true;
        fprintf('Loaded: %s\n', candidate_files{ci});
        break;
    end
end

if ~data_loaded
    % Fallback: try to reconstruct from LOOCV results
    ga_dir = fullfile(out_dir, 'ga_feature');
    loocv_files = dir(fullfile(ga_dir, 'results_LOOCV_*.mat'));
    if isempty(loocv_files)
        error(['Temporal evolution data not found. Please check output directory: %s\n' ...
               'Expected file: temporal_evolution_data.mat'], out_dir);
    end
    [~, idx] = max([loocv_files.datenum]);
    load(fullfile(loocv_files(idx).folder, loocv_files(idx).name), 'results');
    % Extract per-subject CoP RMS if stored in results
    if isfield(results, 'temporal_cop_rms')
        cop_rms_per_trial = results.temporal_cop_rms; % [N_subs x N_trials]
    else
        error(['Cannot find temporal evolution data.\n' ...
               'Please run the main pipeline first to generate temporal_evolution_data.mat']);
    end
end

%% --- Compute per-subject linear slopes ---
% cop_rms_per_trial: [N_subjects x N_unexpected_trials]
% Each row = one subject's CoP RMS across the sequence of unexpected trials

[N_subs, N_trials] = size(cop_rms_per_trial);
trial_seq = (1:N_trials)';

slopes = NaN(N_subs, 1);
for s = 1:N_subs
    y = cop_rms_per_trial(s, :)';
    valid = ~isnan(y);
    if sum(valid) > 2
        p = polyfit(trial_seq(valid), y(valid), 1);
        slopes(s) = p(1);
    end
end
slopes = slopes(~isnan(slopes));
n = length(slopes);

group_slope = mean(slopes);
group_sd    = std(slopes);
group_se    = group_sd / sqrt(n);
t_obs       = group_slope / group_se;
p_two_sided = 2 * tcdf(-abs(t_obs), n-1);

fprintf('\n=== MC-3: Bayesian Equivalence Analysis ===\n');
fprintf('N participants with valid slope: %d\n', n);
fprintf('Mean slope: %.4f CoP-RMS/trial\n', group_slope);
fprintf('SD:         %.4f\n', group_sd);
fprintf('t(%d) = %.3f,  p (two-sided) = %.4f\n\n', n-1, t_obs, p_two_sided);

%% --- TOST Equivalence Test ---
% Equivalence bound: smallest meaningful slope (e.g., 0.1 SD = small effect)
% A slope smaller than this in absolute value is considered "equivalent to 0"
equiv_bound = 0.1 * group_sd; % Cohen's d = 0.1 (negligible effect)

% TOST: two one-sided t-tests
% H0a: slope <= -delta  (test H1a: slope > -delta)
% H0b: slope >= +delta  (test H1b: slope < +delta)
t_lower = (group_slope - (-equiv_bound)) / group_se;
t_upper = (group_slope -   equiv_bound)  / group_se;
p_lower = tcdf( t_lower, n-1, 'upper'); % p for H0a
p_upper = tcdf(-t_upper, n-1, 'upper'); % p for H0b ... wait
p_lower = 1 - tcdf(t_lower, n-1);
p_upper =     tcdf(t_upper, n-1);
p_tost  = max(p_lower, p_upper); % TOST p = larger of the two

fprintf('--- TOST Equivalence Test ---\n');
fprintf('Equivalence bound (delta): ±%.4f (0.1 SD, negligible effect)\n', equiv_bound);
fprintf('p_TOST = %.4f\n', p_tost);
if p_tost < 0.05
    fprintf('CONCLUSION: Slope is EQUIVALENT to zero (p_TOST < 0.05)\n');
    fprintf('→ The trajectory is statistically flat within the negligible-effect band.\n\n');
else
    fprintf('CONCLUSION: Cannot confirm equivalence to zero (p_TOST >= 0.05)\n');
    fprintf('→ Data insufficient to prove flatness; state as "consistent with stability."\n\n');
end

%% --- Bayes Factor (JZS prior, r=0.707) ---
% BF01: evidence for H0 (flat) over H1 (non-zero slope)
% Using Rouder et al. (2009) formula for one-sample t-test BF
r_cauchy = sqrt(2)/2; % JZS default scale = 0.707
integrand = @(g) (1 + n*g*r_cauchy^2).^(-0.5) .* ...
                 (1 + t_obs^2 ./ ((1 + n*g*r_cauchy^2)*(n-1))).^(-n/2) .* ...
                 (2*pi).^(-0.5) .* g.^(-3/2) .* exp(-1./(2*g));

BF10_num = integral(integrand, 0, Inf, 'RelTol', 1e-6, 'AbsTol', 1e-10);
BF10_den = (1 + t_obs^2/(n-1))^(-n/2);
BF10     = BF10_num / BF10_den;
BF01     = 1 / BF10;

fprintf('--- Bayes Factor (JZS prior, r=0.707) ---\n');
fprintf('BF10 = %.3f  (evidence for H1: non-zero slope)\n', BF10);
fprintf('BF01 = %.3f  (evidence for H0: slope = 0)\n', BF01);
if BF01 > 3
    fprintf('CONCLUSION: Moderate-to-strong evidence for H0 (flat slope).\n');
    fprintf('→ Supports: "the trajectory is consistent with a stable response."\n\n');
elseif BF01 > 1
    fprintf('CONCLUSION: Weak evidence for H0 — inconclusive.\n\n');
else
    fprintf('CONCLUSION: Evidence favours H1 (non-zero slope) — learning may be present.\n\n');
end

%% --- Figure: Slope distribution + equivalence bounds ---
fig = figure('Color','w','Position',[100 100 700 400]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Panel 1: Histogram of per-subject slopes
nexttile;
histogram(slopes, 8, 'FaceColor',[0.3 0.6 0.9], 'EdgeColor','w');
hold on;
% Keep explicit handles. Passing a bare label list here attached "Group mean"
% to the lower +/-delta line and left the actual group-mean line unlabelled,
% so the legend showed the group mean as a red dashed line while the plot drew
% it solid blue. Fixed 2026-07-25.
h_zero  = xline(0, 'k--', 'LineWidth', 1.5);
h_bound = xline( equiv_bound, 'r--', 'LineWidth', 1.2);
          xline(-equiv_bound, 'r--', 'LineWidth', 1.2);
h_mean  = xline(group_slope, 'b-', 'LineWidth', 2.5);
hold off;
xlabel('Linear slope (CoP RMS / trial)', 'FontSize', 11);
ylabel('Number of participants', 'FontSize', 11);
title(sprintf('Per-participant slopes\nMean = %.3f, p = %.3f', group_slope, p_two_sided), ...
    'FontSize', 10);
legend([h_zero h_bound h_mean], {'Zero', '±δ bound', 'Group mean'}, ...
       'Location','northwest','FontSize',8);
set(gca, 'FontSize', 10);

% Panel 2: Bayes Factor interpretation chart
nexttile;
bf_vals   = [0.01, 0.1, 1/3, 1, 3, 10, 100];
bf_labels = {'100:1 for H1','10:1 for H1','3:1 for H1','No evidence','3:1 for H0','10:1 for H0','100:1 for H0'};
bf_colors = [0.8 0.2 0.2; 0.9 0.5 0.2; 0.9 0.8 0.3; 0.7 0.7 0.7; 0.5 0.8 0.5; 0.2 0.7 0.3; 0.1 0.5 0.1];
for k = 1:length(bf_vals)
    barh(k, log10(bf_vals(k)), 'FaceColor', bf_colors(k,:), 'EdgeColor','w');
    hold on;
end
xline(log10(BF01), 'k-', 'LineWidth', 3);
text(log10(BF01)+0.05, 4, sprintf(' BF_{01}=%.2f', BF01), 'FontSize', 11, 'FontWeight','bold');
yticks(1:7); yticklabels(bf_labels); xlim([-2.5 2.5]);
xlabel('log_{10}(BF_{01})', 'FontSize', 11);
title('Bayes Factor (H0: flat slope)', 'FontSize', 10);
set(gca, 'FontSize', 9);
hold off;

% "MC-3" is an internal reviewer-response label and must not appear in a
% published figure.
title(tl, 'Bayesian equivalence analysis: temporal learning slope', 'FontSize', 12);
exportgraphics(fig, fullfile(paper_fig_dir, 'bayesian_equivalence_learning.png'), 'Resolution', 300);
fprintf('Figure saved to: %s\n', fullfile(paper_fig_dir, 'bayesian_equivalence_learning.png'));

%% --- Save text results to file (for Claude Code to read) ---
results_txt = fullfile(paper_fig_dir, 'bayesian_equivalence_results.txt');
fid = fopen(results_txt, 'w');
fprintf(fid, '=== MC-3: Bayesian Equivalence Analysis Results ===\n');
fprintf(fid, 'Date: %s\n\n', datestr(now));
fprintf(fid, 'N subjects: %d\n', n);
fprintf(fid, 'N trials per subject: %d\n', N_trials);
fprintf(fid, 'Mean slope: %.6f CoP-RMS/trial\n', group_slope);
fprintf(fid, 'SD: %.6f\n', group_sd);
fprintf(fid, 'SE: %.6f\n', group_se);
fprintf(fid, 't(%d) = %.4f, p (two-sided) = %.4f\n\n', n-1, t_obs, p_two_sided);
fprintf(fid, '--- TOST Equivalence Test ---\n');
fprintf(fid, 'Equivalence bound (delta): +/-%.6f (0.1 SD)\n', equiv_bound);
fprintf(fid, 'p_lower = %.4f\n', p_lower);
fprintf(fid, 'p_upper = %.4f\n', p_upper);
fprintf(fid, 'p_TOST = %.4f\n', p_tost);
fprintf(fid, 'TOST conclusion: %s\n\n', ...
    ternary(p_tost < 0.05, 'EQUIVALENT to zero', 'Cannot confirm equivalence'));
fprintf(fid, '--- Bayes Factor (JZS prior, r=0.707) ---\n');
fprintf(fid, 'BF10 = %.4f (evidence for non-zero slope)\n', BF10);
fprintf(fid, 'BF01 = %.4f (evidence for flat slope)\n', BF01);
if BF01 > 10
    bf_interp = 'Strong evidence for H0 (flat)';
elseif BF01 > 3
    bf_interp = 'Moderate evidence for H0 (flat)';
elseif BF01 > 1
    bf_interp = 'Weak/anecdotal evidence for H0';
else
    bf_interp = 'Evidence favours H1 (non-zero slope)';
end
fprintf(fid, 'Interpretation: %s\n', bf_interp);
fclose(fid);
fprintf('Results saved to: %s\n', results_txt);

%% --- Save .mat results ---
save(fullfile(out_dir, 'bayesian_equivalence_results.mat'), ...
    'slopes', 'group_slope', 'group_sd', 'group_se', ...
    't_obs', 'p_two_sided', 'equiv_bound', 'p_tost', ...
    'BF10', 'BF01', 'n', 'N_trials', '-v7.3');
fprintf('MAT results saved to: %s\n', fullfile(out_dir, 'bayesian_equivalence_results.mat'));

fprintf('=== DONE ===\n');

%% --- Helper ---
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
