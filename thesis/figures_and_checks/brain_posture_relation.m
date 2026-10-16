% Section 5.6: did participants with a larger Surprise ERD also show a larger late postural response?
% Surprise ERD = each participant's value at CPz (erd_per_participant_cpz.csv, from fig_erd_per_participant.m).
% Late postural response = each participant's mean over failure trials of the measure in Figure 5.5:
% foot data as saved by preprocess_eeg_and_foot.m step 4 (trials kept by foot_cleaning.m, 60 Hz low-pass, minus the mean
% of the first second of the trial), root mean square over 1-1.5 s after the trigger.
clear; clc;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % repository root
addpath(fullfile(root,'analysis','functions','preprocessing'));
config = configure_parameters();
E = readtable(fullfile(root,'thesis','figures_and_checks','erd_per_participant_cpz.csv'));   % written by fig_erd_per_participant.m
FF = load(fullfile(config.out_dir,'foot','foot_cleaned_session_combined_filtered.mat'),'foot_filtered');
FF = FF.foot_filtered;
tt = config.epoch_start : 1/config.eeg_fs : config.epoch_end; lw = tt >= 1.0 & tt <= 1.5;
late = nan(height(E),1);
for k = 1:height(E)
    X = FF.(sprintf('subject%d', E.participant(k))).epochs_unexpected;
    late(k) = mean(sqrt(mean(X(:,lw).^2, 2)));
end
[r, p, lo, hi] = corrcoef(E.surprise_erd_cpz_db, late);
fprintf('N = %d participants\n', height(E));
fprintf('Surprise ERD vs late postural response (failure trials): r = %.3f, p = %.3f, 95%% CI %.2f to %.2f\n', ...
    r(1,2), p(1,2), lo(1,2), hi(1,2));
