% One-participant walkthrough figures for the thesis (Figures 5.1 and 5.2): every processing stage,
% from the EMG trigger to the brain and foot measures.
% Example participant: the one whose Surprise ERD at CPz is closest to the group mean (not the
% clearest one). Uses the analysis pipeline's own functions and saved outputs.
clear; clc;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % repository root
addpath(fullfile(root,'analysis','functions','preprocessing'));
addpath(fullfile(root,'analysis'));
config = configure_parameters();
fig_dir = fullfile(root,'thesis','revised','images','results');
here = fullfile(root,'thesis','figures_and_checks');

% ---------- choose the example participant
T = readtable(fullfile(here,'erd_per_participant_cpz.csv'));
[~, i] = min(abs(T.surprise_erd_cpz_db - mean(T.surprise_erd_cpz_db)));
P = T.participant(i); S = 1;
fprintf('Example participant %d (Surprise ERD %.2f dB, group mean %.2f), session %d\n', ...
    P, T.surprise_erd_cpz_db(i), mean(T.surprise_erd_cpz_db), S);

% ---------- saved epochs and single-trial ERSP for this session
E = load(fullfile(config.out_dir,'epochs',sprintf('subject_%d___session%d___epoch_data.mat',P,S)));
E = E.epoch_data;
fs = config.eeg_fs; t_ep = (0:size(E.epochs_expected.eeg,2)-1)/fs + config.epoch_start;   % -2..+4 s
cpz = 8;                                  % 9-channel order FC1 FCz FC2 C1 Cz C2 CP1 CPz CP2
coi = config.coi;                         % raw-file column of each analysed channel

sl = load(fullfile(config.out_dir,'subject_number_list.mat')); sl = sl.subject_number_list(:);
k = find(sl == P); k = k(S);                     % ersp_epoch_ave_sub_<k>.mat is this session
d = load(fullfile(config.out_dir,'ersp_ave',sprintf('ersp_epoch_ave_sub_%d.mat',k)));
key = sprintf('subject%02d_session%02d', P, S); R = [];
if isfield(d.ersp_epoch, key), R = d.ersp_epoch.(key); end
assert(~isempty(R), 'ERSP for this session not found');
Rx = R.expected.data_ersp; Ru = R.unexpected.data_ersp;
fchan = sprintf('all_trial_ersp_db_single_ch_%d', cpz);
assert(size(Rx.(fchan),1) == size(E.epochs_expected.eeg,1), 'Trial counts differ between epochs and ERSP');

% ---------- raw session data through the pipeline's own loader (for raw EEG and raw pressure cells)
subj = fullfile(config.baseDir, sprintf('Sub%02d', P));
files = all_files(subj);
[mf, td] = eeg_emg_timeoffset(files, S, config);
D = load_data_optimized(subj, mf, config);
D = footdist(D); D = footCompute_optimized(D); D = foot_load_average(D, config);
offset = seconds(td);

% one expected trial from the middle of the session, and one failure trial
tr = round(size(E.idxs_expected,1)/2);
ix = E.idxs_expected(tr,1):E.idxs_expected(tr,2);
iu = E.idxs_unexpected(ceil(end/2),1):E.idxs_unexpected(ceil(end/2),2);
t0 = E.time_master(ix(1)) - config.epoch_start;      % trigger time on the master (EMG) clock

raw_cpz = double(D.eeg_raw(ix, cpz));          % loader keeps the 9 analysed channels
flt_cpz = double(E.data_eeg(ix, cpz));
fp = design_filter(config); Fr = filter_eeg_vectorized(D.eeg_raw, fp, config);
fprintf('check: pipeline filter on raw vs saved filtered EEG at CPz, r = %.4f\n', corr(double(Fr(ix,cpz)), flt_cpz));

% ================= EEG figure
% top row: one trial at each stage (a-c); bottom row: normal | failure | difference side by side (d-f)
fig = figure('Color','w','Units','centimeters','Position',[1 1 24 13], ...
    'DefaultAxesFontSize',11);                          % shrunk to the text width in print; 11 pt keeps labels readable
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
lim = [-2 4]; fl = [4 40];

nexttile; hold on;
plot(t_ep, E.epochs_expected.emg(tr,:), 'Color',[0 0.45 0.74]);
plot(t_ep, E.epochs_unexpected.emg(ceil(end/2),:), 'Color',[0.85 0.2 0.1]);
xline(0,'k--'); xlim(lim); box off;
xlabel('Time from trigger (s)'); ylabel('Biceps EMG (raw units)');
lg = legend({'normal trial','failure trial'},'Box','off','Location','northeast');   % top right is empty
lg.ItemTokenSize = [12 18];                                                          % short sample lines keep it right of the peak at t = 0
title('(a) EMG trigger at t = 0','FontWeight','normal');

nexttile; hold on;
rr = raw_cpz - mean(raw_cpz);
plot(t_ep(1:numel(rr)), rr, 'Color',[0.6 0.6 0.6]);
plot(t_ep(1:numel(flt_cpz)), flt_cpz - 15, 'Color',[0 0 0]);
xline(0,'k--'); xlim(lim); box off;
xlabel('Time from trigger (s)'); ylabel('EEG at CPz (\muV)');
lg = legend({'raw','filtered, shifted down'},'Location','east', ...   % in the empty band between the two traces;
    'Color','w','EdgeColor','w');                                                     % white so the t = 0 line does not cross the text
lg.ItemTokenSize = [12 18];                                                           % short sample lines keep the box clear of the y-axis
title('(b) Raw and filtered EEG','FontWeight','normal');

tf = Rx.time; ff = Rx.freq(:); fi = ff >= fl(1) & ff <= fl(2);
nexttile;
imagesc(tf, ff(fi), squeeze(Rx.(fchan)(tr,fi,:))); axis xy; xline(0,'k--');
xlim(lim); caxis([-10 10]); colormap(gca, turbo); cb = colorbar; cb.Label.String = 'dB';
xlabel('Time from trigger (s)'); ylabel('Frequency (Hz)');
title('(c) Power of this trial','FontWeight','normal');

mx = squeeze(mean(Rx.(fchan)(:,fi,:),1)); mu = squeeze(mean(Ru.(fchan)(:,fi,:),1));
cl = max(abs([mx(:); mu(:)])); cl = min(cl, 4);
nexttile; imagesc(tf, ff(fi), mx); axis xy; xline(0,'k--'); xlim(lim); caxis([-cl cl]);
colormap(gca, redblue()); cb = colorbar; cb.Label.String = 'dB';
xlabel('Time from trigger (s)'); ylabel('Frequency (Hz)');
title(sprintf('(d) Average of %d normal trials', size(Rx.(fchan),1)),'FontWeight','normal');   % same colour scale in d-f

nexttile; imagesc(tf, ff(fi), mu); axis xy; xline(0,'k--'); xlim(lim); caxis([-cl cl]);
colormap(gca, redblue()); cb = colorbar; cb.Label.String = 'dB';
xlabel('Time from trigger (s)'); ylabel('Frequency (Hz)');
title(sprintf('(e) Average of %d failure trials', size(Ru.(fchan),1)),'FontWeight','normal');

dm = mu - mx; bi = ff(fi) >= 13 & ff(fi) <= 30; ti = tf >= 0 & tf <= 1;
val = mean(mean(dm(bi,ti)));
nexttile; imagesc(tf, ff(fi), dm); axis xy; hold on; xline(0,'k--');
rectangle('Position',[0 13 1 17],'EdgeColor','k','LineWidth',1.2);
xlim(lim); caxis([-cl cl]); colormap(gca, redblue()); cb = colorbar; cb.Label.String = 'dB';
xlabel('Time from trigger (s)'); ylabel('Frequency (Hz)');
title(sprintf('(f) Failure - normal (box: %.2f dB)', val),'FontWeight','normal');
exportgraphics(fig, fullfile(fig_dir,'walkthrough_eeg.png'), 'Resolution', 250);
fprintf('Session %d Surprise ERD (box 13-30 Hz, 0-1 s): %.3f dB\n', S, val);

% ================= Foot figure
fig2 = figure('Color','w','Units','centimeters','Position',[1 1 18 14]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
tfoot = D.foot_time_dist_Both(:);
w = tfoot >= t0 + lim(1) & tfoot <= t0 + lim(2);
nexttile; hold on;
cl5 = lines(5);
for c = 1:5
    plot(tfoot(w)-t0, D.foot_raw_dist_L(w,c), '-',  'Color', cl5(c,:));
    plot(tfoot(w)-t0, D.foot_raw_dist_R(w,c), '--', 'Color', cl5(c,:));
end
xline(0,'k--'); xlim(lim); box off;
xlabel('Time from trigger (s)'); ylabel('Pressure (sensor units)');
title('(a) Raw pressure, 10 cells','FontWeight','normal');

nexttile;
plot(t_ep, E.epochs_expected.foot(tr,:), 'k'); xline(0,'k--'); xlim(lim); box off;
xlabel('Time from trigger (s)');
ylabel({'Side-to-side position (a.u.)', 'left \leftarrow    \rightarrow right'});   % left-foot cells weighted -6..-1, right +1..+6 (configure_parameters.m)
title('(b) Side-to-side position','FontWeight','normal');

bw = t_ep >= -0.5 & t_ep <= 0;
Fx = E.epochs_expected.foot - mean(E.epochs_expected.foot(:,bw),2);
Fu = E.epochs_unexpected.foot - mean(E.epochs_unexpected.foot(:,bw),2);
nexttile; hold on;
plot(t_ep, Fx', 'Color',[0 0.45 0.74 0.15]); plot(t_ep, Fu', 'Color',[0.85 0.2 0.1 0.3]);
h1 = plot(t_ep, mean(Fx,1), 'Color',[0 0.3 0.6], 'LineWidth',2);
h2 = plot(t_ep, mean(Fu,1), 'Color',[0.7 0.1 0.05], 'LineWidth',2);
xline(0,'k--'); xlim(lim); box off;
xlabel('Time from trigger (s)'); ylabel('Change in position (a.u.)');
legend([h1 h2], {'normal (mean)','failure (mean)'}, 'Box','off','Location','best');
title('(c) All trials, baseline removed','FontWeight','normal');

% (d) the late response exactly as in the session analysis (Figure 5.5; preprocess_eeg_and_foot.m step 4
% and posture_across_session.m): trials kept by foot_cleaning.m's rule, 60 Hz low-pass, minus the
% mean of the first second of the trial, root mean square over 1-1.5 s. (The trial table's
% cop_LateComp_rms column is NOT this: it has no baseline removed, so it mostly measures where the
% person stands.) Trials of both types are put in time order within each session.
tt = config.epoch_start : 1/config.eeg_fs : config.epoch_end; lw = tt >= 1.0 & tt <= 1.5;
late = @(x) sqrt(mean(x(lw).^2));
rows = [];                                        % [session, onset sample, type (1 normal, 2 failure), late RMS]
chk = {[], []};                                   % values per type, in the order preprocess_eeg_and_foot.m combines them
types = {'epochs_expected','epochs_unexpected'}; idxf = {'idxs_expected','idxs_unexpected'};
for s = 1:3
    Es = load(fullfile(config.epochDir, sprintf('subject_%d___session%d___epoch_data.mat',P,s))); Es = Es.epoch_data;
    for ty = 1:2
        X = Es.(types{ty}).foot; on = Es.(idxf{ty})(:,1);
        for i = 1:size(X,1)
            if ~keep_foot_trial(X(i,:)), continue; end
            y = lowpass(X(i,:), config.foot_fpass(2), config.eeg_fs);
            y = y - mean(y(1:1200));
            v = late(y);
            rows(end+1,:) = [s, on(i), ty, v]; %#ok<AGROW>
            chk{ty}(end+1,1) = v;
        end
    end
end
FF = load(fullfile(config.out_dir,'foot','foot_cleaned_session_combined_filtered.mat'), 'foot_filtered');
ref = FF.foot_filtered.(sprintf('subject%d',P));
for ty = 1:2
    r = arrayfun(@(i) late(ref.(types{ty})(i,:)), (1:size(ref.(types{ty}),1))');
    assert(isequal(size(r), size(chk{ty})) && max(abs(r - chk{ty})) < 1e-9, ...
        'Panel (d) values differ from the data behind Figure 5.5');
end
fprintf('check: panel (d) matches the session-analysis data (%d normal, %d failure trials)\n', numel(chk{1}), numel(chk{2}));
rows = sortrows(rows, [1 2]); ord = (1:size(rows,1))';
nexttile; hold on;
a = rows(:,3) == 1; b = rows(:,3) == 2;
scatter(ord(a), rows(a,4), 10, [0 0.45 0.74], 'filled', 'MarkerFaceAlpha', 0.5);
scatter(ord(b), rows(b,4), 22, [0.85 0.2 0.1], 'filled');
sb = find(diff(rows(:,1)) ~= 0) + 0.5; for k = sb', xline(k, ':', 'Color', [0.5 0.5 0.5]); end
box off; xlabel('Trial number (all three sessions)'); ylabel('Late response, 1-1.5 s (RMS, a.u.)');
legend({'normal','failure'},'Box','off','Location','best');
title('(d) Late response on every trial','FontWeight','normal'); xlim([0 numel(ord)+1]);
exportgraphics(fig2, fullfile(fig_dir,'walkthrough_foot.png'), 'Resolution', 250);
fprintf('Done. Participant %d: %d normal and %d failure trials in session %d.\n', P, ...
    size(E.epochs_expected.eeg,1), size(E.epochs_unexpected.eeg,1), S);

function ok = keep_foot_trial(x)
% The trial-keeping rule of analysis/functions/foot_pressure/foot_cleaning.m: drop a trial whose foot signal has
% a stretch of 1000 or more samples with no change in slope (a frozen or extrapolated sensor).
ddf = diff(diff(x)); run = 0; mx = 0;
for k = 1:numel(ddf)
    if abs(ddf(k)) <= 0.001, run = run + 1; mx = max(mx, run); else, run = 0; end
end
rs = [1, find(diff(ddf) ~= 0) + 1]; rl = diff([rs, numel(ddf) + 1]);
ok = max(rl) < 1000 && mx < 1000;
end

function c = redblue()
    n = 128; r = [linspace(0,1,n) ones(1,n)]'; g = [linspace(0,1,n) linspace(1,0,n)]'; b = [ones(1,n) linspace(1,0,n)]';
    c = [r g b];
end
