% One-participant walkthrough figures for the thesis (examiner: "graphs showing how each stage of
% the analysis affects the data"; "examples ... from the EMG signals through to the ERD signal").
% Example participant: the one whose Surprise ERD at CPz is closest to the group mean (not the
% clearest one). Uses the analysis pipeline's own functions and saved outputs.
clear; clc;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % repository root
addpath(fullfile(root,'codes_improved','functions','for_obtainERD'));
addpath(fullfile(root,'codes_improved'));
config = configure_parameters();
fig_dir = fullfile(root,'thesis','revised','images','results');
here = fullfile(root,'thesis','analysis');

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
fig = figure('Color','w','Units','centimeters','Position',[1 1 18 20]);
tl = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
lim = [-2 4]; fl = [4 40];

nexttile; hold on;
plot(t_ep, E.epochs_expected.emg(tr,:), 'Color',[0 0.45 0.74]);
plot(t_ep, E.epochs_unexpected.emg(ceil(end/2),:), 'Color',[0.85 0.2 0.1]);
xline(0,'k--'); xlim(lim); box off;
xlabel('Time from trigger (s)'); ylabel('Biceps EMG (raw units)');
legend({'normal trial','failure trial'},'Box','off','Location','northeast');
title('(a) EMG trigger at t = 0','FontWeight','normal');

nexttile; hold on;
rr = raw_cpz - mean(raw_cpz);
plot(t_ep(1:numel(rr)), rr, 'Color',[0.6 0.6 0.6]);
plot(t_ep(1:numel(flt_cpz)), flt_cpz - 15, 'Color',[0 0 0]);
xline(0,'k--'); xlim(lim); box off;
xlabel('Time from trigger (s)'); ylabel('EEG at CPz (\muV)');
legend({'raw (mean removed)','filtered 1-40 Hz (shifted down)'},'Box','off','Location','northeast');
title('(b) Raw and filtered EEG','FontWeight','normal');

tf = Rx.time; ff = Rx.freq(:); fi = ff >= fl(1) & ff <= fl(2);
nexttile;
imagesc(tf, ff(fi), squeeze(Rx.(fchan)(tr,fi,:))); axis xy; xline(0,'k--');
xlim(lim); caxis([-10 10]); colormap(gca, turbo); cb = colorbar; cb.Label.String = 'dB';
xlabel('Time from trigger (s)'); ylabel('Frequency (Hz)');
title('(c) Time-frequency power','FontWeight','normal');

mx = squeeze(mean(Rx.(fchan)(:,fi,:),1)); mu = squeeze(mean(Ru.(fchan)(:,fi,:),1));
cl = max(abs([mx(:); mu(:)])); cl = min(cl, 4);
nexttile; imagesc(tf, ff(fi), mx); axis xy; xline(0,'k--'); xlim(lim); caxis([-cl cl]);
colormap(gca, redblue()); cb = colorbar; cb.Label.String = 'dB';
xlabel('Time from trigger (s)'); ylabel('Frequency (Hz)');
title(sprintf('(d) Average of %d normal trials', size(Rx.(fchan),1)),'FontWeight','normal');

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
title(sprintf('(f) Difference: %.2f dB in box', val),'FontWeight','normal');
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
xlabel('Time from trigger (s)'); ylabel('Side-to-side position (a.u.)');
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

TT = readtable(fullfile(config.out_dir,'trial_level','trial_table.csv'));
q = TT(TT.subject == P & ismember(TT.cond,[1 2]), :);
q = sortrows(q, {'session','onset'});
ord = (1:height(q))';
nexttile; hold on;
a = q.cond == 1; b = q.cond == 2;
scatter(ord(a), q.cop_LateComp_rms(a), 10, [0 0.45 0.74], 'filled', 'MarkerFaceAlpha', 0.5);
scatter(ord(b), q.cop_LateComp_rms(b), 22, [0.85 0.2 0.1], 'filled');
sb = find(diff(q.session) ~= 0) + 0.5; for k = sb', xline(k, ':', 'Color', [0.5 0.5 0.5]); end
box off; xlabel('Trial number (all three sessions)'); ylabel('Late response, 1-1.5 s (RMS, a.u.)');
legend({'normal','failure'},'Box','off','Location','best');
title('(d) Late response on every trial','FontWeight','normal'); xlim([0 height(q)+1]);
exportgraphics(fig2, fullfile(fig_dir,'walkthrough_foot.png'), 'Resolution', 250);
fprintf('Done. Participant %d: %d normal and %d failure trials in session %d.\n', P, ...
    size(E.epochs_expected.eeg,1), size(E.epochs_unexpected.eeg,1), S);

function c = redblue()
    n = 128; r = [linspace(0,1,n) ones(1,n)]'; g = [linspace(0,1,n) linspace(1,0,n)]'; b = [ones(1,n) linspace(1,0,n)]';
    c = [r g b];
end
