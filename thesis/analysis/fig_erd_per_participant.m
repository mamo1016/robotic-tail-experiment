% Per-participant Surprise ERD at CPz, for the thesis (Fig. erd_per_participant).
% The calculation is copied from codes_improved/stats_beta_erd_ttest_N25.m (STEP B, N = 25),
% and the script stops if the group values do not reproduce the published CSV.
clear; clc;
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'codes_improved', 'functions', 'for_obtainERD'));
config  = configure_parameters();
out_dir = config.out_dir;
fig_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');

beta_range = [13 30]; analysis_window = [0.0 1.0]; n_ch = 9;
ch_labels = {'FC1','FCz','FC2','C1','Cz','C2','CP1','CPz','CP2'};
cpz = find(strcmp(ch_labels,'CPz'));

tmp = load(fullfile(out_dir,'subject_number_list.mat'),'subject_number_list');
subject_number_list = tmp.subject_number_list;
unique_num = unique(subject_number_list);
st = load(fullfile(out_dir,'ga_feature','sub_tile.mat'),'sub_tile');
st_pos = sort(cellfun(@(x) sscanf(x,'sub_%d'), fieldnames(st.sub_tile)));
keep_ids = sort(reshape(unique_num(st_pos),1,[]));

files = dir(fullfile(out_dir,'ersp_ave','ersp_epoch_ave_sub_*.mat'));
fnums = cellfun(@(x) sscanf(x,'ersp_epoch_ave_sub_%d.mat'), {files.name});
[~,o] = sort(fnums); files = files(o);

sub_data = containers.Map('KeyType','char','ValueType','any');
for f = 1:min(numel(files), numel(subject_number_list))
    sid = num2str(subject_number_list(f));
    d = load(fullfile(files(f).folder, files(f).name));
    ss = fieldnames(d.ersp_epoch);
    for s = 1:numel(ss)
        e = d.ersp_epoch.(ss{s});
        if ~isfield(e,'unexpected') || ~isfield(e,'expected'), continue; end
        u = e.unexpected.data_ersp; x = e.expected.data_ersp;
        bi = u.freq >= beta_range(1) & u.freq <= beta_range(2);
        ti = u.time >= analysis_window(1) & u.time <= analysis_window(2);
        v = NaN(1,n_ch);
        for ch = 1:n_ch
            fn = sprintf('all_trial_ersp_db_single_ch_%d', ch);
            dm = squeeze(mean(u.(fn),1)) - squeeze(mean(x.(fn),1));
            v(ch) = mean(mean(dm(bi,ti)));
        end
        if isKey(sub_data,sid), sub_data(sid) = [sub_data(sid); v]; else, sub_data(sid) = v; end
    end
end
ids = cellfun(@str2double, keys(sub_data)); k = keys(sub_data);
[ids,o] = sort(ids); k = k(o);
M = NaN(numel(ids), n_ch);
for i = 1:numel(ids), M(i,:) = mean(sub_data(k{i}),1,'omitnan'); end
keep = ismember(ids(:), keep_ids(:));
vals = M(keep, cpz); pid = ids(keep);

% --- check against the published result (stats_beta_erd_ttest_results_N25.csv: CPz -0.9361, SD 2.7294)
m = mean(vals); sd = std(vals); n = numel(vals);
fprintf('N = %d, CPz mean = %.4f dB, SD = %.4f\n', n, m, sd);
assert(n == 25 && abs(m - (-0.9361)) < 5e-4 && abs(sd - 2.7294) < 5e-4, 'Does not reproduce the published values');

writetable(table(pid(:), vals(:), 'VariableNames', {'participant','surprise_erd_cpz_db'}), ...
    fullfile(fileparts(mfilename('fullpath')), 'erd_per_participant_cpz.csv'));

% --- figure: one dot per participant, sorted, with the mean and its 95% confidence interval
[vs, oi] = sort(vals); ci = tinv(0.975, n-1) * sd / sqrt(n);
fig = figure('Color','w','Units','centimeters','Position',[2 2 16 8]);
hold on;
fill([0.5 n+0.5 n+0.5 0.5], [m-ci m-ci m+ci m+ci], [0.85 0.9 1], 'EdgeColor','none');
yline(0, 'k-', 'LineWidth', 0.8);
yline(m, '-', 'Color', [0.1 0.3 0.8], 'LineWidth', 1.5);
scatter(1:n, vs, 36, [0.1 0.3 0.8], 'filled', 'MarkerEdgeColor', 'w');
xlim([0.5 n+0.5]); box off;
xticks(1:n); xticklabels(string(pid(oi)));
xlabel('Participant (sorted by Surprise ERD)');
ylabel('Surprise ERD at CPz (dB)');
text(n-0.2, m+ci+0.35, sprintf('mean %.2f dB, 95%% CI %.2f to %.2f', m, m-ci, m+ci), ...
    'HorizontalAlignment','right','Color',[0.1 0.3 0.8],'FontSize',8);
set(gca,'FontSize',9,'TickDir','out');
exportgraphics(fig, fullfile(fig_dir,'erd_per_participant_cpz.png'), 'Resolution', 300);
fprintf('Saved figure. Participants below -1 dB: %d; within +/-0.5 dB: %d\n', sum(vals < -1), sum(abs(vals) <= 0.5));
fprintf('min %.2f, max %.2f\n', min(vals), max(vals));
