% Plot the per-participant Surprise ERD at CPz from erd_per_participant_cpz.csv
% (the CSV is made by fig_erd_per_participant.m, which checks it against the published values).
here = fileparts(mfilename('fullpath'));
T = readtable(fullfile(here, 'erd_per_participant_cpz.csv'));
vals = T.surprise_erd_cpz_db; pid = T.participant; n = numel(vals);
m = mean(vals); sd = std(vals); ci = tinv(0.975, n-1) * sd / sqrt(n);
[vs, oi] = sort(vals);
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
set(gca,'FontSize',9,'TickDir','out');
exportgraphics(fig, fullfile(fileparts(mfilename('fullpath')), 'figures', 'erd_per_participant_cpz.png'), 'Resolution', 300);
fprintf('mean %.2f, 95%% CI %.2f to %.2f, median %.2f, negative %d of %d\n', m, m-ci, m+ci, median(vals), sum(vals<0), n);
