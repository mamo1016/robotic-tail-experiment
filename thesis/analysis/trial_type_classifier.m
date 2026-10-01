% Can foot pressure predict the trial type (tail moved vs failure)?  Examiner request (Ch 6).
% Method fixed before running:
%   trials   : expected (cond 1) and unexpected (cond 2), participants 1-25, frac_extrap == 0
%   features : foot side-to-side signal, baseline -0.5..0 s removed, mean in 100 ms bins from 0 to 3 s
%              (30 features), z-scored within each participant (no labels used)
%   model    : logistic regression, classes weighted equally, leave-one-participant-out
%   score    : AUC per held-out participant (0.5 = chance) and balanced accuracy
%   chance   : 200 runs with labels shuffled within each participant
%   control  : the same with one feature only, overall response size (RMS 0-4 s)
clear; clc; rng(1);
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % repository root
addpath(fullfile(root,'codes_improved','functions','for_obtainERD'));
config = configure_parameters();
here = fullfile(root,'thesis','analysis');
fig_dir = fullfile(root,'thesis','revised','images','results');

V = readtable(fullfile(config.out_dir,'trial_level','cop_corrected.csv'));
fs = config.eeg_fs; t = (0:7200)/fs + config.epoch_start;
bw = t >= -0.5 & t <= 0; edges = 0:0.1:3; nb = numel(edges)-1;
X = []; Xs = []; y = []; g = [];
for P = 1:25
    for S = 1:3
        f = fullfile(config.out_dir,'epochs',sprintf('subject_%d___session%d___epoch_data.mat',P,S));
        if ~exist(f,'file'), continue; end
        E = load(f); E = E.epoch_data;
        for c = 1:2
            if c == 1, F = E.epochs_expected.foot; else, F = E.epochs_unexpected.foot; end
            v = V(V.subject==P & V.session==S & V.cond==c, :);
            assert(height(v) == size(F,1), 'trial count mismatch P%d S%d c%d', P, S, c);
            ok = v.frac_extrap == 0 & all(isfinite(F),2);
            F = F(ok,:) - mean(F(ok,bw),2);
            B = zeros(size(F,1), nb);
            for k = 1:nb, B(:,k) = mean(F(:, t >= edges(k) & t < edges(k+1)), 2); end
            X = [X; B]; Xs = [Xs; sqrt(mean(F(:, t >= 0 & t <= 4).^2, 2))];
            y = [y; repmat(c == 2, size(F,1), 1)]; g = [g; repmat(P, size(F,1), 1)];
        end
    end
end
y = logical(y);
fprintf('Trials: %d (%d failure) from %d participants\n', numel(y), sum(y), numel(unique(g)));

% z-score within participant (unsupervised)
for P = unique(g)'
    m = g == P;
    X(m,:) = (X(m,:) - mean(X(m,:))) ./ std(X(m,:));
    Xs(m)  = (Xs(m) - mean(Xs(m))) ./ std(Xs(m));
end

[auc, bacc, sc] = lopo(X, y, g);
[auc_s, bacc_s] = lopo(Xs, y, g);
fprintf('\nTime course: mean AUC %.3f (median %.3f), balanced accuracy %.3f\n', mean(auc), median(auc), mean(bacc));
fprintf('Size only  : mean AUC %.3f, balanced accuracy %.3f\n', mean(auc_s), mean(bacc_s));
fprintf('Participants with AUC > 0.5: %d of %d\n', sum(auc > 0.5), numel(auc));

nperm = 200; null = zeros(nperm,1);
for r = 1:nperm
    yp = y;
    for P = unique(g)', m = find(g == P); yp(m) = y(m(randperm(numel(m)))); end
    null(r) = mean(lopo(X, yp, g));
end
p = (1 + sum(null >= mean(auc))) / (1 + nperm);
fprintf('Shuffled labels: mean AUC %.3f (95th percentile %.3f); p = %.4f\n', mean(null), prctile(null,95), p);

% which time bins carry the information: difference of the within-participant z-scored means
dif = mean(X(y==1,:)) - mean(X(y==0,:));

save(fullfile(here,'trial_type_classifier_results.mat'), 'auc','bacc','auc_s','bacc_s','null','p','dif','edges');

fig = figure('Color','w','Units','centimeters','Position',[2 2 17 7]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
[as, o] = sort(auc);
scatter(1:numel(as), as, 30, [0.1 0.3 0.8], 'filled', 'MarkerEdgeColor','w');
yline(0.5, 'k--'); yline(mean(auc), '-', 'Color', [0.1 0.3 0.8], 'LineWidth', 1.2);
ylim([0 1]); xlim([0.5 numel(as)+0.5]); box off;
xlabel('Held-out participant (sorted)'); ylabel('AUC (0.5 = chance)');
title('(a) Prediction for each participant','FontWeight','normal');
nexttile; hold on;
ctr = edges(1:end-1) + 0.05;
bar(ctr, dif, 1, 'FaceColor', [0.85 0.2 0.1], 'EdgeColor', 'none');
yline(0,'k-'); box off; xlim([0 3]);
xlabel('Time from trigger (s)'); ylabel('Failure minus normal (z)');
title('(b) Where the trial types differ','FontWeight','normal');
exportgraphics(fig, fullfile(fig_dir,'trial_type_classifier.png'), 'Resolution', 250);

function [auc, bacc, sc] = lopo(X, y, g)
    ps = unique(g)'; auc = zeros(numel(ps),1); bacc = auc; sc = zeros(size(y));
    for i = 1:numel(ps)
        te = g == ps(i); tr = ~te;
        w = ones(sum(tr),1); yt = y(tr);
        w(yt) = 0.5 / mean(yt); w(~yt) = 0.5 / mean(~yt);
        mdl = fitclinear(X(tr,:), yt, 'Learner','logistic', 'Weights', w, 'Regularization','ridge', 'Lambda', 1e-3);
        [~, s] = predict(mdl, X(te,:)); s = s(:,2); sc(te) = s;
        [~,~,~,auc(i)] = perfcurve(y(te), s, true);
        pr = s > 0.5; % scores are probabilities; classes were weighted equally
        bacc(i) = 0.5 * (mean(pr(y(te))) + mean(~pr(~y(te))));
    end
end
