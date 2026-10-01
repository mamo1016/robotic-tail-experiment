function LOOCV_Plot(save_filename)
% LOOCV_Plot  Loads a saved LOOCV results file and generates all figures.
%
% Step 6b — fully self-contained, no workspace variables needed.
%
% Usage:
%   LOOCV_Plot([])          % auto-loads most recent results file
%   LOOCV_Plot(full_path)   % loads a specific results file

out_dir = 'output/ga_feature/';

if isempty(save_filename)
    files = dir(fullfile(out_dir, 'results_LOOCV_*.mat'));
    if isempty(files)
        error('No results file found in %s', out_dir);
    end
    [~, idx]  = max([files.datenum]);
    save_filename = fullfile(out_dir, files(idx).name);
end

fprintf('=== STEP 6b: LOOCV Plotting ===\n');
fprintf('Loading: %s\n', save_filename);
load(save_filename, 'results', 'run_info', 'eeg_data', 'foot_data', 'foot_labels');

[~, fname, ~] = fileparts(save_filename);
timestamp = strrep(fname, 'results_LOOCV_', '');

methods = fieldnames(results);
tile_t  = run_info.tile_t;
tile_f  = run_info.tile_f;

eeg_mean_diff  = mean(eeg_data, 1);
[h_eeg_all, ~] = ttest(eeg_data, 0, 'Alpha', 0.05);

phases        = {'Expected vs Expected', 'Unexpected vs Expected', 'Right After vs Right Before'};
phase_indices = {1:400, 401:800, 801:1200};

% -----------------------------------------------------------------------
% Fig A: Per-method ERSP maps
% -----------------------------------------------------------------------
for m = 1:length(methods)
    current_method = methods{m};
    mask_x = results.(current_method).mask_x;
    mask_y = results.(current_method).mask_y;
    mask_x_counts_method = results.(current_method).mask_x_counts;

    fig_vis      = figure('Color','w','Position',[100 50 600 1200],'Name',['EEG: ' current_method]);
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    fig_vis_only = figure('Color','w','Position',[750 50 600 1200],'Name',['EEG Sig: ' current_method]);
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    fig_counts   = figure('Color','w','Position',[750 50 600 1200],'Name',['EEG Counts: ' current_method]);
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    for p_idx = 1:3
        phase_name  = phases{p_idx};
        indices     = phase_indices{p_idx};
        vis_mask_x  = mask_x(indices);
        eeg_map_2d  = reshape(eeg_mean_diff(indices), [length(tile_f) length(tile_t)]);
        sig_mask_2d = reshape(h_eeg_all(indices),     [length(tile_f) length(tile_t)]);
        eeg_sig     = eeg_map_2d;
        eeg_sig(sig_mask_2d == 0) = 0;
        [sel_t, sel_f] = mask_visualiser(vis_mask_x, tile_f, tile_t);

        figure(fig_vis);
        visualisation_ersp(eeg_map_2d, fig_vis, tile_t, tile_f, [-0.5 0.5]);
        ylim([5 35]); hold on;
        plot(sel_t, sel_f, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8, 'LineWidth', 1.5);
        title(['EEG: ' phase_name], 'FontSize', 12);
        if p_idx==1, subtitle(['Method: ' current_method], 'FontSize', 10); end
        hold off;

        figure(fig_vis_only);
        visualisation_ersp(eeg_sig, fig_vis_only, tile_t, tile_f, [-0.5 0.5]);
        ylim([5 35]); hold on;
        plot(sel_t, sel_f, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8, 'LineWidth', 1.5);
        title(['EEG (Sig Only): ' phase_name], 'FontSize', 12);
        if p_idx==1, subtitle(['Method: ' current_method], 'FontSize', 10); end
        hold off;

        % Count-annotated figure
        p_counts   = mask_x_counts_method(indices);
        figure(fig_counts);
        visualisation_ersp(eeg_sig, fig_counts, tile_t, tile_f, [-0.5 0.5]);
        ylim([5 35]); hold on;
        plot(sel_t, sel_f, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 10, 'LineWidth', 1.5);
        p_tile_idx = find(vis_mask_x == 1);
        for kk = 1:length(p_tile_idx)
            tidx    = p_tile_idx(kk);
            cnt_val = p_counts(tidx);
            [fi, ti] = ind2sub([length(tile_f), length(tile_t)], tidx);
            text(tile_t(ti)+0.05, tile_f(fi)+0.8, sprintf('%d', cnt_val), ...
                'Color','black','FontSize',9,'FontWeight','bold', ...
                'HorizontalAlignment','left','BackgroundColor',[1 1 0.6]);
        end
        title(['EEG (Sig+Counts): ' phase_name], 'FontSize', 12);
        if p_idx==1, subtitle(['Method: ' current_method ' | Number = folds selected'], 'FontSize',9); end
        hold off;
    end

    exportgraphics(fig_vis, sprintf('%sselected_eeg_%s_%s.png', out_dir, current_method, timestamp), 'Resolution', 300);
    close(fig_vis);
    exportgraphics(fig_vis_only, sprintf('%sselected_eeg_ONLY_%s_%s.png', out_dir, current_method, timestamp), 'Resolution', 300);
    close(fig_vis_only);
    exportgraphics(fig_counts, sprintf('%sselected_eeg_COUNTS_%s_%s.png', out_dir, current_method, timestamp), 'Resolution', 300);
    close(fig_counts);
    fprintf('[%s] figures saved.\n', upper(current_method));

    fprintf('--- [%s] Consensus FOOT Features ---\n', upper(current_method));
    sel_foot = find(mask_y == 1);
    for k = 1:length(sel_foot)
        if sel_foot(k) <= length(foot_labels)
            fprintf('  %d. %s\n', k, foot_labels{sel_foot(k)});
        end
    end
    fprintf('\n');
end

% -----------------------------------------------------------------------
% Fig B: Clean ERSP (no dots)
% -----------------------------------------------------------------------
fprintf('Generating clean ERSP maps...\n');
fig_clean = figure('Color','w','Position',[100 50 600 1200]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
clean_titles = {'EEG: Expected vs Expected','EEG: Unexpected vs Expected','EEG: Right After vs Right Before'};
for p_idx = 1:3
    eeg_map_2d = reshape(eeg_mean_diff(phase_indices{p_idx}), [length(tile_f) length(tile_t)]);
    nexttile;
    imagesc(tile_t, tile_f, eeg_map_2d, [-0.5 0.5]); axis xy;
    cmap = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];
    colormap(cmap); colorbar;
    xlabel('Time (s)','FontSize',14); ylabel('Frequency (Hz)','FontSize',14);
    ylim([5 35]); hold on; xline(0,'k--','LineWidth',0.8); hold off;
    title(clean_titles{p_idx},'FontSize',16,'FontWeight','bold');
    set(gca,'FontSize',12);
end
exportgraphics(fig_clean, sprintf('%sersp_clean_%s.png', out_dir, timestamp), 'Resolution', 300);
close(fig_clean);

% -----------------------------------------------------------------------
% Fig C: Figure 4 — Foot bar + 3 scatter plots
% -----------------------------------------------------------------------
fprintf('Generating Figure 4...\n');
show_methods = {'spearman','ridge','svr'};
show_colors  = {[0.2 0.4 0.8],[0.3 0.7 0.9],[0.8 0.1 0.5]};
fig4 = figure('Color','w','Position',[50 50 1600 800]);
tiledlayout(3,4,'TileSpacing','compact','Padding','compact');

nexttile([3 2]);
n_foot   = length(foot_labels);
foot_vals = mean(foot_data, 1);
foot_mask = zeros(1, n_foot);
if isfield(results,'svr'), foot_mask = results.svr.mask_y; end
bar_colors = repmat([0.75 0.75 0.75], n_foot, 1);
sel_idx = find(foot_mask == 1);
for k = 1:length(sel_idx), bar_colors(sel_idx(k),:) = [0.9 0.7 0.1]; end
b = bar(foot_vals,'FaceColor','flat'); b.CData = bar_colors;
hold on;
for k = 1:length(sel_idx)
    plot(sel_idx(k), foot_vals(sel_idx(k)), 'p', 'MarkerSize', 12, 'MarkerFaceColor', [0.9 0.2 0.1], 'MarkerEdgeColor', 'k');
end
h1 = bar(NaN,'FaceColor',[0.75 0.75 0.75]);
h2 = bar(NaN,'FaceColor',[0.9 0.7 0.1]);
h3 = plot(NaN,NaN,'p','MarkerSize',12,'MarkerFaceColor',[0.9 0.2 0.1],'MarkerEdgeColor','k');
legend([h1 h2 h3],{'Not Selected','GA Selected','Selected Marker'},'Location','northwest','FontSize',9);
hold off;
ylabel('Feature Value (normalized)','FontSize',12); xlabel('Foot Feature Index','FontSize',12);
title('Foot CoP Features','FontSize',14,'FontWeight','bold');

for sm = 1:3
    nexttile([1 2]);
    method_name = show_methods{sm};
    p_scores = results.(method_name).predicted;
    a_scores = results.(method_name).actual;
    [r_val, p_val] = corr(p_scores, a_scores, 'Type', 'Spearman');
    scatter(p_scores, a_scores, 80, 'filled', 'MarkerFaceColor', show_colors{sm});
    hold on; lsline; hold off;
    xlabel('Predicted (Brain)','FontSize',11); ylabel('Actual (Foot)','FontSize',11);
    if p_val < 0.001
        title(sprintf('%s  r=%.3f, p<0.001', upper(method_name), r_val), 'FontSize',13,'FontWeight','bold');
    else
        title(sprintf('%s  r=%.3f, p=%.3f', upper(method_name), r_val, p_val), 'FontSize',13,'FontWeight','bold');
    end
    grid on; set(gca,'FontSize',10);
end
exportgraphics(fig4, sprintf('%sfigure4_brain_body_%s.png', out_dir, timestamp), 'Resolution', 300);


% -----------------------------------------------------------------------
% Fig D: GA Fitness Evolution across iterations (averaged over LOOCV folds)
% -----------------------------------------------------------------------
fprintf('Generating GA Fitness Evolution figure...\n');
methods_list = fieldnames(results);
n_m = length(methods_list);
fig_ga = figure('Color', 'w', 'Position', [50 50 min(1400, n_m*250) 450]);
tiledlayout(1, n_m, 'TileSpacing', 'compact', 'Padding', 'compact');
cmap_ga = lines(n_m);
for mg = 1:n_m
    mname = methods_list{mg};
    if ~isfield(results.(mname), 'ga_fitness_history'), continue; end
    hist_mat = results.(mname).ga_fitness_history;
    mean_fit = mean(hist_mat, 1);
    std_fit  = std(hist_mat, 0, 1);
    n_gens = length(mean_fit);
    gv = 1:n_gens;
    nexttile;
    fill([gv, fliplr(gv)], [mean_fit+std_fit, fliplr(mean_fit-std_fit)], cmap_ga(mg,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold on;
    plot(gv, mean_fit, '-', 'Color', cmap_ga(mg,:), 'LineWidth', 2.5);
    xlabel('GA Generation', 'FontSize', 11);
    ylabel('Best Fitness', 'FontSize', 11);
    title(sprintf('%s | Final: %.3f +/- %.3f', upper(mname), mean_fit(end), std_fit(end)), 'FontSize', 10, 'FontWeight', 'bold');
    ylim([0 1]); grid on; hold off;
end
n_folds_plot = size(results.(methods_list{1}).ga_fitness_history, 1);
sgtitle(sprintf('GA Fitness Evolution (mean +/- std, %d folds)', n_folds_plot), 'FontSize', 13, 'FontWeight', 'bold');
exportgraphics(fig_ga, sprintf('%sga_fitness_evolution_%s.png', out_dir, timestamp), 'Resolution', 300);
close(fig_ga);
fprintf('Saved GA fitness figure.\n');

% -----------------------------------------------------------------------

% -----------------------------------------------------------------------
% Fig E: GA Internal Training Scatter - one figure per fold (25 figures)
% Each figure shows all methods side by side for that fold
% -----------------------------------------------------------------------
fprintf('Generating per-fold GA internal scatter figures (25 folds)...\n');
methods_e = fieldnames(results);
n_me = length(methods_e);
mname0 = methods_e{1};
n_subs_p = length(results.(mname0).predicted);
n_train_p = n_subs_p - 1;
out_dir_folds = sprintf('%sga_internal_folds/', out_dir);
if ~exist(out_dir_folds, 'dir'), mkdir(out_dir_folds); end
for kf = 1:n_subs_p
    fold_range = (kf-1)*n_train_p + 1 : kf*n_train_p;
    fig_f = figure('Color', 'w', 'Position', [30 30 min(1600, n_me*260) 450]);
    tiledlayout(2, n_me, 'TileSpacing', 'compact', 'Padding', 'compact');
    cmap_f = lines(n_me);
    for me = 1:n_me
        mname = methods_e{me};
        if ~isfield(results.(mname),'ga_train_brain'), continue; end
        gb_f = results.(mname).ga_train_brain(fold_range);
        gf_f = results.(mname).ga_train_foot(fold_range);
        valid_f = gb_f ~= 0 | gf_f ~= 0;
        gb_f = gb_f(valid_f); gf_f = gf_f(valid_f);
        % Top row: GA internal training scatter for this fold
        nexttile;
        if length(gb_f) > 2
            [r_fi, p_fi] = corr(gb_f, gf_f, 'Type', 'Spearman');
            scatter(gb_f, gf_f, 50, [0.4 0.4 0.4], 'filled');
            hold on; lsline; hold off;
            if p_fi < 0.001
                title(sprintf('[%s] GA Train\nr=%.3f, p<0.001', upper(mname), r_fi), 'FontSize', 9, 'FontWeight', 'bold');
            else
                title(sprintf('[%s] GA Train\nr=%.3f, p=%.3f', upper(mname), r_fi, p_fi), 'FontSize', 9, 'FontWeight', 'bold');
            end
        else
            title(sprintf('[%s] GA Train\n(no valid tiles)', upper(mname)), 'FontSize', 9);
        end
        xlabel('Mean EEG (tiles)', 'FontSize', 8); ylabel('Mean Foot', 'FontSize', 8);
        grid on; set(gca, 'FontSize', 8);
        % Bottom row: LOOCV held-out scatter (full 25 pts), highlight current fold
        nexttile;
        ps = results.(mname).predicted;
        as = results.(mname).actual;
        [r_lo, p_lo] = corr(ps, as, 'Type', 'Spearman');
        other_idx = setdiff(1:n_subs_p, kf);
        scatter(ps(other_idx), as(other_idx), 40, cmap_f(me,:), 'filled', 'MarkerFaceAlpha', 0.3);
        hold on;
        scatter(ps(kf), as(kf), 80, cmap_f(me,:), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        lsline; hold off;
        if p_lo < 0.001
            title(sprintf('[%s] LOOCV\nr=%.3f, p<0.001', upper(mname), r_lo), 'FontSize', 9, 'FontWeight', 'bold');
        else
            title(sprintf('[%s] LOOCV\nr=%.3f, p=%.3f', upper(mname), r_lo, p_lo), 'FontSize', 9, 'FontWeight', 'bold');
        end
        xlabel('Predicted', 'FontSize', 8); ylabel('Actual', 'FontSize', 8);
        grid on; set(gca, 'FontSize', 8);
    end
    sgtitle(sprintf('Fold %02d/%d — GA Internal Training (top) vs LOOCV Prediction (bottom)', kf, n_subs_p), 'FontSize', 11, 'FontWeight', 'bold');
    fname_fold = sprintf('%sga_internal_fold_%02d_%s.png', out_dir_folds, kf, timestamp);
    exportgraphics(fig_f, fname_fold, 'Resolution', 200);
    close(fig_f);
    fprintf('  Saved fold %02d\n', kf);
end
fprintf('All per-fold figures saved to: %s\n', out_dir_folds);
fprintf('=== STEP 6b COMPLETE ===\n');
