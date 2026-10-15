config = configure_parameters();
out_dir = config.out_dir;
save_dir      = fullfile(out_dir, 'ga_feature');

checkpoint_file = fullfile(save_dir, 'permutation_checkpoint.mat');
load(checkpoint_file, 'r_null', 'perm_resume');
n_done = perm_resume - 1;   % perm_resume is the NEXT one to run, so completed = perm_resume-1
fprintf('\n>>> Checkpoint loaded: %d permutations completed <<<\n', n_done);

% r_null is [N_PERM x N_methods] — summarise each method
r_done = r_null(1:n_done, :);   % rows completed so far
n_methods = size(r_done, 2);
method_names = {'spearman','pca','robust','pls','lasso','ridge','svr'};

observed_r = 0.99;   % update this if your observed r differs

fprintf('\n%-10s  %8s  %8s  %8s  %8s\n', 'Method', 'Max_null', 'Med_null', '95pct', 'p_emp');
fprintf('%s\n', repmat('-', 1, 52));
for m = 1:n_methods
    col = r_done(:, m);
    col = col(~isnan(col));
    if isempty(col), continue; end
    p_emp = sum(col >= observed_r) / numel(col);
    if strcmp(method_names{min(m,end)}, method_names{m})
        mname = method_names{m};
    else
        mname = sprintf('method%d', m);
    end
    fprintf('%-10s  %8.4f  %8.4f  %8.4f  %8.4f\n', ...
        mname, max(col), median(col), prctile(col,95), p_emp);
end

% Plot null distribution for each method
n_cols = min(4, n_methods);
n_rows = ceil(n_methods / n_cols);
figure('Color','w','Position',[50 50 300*n_cols 260*n_rows]);
tl = tiledlayout(n_rows, n_cols, 'TileSpacing','compact','Padding','compact');
for m = 1:n_methods
    col = r_done(:, m);
    col = col(~isnan(col));
    if isempty(col), continue; end
    nexttile;
    histogram(col, 30, 'FaceColor',[0.6 0.6 0.8], 'EdgeColor','none');
    hold on;
    xline(observed_r, 'r--', 'LineWidth', 2);
    xline(prctile(col,95), 'k--', 'LineWidth', 1.2);
    hold off;
    p_emp = sum(col >= observed_r) / numel(col);
    if p_emp == 0
        p_str = sprintf('p < %.4f', 1/n_done);
    else
        p_str = sprintf('p = %.4f', p_emp);
    end
    if m <= numel(method_names)
        mname = method_names{m};
    else
        mname = sprintf('method%d', m);
    end
    title(sprintf('%s | %s', upper(mname), p_str), 'FontSize', 9);
    xlabel('Permuted r'); ylabel('Count');
end
title(tl, sprintf('Partial null distribution (%d / 1000 permutations)', n_done), 'FontSize', 11);

%% --- Save figure ---
paper_fig_dir = 'figures';
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM');
fig_path = fullfile(paper_fig_dir, ...
    sprintf('permutation_null_distribution_partial_%dperms_%s.png', n_done, timestamp));
exportgraphics(gcf, fig_path, 'Resolution', 300);
fprintf('\nFigure saved to: %s\n', fig_path);

%% --- Save results .mat ---
% Collect summary stats into a struct for each method
p_empirical = NaN(1, n_methods);
r_max_null  = NaN(1, n_methods);
r_med_null  = NaN(1, n_methods);
r_95_null   = NaN(1, n_methods);
for m = 1:n_methods
    col = r_done(:, m);
    col = col(~isnan(col));
    if isempty(col), continue; end
    p_empirical(m) = sum(col >= observed_r) / numel(col);
    r_max_null(m)  = max(col);
    r_med_null(m)  = median(col);
    r_95_null(m)   = prctile(col, 95);
end

mat_path = fullfile(save_dir, ...
    sprintf('permutation_partial_%dperms_%s.mat', n_done, timestamp));
save(mat_path, 'r_done', 'r_null', 'n_done', 'observed_r', ...
    'p_empirical', 'r_max_null', 'r_med_null', 'r_95_null', ...
    'method_names', '-v7.3');
fprintf('Results saved to: %s\n', mat_path);
fprintf('\n=== DONE ===\n');