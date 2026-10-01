function [results, loocv_save_path] = LOOCV_Corrected(all_sub_diffs, foot_features, foot_labels, config, save_dir, methods_to_run, batch_id)
% LOOCV_Corrected.m




% Performs RIGOROUS Leave-One-Out Cross-Validation (LOOCV) for GA Feature Selection

% Optional: specify custom save directory for results
if nargin < 5 || isempty(save_dir)
    save_dir = fullfile(config.out_dir, 'ga_feature');
end

% Optional: specify which methods to run (default: all 7)
if nargin < 6 || isempty(methods_to_run)
    methods_to_run = {'spearman', 'pca', 'robust', 'pls', 'lasso', 'ridge', 'svr'};
end

% Optional: batch ID to group files from the same run session
if nargin < 7 || isempty(batch_id)
    batch_id = '';
end

fprintf('=== Starting Rigorous LOOCV (Corrected for Data Leakage) ===\n');
fprintf('Save directory: %s\n', save_dir);

results = struct();
if ~exist('all_sub_diffs', 'var')
    error('Please load "all_sub_diffs" into the workspace first.');
end

sub_list = fieldnames(all_sub_diffs);
n_subs = length(sub_list);
tile_t = all_sub_diffs.sub_1.t;
tile_f = all_sub_diffs.sub_1.f;

fprintf('[DIAG] n_subs=%d\n', n_subs);
fprintf('[DIAG] all_sub_diffs fields: %s\n', strjoin(sub_list(1:min(3,end))', ', '));
eeg_fields_check = fieldnames(all_sub_diffs.(sub_list{1}).eeg);
fprintf('[DIAG] all_sub_diffs.sub_1.eeg fields: %s\n', strjoin(eeg_fields_check', ', '));

fprintf('Extracting RAW features from all_sub_diffs...\n');
eeg_feature_raw = zeros(n_subs, 1200);

get_raw_vec = @(sub, domain, field) reshape(all_sub_diffs.(sub).(domain).(field), 1, []);

for i = 1:n_subs
    sub_name = sub_list{i};
    eeg_feature_raw(i, 1:400)    = get_raw_vec(sub_name, 'eeg', 'expected_vs_expected_duplicated');
    eeg_feature_raw(i, 401:800)  = get_raw_vec(sub_name, 'eeg', 'unexpected_vs_expected');
    eeg_feature_raw(i, 801:1200) = get_raw_vec(sub_name, 'eeg', 'right_after_vs_right_before');
end

fprintf('[DIAG] eeg_feature_raw: size=[%s], anyNaN=%d, allZero=%d, std=%.4f\n', ...
    num2str(size(eeg_feature_raw)), any(isnan(eeg_feature_raw(:))), all(eeg_feature_raw(:)==0), std(eeg_feature_raw(:)));
fprintf('[DIAG] foot_features: size=[%s], anyNaN=%d, allZero=%d, std=%.4f\n', ...
    num2str(size(foot_features)), any(isnan(foot_features(:))), all(foot_features(:)==0), std(foot_features(:)));

eeg_data = eeg_feature_raw;
foot_data = foot_features;

eeg_data = remove_unwanted_area(all_sub_diffs, eeg_data, tile_f, [10, 30]);

fprintf('[DIAG] eeg_data after remove_unwanted_area: size=[%s], nonzero_cols=%d/%d, std=%.4f\n', ...
    num2str(size(eeg_data)), sum(any(eeg_data ~= 0)), size(eeg_data,2), std(eeg_data(:)));

% Row-wise normalisation: each subject's features normalised to mean=0, std=1
% using only that subject's own values → no cross-subject leakage.
% This preserves between-subject column variance, so ttest inside the fold
% can detect features whose group mean is significantly non-zero.
eeg_data  = normalize(eeg_data,  2);
foot_data = normalize(foot_data, 2);

fprintf('[DIAG] eeg_data after row-normalise: std=%.4f  foot_data after row-normalise: std=%.4f\n', ...
    std(eeg_data(:)), std(foot_data(:)));

methods = methods_to_run;  % Use parameter passed in (or default to all 7)
results = struct();

fprintf('Running methods: %s\n', strjoin(methods, ', '));
% Colormap for ERSP
cmap_bwr = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ...
            ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];
for m = 1:length(methods)
    current_method = methods{m};
    fprintf('\n\n========== RUNNING LOOCV WITH METHOD: %s ==========\n', current_method);
    
    predicted_scores = zeros(n_subs, 1);
    actual_scores = zeros(n_subs, 1);
    
    % Accumulate GA masks across all folds for consensus
    mask_x_accum = zeros(1, size(eeg_data, 2));
    mask_y_accum = zeros(1, size(foot_data, 2));
    ga_fitness_all_folds = [];  % will grow to [n_subs x max_gens]
    ga_brain_train_all = [];    % concat of 24 training brain scores per fold
    ga_foot_train_all  = [];    % concat of 24 training foot scores per fold
    
    for i = 1:n_subs
        fprintf('Method %s | Fold %d/%d... ', current_method, i, n_subs);
        
        test_idx = i;
        train_idx = setdiff(1:n_subs, i);
        
        eeg_train = eeg_data(train_idx, :);
        foot_train = foot_data(train_idx, :);
        
        eeg_test = eeg_data(test_idx, :);
        foot_test = foot_data(test_idx, :);

        % ttest on row-normalised training data: detects columns where the group mean
        % is significantly non-zero (i.e. features that consistently differ from
        % each subject's own baseline across the training fold).
        std_control = std(eeg_train(:, 1:400));
        std_surprise = std(eeg_train(:, 401:800));
        std_learning = std(eeg_train(:, 801:1200));
        fprintf('[DIAG F%d] eeg_train std by phase: Control=%.4f, Surprise=%.4f, Learning=%.4f\n', i, ...
            std_control, std_surprise, std_learning);

        % Uniform t-test threshold for all phases — no thumb on the scale
        [h_eeg, ~] = ttest(eeg_train, 0, 'Alpha', 0.05);

        [h_foot, ~] = ttest(foot_train, 0, 'Alpha', 0.05);

        % ttest returns NaN for zero-variance columns; treat as non-significant
        h_eeg(isnan(h_eeg)) = 0;
        h_foot(isnan(h_foot)) = 0;

        if i == 1
            fprintf('[DIAG F%d] eeg_train: size=[%s] std=%.4f\n', i, num2str(size(eeg_train)), std(eeg_train(:)));
            fprintf('[DIAG F%d] ttest: eeg sig cols=%d/%d, foot sig cols=%d/%d\n', i, sum(h_eeg), length(h_eeg), sum(h_foot), length(h_foot));
        end

        eeg_train_filtered = eeg_train;
        eeg_train_filtered(:, h_eeg == 0) = 0;

        foot_train_filtered = foot_train;
        foot_train_filtered(:, h_foot == 0) = 0;

        if i == 1
            fprintf('[DIAG F%d] ttest: Control(α=0.001)=%d/400, Surprise(α=0.05)=%d/400, Learning(α=0.05)=%d/400, foot=%d/40\n', i, ...
                sum(h_eeg(1:400)), sum(h_eeg(401:800)), sum(h_eeg(801:1200)), sum(h_foot));
            fprintf('[DIAG F%d] eeg_train_filtered nonzero cols=%d\n', i, sum(any(eeg_train_filtered ~= 0)));
            fprintf('[DIAG F%d] foot_train_filtered nonzero cols=%d\n', i, sum(any(foot_train_filtered ~= 0)));

            % --- PLOT RAW EEG DATA + T-TEST RESULTS BY PHASE ---
            n_tile_f = length(tile_f);
            n_tile_t = length(tile_t);

            phase_names = {'Control', 'Surprise', 'Learning'};
            phase_ranges = {1:400, 401:800, 801:1200};

            fig = figure('Color', 'w', 'Position', [100 100 2100 800]);
            tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

            % Row 1: Raw EEG mean across training fold
            for p = 1:3
                nexttile;
                phase_range = phase_ranges{p};
                eeg_phase = eeg_train(:, phase_range);  % [24 subjects × 400 tiles]
                eeg_mean = mean(eeg_phase, 1);  % [1 × 400] average across subjects
                eeg_matrix = reshape(eeg_mean, n_tile_f, n_tile_t);  % Reshape to 2D grid

                imagesc(tile_t, tile_f, eeg_matrix);
                colormap(gca, cmap_bwr);
                caxis_lim = max(abs(eeg_matrix(:)));
                caxis([-caxis_lim, caxis_lim]);  % Center colormap at 0
                colorbar;

                title(sprintf('%s (Raw Mean): min=%.4f, max=%.4f, std=%.4f', phase_names{p}, min(eeg_mean), max(eeg_mean), std(eeg_mean)));
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
            end

            % Row 2: T-test significance maps
            for p = 1:3
                nexttile;
                phase_range = phase_ranges{p};
                h_phase = h_eeg(phase_range);  % Binary: 1 if significant, 0 otherwise
                h_matrix = reshape(h_phase, n_tile_f, n_tile_t);  % Reshape to 2D grid

                imagesc(tile_t, tile_f, h_matrix);
                colormap(gca, 'gray');
                caxis([0 1]);
                colorbar;

                title(sprintf('%s (T-test): %d/%d tiles significant (α=0.001)', phase_names{p}, sum(h_phase), length(h_phase)));
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
            end

            ttest_fig_dir = fullfile(config.out_dir, 'all_in_one_pic', 'ttest_maps');
            if ~exist(ttest_fig_dir, 'dir')
                mkdir(ttest_fig_dir);
            end
            ttest_fig_path = fullfile(ttest_fig_dir, sprintf('raw_vs_ttest_phase_maps_fold_%02d.png', i));
            exportgraphics(fig, ttest_fig_path, 'Resolution', 150);
            close(fig);

            fprintf('[TTEST VIZ] Saved raw EEG + t-test phase maps to: %s\n', ttest_fig_path);
        end

        [best_fitness_history, cfg_ga, best_solution] = ga_iteration(eeg_train_filtered, foot_train_filtered, current_method);
        ga_fitness_all_folds = [ga_fitness_all_folds; best_fitness_history];
        
        num_genes_x = cfg_ga.num_genes_x;
        ga_mask_x = best_solution(1 : num_genes_x);
        ga_mask_y = best_solution(num_genes_x + 1 : end);

        if i == 1
            fprintf('[DIAG F%d] GA best_fitness=%.4f, mask_x selected=%d/%d, mask_y selected=%d/%d\n', ...
                i, cfg_ga.max_fitness_ever, sum(ga_mask_x), length(ga_mask_x), sum(ga_mask_y), length(ga_mask_y));
        end
        % Capture GA training scatter (what GA's correlation actually looks like)
        if sum(ga_mask_x) > 0 && sum(ga_mask_y) > 0
            ga_brain_train_all = [ga_brain_train_all; mean(eeg_train_filtered(:, ga_mask_x==1), 2)];
            ga_foot_train_all  = [ga_foot_train_all;  mean(foot_train_filtered(:, ga_mask_y==1), 2)];
        else
            ga_brain_train_all = [ga_brain_train_all; zeros(length(train_idx), 1)];
            ga_foot_train_all  = [ga_foot_train_all;  zeros(length(train_idx), 1)];
        end
        
        eeg_test_filtered = eeg_test;
        eeg_test_filtered(:, h_eeg == 0) = 0;
        foot_test_filtered = foot_test;
        foot_test_filtered(:, h_foot == 0) = 0;
        
        if strcmp(current_method, 'pca')
            if sum(ga_mask_x) > 0
                train_X = eeg_train_filtered(:, ga_mask_x == 1);
                [coeff_X, ~, ~, ~, ~, mu_X] = pca(train_X, 'NumComponents', 1);
                test_X = eeg_test_filtered(:, ga_mask_x == 1);  % FIX: Use (:, mask) to preserve row vector
                if ~isempty(coeff_X)
                    brain_score = (test_X - mu_X) * coeff_X;
                else
                    brain_score = 0;
                end
            else
                brain_score = 0;
            end

            if sum(ga_mask_y) > 0
                train_Y = foot_train_filtered(:, ga_mask_y == 1);
                [coeff_Y, ~, ~, ~, ~, mu_Y] = pca(train_Y, 'NumComponents', 1);
                test_Y = foot_test_filtered(:, ga_mask_y == 1);  % FIX: Use (:, mask) to preserve row vector
                if ~isempty(coeff_Y)
                    foot_score = (test_Y - mu_Y) * coeff_Y;
                else
                    foot_score = 0;
                end
            else
                foot_score = 0;
            end
            
        elseif strcmp(current_method, 'pls')
            if sum(ga_mask_x) > 0 && sum(ga_mask_y) > 0
                train_X = eeg_train_filtered(:, ga_mask_x == 1);
                train_Y = foot_train_filtered(:, ga_mask_y == 1);
                [~, Y_loadings, ~, ~, ~, ~, ~, stats] = plsregress(train_X, train_Y, 1);

                test_X = eeg_test_filtered(:, ga_mask_x == 1);  % FIX: Use (:, mask) to preserve row vector
                if ~isempty(stats.W)
                    X_mean = mean(train_X);
                    brain_score = (test_X - X_mean) * stats.W;
                else
                    brain_score = 0;
                end

                test_Y = foot_test_filtered(:, ga_mask_y == 1);  % FIX: Use (:, mask) to preserve row vector
                if ~isempty(Y_loadings)
                    Y_mean = mean(train_Y);
                    u_weights = pinv(Y_loadings');
                    foot_score = (test_Y - Y_mean) * u_weights;
                else
                    foot_score = 0;
                end
            else
                brain_score = 0;
                foot_score = 0;
            end
            
        elseif strcmp(current_method, 'ridge')
            if sum(ga_mask_x) > 0 && sum(ga_mask_y) > 0
                train_X = eeg_train_filtered(:, ga_mask_x == 1);
                train_Y_1D = mean(foot_train_filtered(:, ga_mask_y == 1), 2);
                try
                    weights = ridge(train_Y_1D, train_X, 0.1, 0);
                    % ridge(y, X, k, 0) returns [intercept; coeff_1; ...; coeff_p]
                    if isrow(weights), weights = weights'; end
                    intercept = weights(1);
                    feature_weights = weights(2:end);

                    test_X = eeg_test_filtered(:, ga_mask_x == 1);
                    test_Y_1D = mean(foot_test_filtered(:, ga_mask_y == 1), 2);
                    brain_score = test_X * feature_weights + intercept;
                    foot_score = test_Y_1D;
                catch
                    brain_score = 0;
                    foot_score = 0;
                end
            else
                brain_score = 0; foot_score = 0;
            end
            
        elseif strcmp(current_method, 'lasso')
            if sum(ga_mask_x) > 0 && sum(ga_mask_y) > 0
                train_X = eeg_train_filtered(:, ga_mask_x == 1);
                train_Y_1D = mean(foot_train_filtered(:, ga_mask_y == 1), 2);
                try
                    [weights, fitInfo] = lasso(train_X, train_Y_1D, 'Lambda', 0.05);
                    test_X = eeg_test_filtered(:, ga_mask_x == 1);
                    test_Y_1D = mean(foot_test_filtered(:, ga_mask_y == 1), 2);
                    brain_score = test_X * weights + fitInfo.Intercept;
                    foot_score = test_Y_1D;
                catch
                    brain_score = 0;
                    foot_score = 0;
                end
            else
                brain_score = 0; foot_score = 0;
            end
            
        elseif strcmp(current_method, 'svr')
            if sum(ga_mask_x) > 0 && sum(ga_mask_y) > 0
                train_X = eeg_train_filtered(:, ga_mask_x == 1);
                train_Y_1D = mean(foot_train_filtered(:, ga_mask_y == 1), 2);
                try
                    svm_mdl = fitrsvm(train_X, train_Y_1D, 'KernelFunction', 'linear', 'Standardize', true);
                    test_X = eeg_test_filtered(:, ga_mask_x == 1);
                    test_Y_1D = mean(foot_test_filtered(:, ga_mask_y == 1), 2);
                    brain_score = predict(svm_mdl, test_X);
                    foot_score = test_Y_1D;
                catch
                    brain_score = 0;
                    foot_score = 0;
                end
            else
                brain_score = 0; foot_score = 0;
            end
            
        else
            if sum(ga_mask_x) > 0
                brain_score = mean(eeg_test_filtered(:, ga_mask_x == 1), 2);
            else
                brain_score = 0;
            end

            if sum(ga_mask_y) > 0
                foot_score = mean(foot_test_filtered(:, ga_mask_y == 1), 2);
            else
                foot_score = 0;
            end
        end
        fprintf('[DIAG F%d] brain_score size=[%s] val=%s | foot_score size=[%s] val=%s\n', ...
            i, num2str(size(brain_score)), num2str(brain_score(:)'), num2str(size(foot_score)), num2str(foot_score(:)'));

        predicted_scores(i) = brain_score;
        actual_scores(i) = foot_score;

        % Accumulate GA mask from every fold
        mask_x_accum = mask_x_accum + ga_mask_x;
        mask_y_accum = mask_y_accum + ga_mask_y;

        fprintf('FIT: %.3f\n', cfg_ga.max_fitness_ever);
    end
    
    fprintf('[DIAG] %s predicted_scores: %s\n', current_method, num2str(predicted_scores(:)', '%.3f '));
    fprintf('[DIAG] %s actual_scores:    %s\n', current_method, num2str(actual_scores(:)', '%.3f '));
    fprintf('[DIAG] %s std(pred)=%.4f std(actual)=%.4f\n', current_method, std(predicted_scores), std(actual_scores));

    results.(current_method).predicted = predicted_scores;
    results.(current_method).actual = actual_scores;
    
    % Compute consensus mask: tiles selected in >=50% of folds
    results.(current_method).mask_x = double(mask_x_accum >= (n_subs / 2));
    results.(current_method).mask_y = double(mask_y_accum >= (n_subs / 2));
    results.(current_method).mask_x_counts = mask_x_accum;
    results.(current_method).mask_y_counts = mask_y_accum;
    results.(current_method).ga_fitness_history = ga_fitness_all_folds;  % [n_folds x max_gens]
    results.(current_method).ga_train_brain = ga_brain_train_all;  % [n_folds*(n_subs-1) x 1]
    results.(current_method).ga_train_foot  = ga_foot_train_all;   % [n_folds*(n_subs-1) x 1]
    fprintf('Consensus mask: %d/%d EEG tiles, %d/%d Foot tiles selected in >=50%% of folds\n', ...
        sum(results.(current_method).mask_x), length(mask_x_accum), ...
        sum(results.(current_method).mask_y), length(mask_y_accum));
end

num_methods = length(methods);
fig = figure('Color', 'w', 'Position', [50 50 1800 900]);
tiledlayout(2, ceil(num_methods/2), 'Padding', 'compact', 'TileSpacing', 'compact'); 

plot_scatter = @(res, name, col) ...
    (scatter(res.predicted, res.actual, 80, 'filled', 'MarkerFaceColor', col));

colors = lines(num_methods);

for m = 1:num_methods
    method_name = methods{m};
    nexttile;
    p_scores = results.(method_name).predicted;
    a_scores = results.(method_name).actual;
    [r, p] = corr(p_scores, a_scores, 'Type', 'Spearman');
    plot_scatter(results.(method_name), '', colors(m,:)); lsline;
    title({sprintf('%d. %s', m, upper(method_name)), sprintf('r=%.3f, p=%.3f', r, p)});
    xlabel('Predicted'); ylabel('Actual'); grid on;
end

timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm'));
fprintf('Saving results with timestamp: %s\n', timestamp);

run_info = struct();
run_info.timestamp = timestamp;
run_info.methods = methods;
run_info.reduced_freq_range_eeg = [10 30]; 
run_info.reduced_freq_range_foot = [0.1 15];
run_info.feature_phases = {'Expected', 'Unexpected', 'Adaptation'};
run_info.ga_config_example = cfg_ga; 

if isempty(batch_id)
    save_filename = fullfile(save_dir, sprintf('results_LOOCV_%s.mat', timestamp));
    fig_filename  = fullfile(save_dir, sprintf('loocv_comparison_methods_%s.png', timestamp));
else
    save_filename = fullfile(save_dir, sprintf('results_LOOCV_%s_%s.mat', timestamp, batch_id));
    fig_filename  = fullfile(save_dir, sprintf('loocv_comparison_methods_%s_%s.png', timestamp, batch_id));
end
run_info.tile_t = tile_t; run_info.tile_f = tile_f;
run_info.batch_id = batch_id;
save(save_filename, 'results', 'run_info', 'eeg_data', 'foot_data', 'foot_labels', 'config');
loocv_save_path = save_filename;
fprintf('Analysis saved to: %s\n', save_filename);
exportgraphics(fig, fig_filename, 'Resolution', 300);

end
