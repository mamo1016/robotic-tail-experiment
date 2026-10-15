function results = run_loocv_core(eeg_data, foot_data)
% run_loocv_core — Lightweight LOOCV for permutation testing.
%
% Runs the full 7-method GA-LOOCV on pre-normalized data.
% NO figures, NO file saves, NO clc. Used by permutation_test_loocv.m.
%
% Inputs:
%   eeg_data  — [N_subs x 1200] pre-normalized EEG feature matrix
%   foot_data — [N_subs x 40]   pre-normalized foot feature matrix
%                (already row-shuffled by the caller for permutation)
%
% Output:
%   results — struct with one field per method containing .predicted and .actual

n_subs  = size(eeg_data, 1);
methods = {'spearman', 'pca', 'robust', 'pls', 'lasso', 'ridge', 'svr'};
results = struct();

for m = 1:length(methods)
    current_method   = methods{m};
    predicted_scores = zeros(n_subs, 1);
    actual_scores    = zeros(n_subs, 1);

    for i = 1:n_subs
        train_idx = setdiff(1:n_subs, i);

        eeg_train  = eeg_data(train_idx, :);
        foot_train = foot_data(train_idx, :);
        eeg_test   = eeg_data(i, :);
        foot_test  = foot_data(i, :);

        % T-test filter (same as LOOCV_Corrected)
        [h_eeg,  ~] = ttest(eeg_train,  0, 'Alpha', 0.05);
        [h_foot, ~] = ttest(foot_train, 0, 'Alpha', 0.05);

        eeg_train_f  = eeg_train;  eeg_train_f(:, h_eeg  == 0) = 0;
        foot_train_f = foot_train; foot_train_f(:, h_foot == 0) = 0;

        % GA feature selection
        [~, cfg_ga, best_solution] = ga_iteration(eeg_train_f, foot_train_f, current_method);

        num_genes_x = cfg_ga.num_genes_x;
        ga_mask_x   = best_solution(1 : num_genes_x);
        ga_mask_y   = best_solution(num_genes_x + 1 : end);

        % Apply t-test filter to test subject
        eeg_test_f  = eeg_test;  eeg_test_f(h_eeg  == 0) = 0;
        foot_test_f = foot_test; foot_test_f(h_foot == 0) = 0;

        % --- Method-specific prediction (mirrors LOOCV_Corrected exactly) ---
        if strcmp(current_method, 'pca')
            if sum(ga_mask_x) > 0
                train_X = eeg_train_f(:, ga_mask_x == 1);
                [coeff_X, ~, ~, ~, ~, mu_X] = pca(train_X, 'NumComponents', 1);
                test_X = eeg_test_f(ga_mask_x == 1);
                if ~isempty(coeff_X)
                    brain_score = (test_X - mu_X) * coeff_X;
                else
                    brain_score = 0;
                end
            else
                brain_score = 0;
            end
            if sum(ga_mask_y) > 0
                train_Y = foot_train_f(:, ga_mask_y == 1);
                [coeff_Y, ~, ~, ~, ~, mu_Y] = pca(train_Y, 'NumComponents', 1);
                test_Y = foot_test_f(ga_mask_y == 1);
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
                train_X = eeg_train_f(:, ga_mask_x == 1);
                train_Y = foot_train_f(:, ga_mask_y == 1);
                [~, Y_loadings, ~, ~, ~, ~, ~, stats] = plsregress(train_X, train_Y, 1);
                test_X = eeg_test_f(ga_mask_x == 1);
                if ~isempty(stats.W)
                    X_mean = mean(train_X);
                    brain_score = (test_X - X_mean) * stats.W;
                else
                    brain_score = 0;
                end
                test_Y = foot_test_f(ga_mask_y == 1);
                if ~isempty(Y_loadings)
                    Y_mean = mean(train_Y);
                    u_weights = pinv(Y_loadings');
                    foot_score = (test_Y - Y_mean) * u_weights;
                else
                    foot_score = 0;
                end
            else
                brain_score = 0; foot_score = 0;
            end

        elseif strcmp(current_method, 'ridge')
            if sum(ga_mask_x) > 0 && sum(ga_mask_y) > 0
                train_X    = eeg_train_f(:, ga_mask_x == 1);
                train_Y_1D = mean(foot_train_f(:, ga_mask_y == 1), 2);
                try
                    weights     = ridge(train_Y_1D, train_X, 0.1, 0);
                    test_X      = eeg_test_f(ga_mask_x == 1);
                    brain_score = [1 test_X] * weights;
                    foot_score  = mean(foot_test_f(ga_mask_y == 1));
                catch
                    brain_score = 0; foot_score = 0;
                end
            else
                brain_score = 0; foot_score = 0;
            end

        elseif strcmp(current_method, 'lasso')
            if sum(ga_mask_x) > 0 && sum(ga_mask_y) > 0
                train_X    = eeg_train_f(:, ga_mask_x == 1);
                train_Y_1D = mean(foot_train_f(:, ga_mask_y == 1), 2);
                try
                    [weights, fitInfo] = lasso(train_X, train_Y_1D, 'Lambda', 0.05);
                    test_X      = eeg_test_f(ga_mask_x == 1);
                    brain_score = test_X * weights + fitInfo.Intercept;
                    foot_score  = mean(foot_test_f(ga_mask_y == 1));
                catch
                    brain_score = 0; foot_score = 0;
                end
            else
                brain_score = 0; foot_score = 0;
            end

        elseif strcmp(current_method, 'svr')
            if sum(ga_mask_x) > 0 && sum(ga_mask_y) > 0
                train_X    = eeg_train_f(:, ga_mask_x == 1);
                train_Y_1D = mean(foot_train_f(:, ga_mask_y == 1), 2);
                try
                    svm_mdl     = fitrsvm(train_X, train_Y_1D, 'KernelFunction', 'linear', 'Standardize', true);
                    test_X      = eeg_test_f(ga_mask_x == 1);
                    brain_score = predict(svm_mdl, test_X);
                    foot_score  = mean(foot_test_f(ga_mask_y == 1));
                catch
                    brain_score = 0; foot_score = 0;
                end
            else
                brain_score = 0; foot_score = 0;
            end

        else  % spearman, robust, and any other mean-based method
            if sum(ga_mask_x) > 0
                brain_score = mean(eeg_test_f(ga_mask_x == 1));
            else
                brain_score = 0;
            end
            if sum(ga_mask_y) > 0
                foot_score = mean(foot_test_f(ga_mask_y == 1));
            else
                foot_score = 0;
            end
        end

        predicted_scores(i) = brain_score;
        actual_scores(i)    = foot_score;
    end

    results.(current_method).predicted = predicted_scores;
    results.(current_method).actual    = actual_scores;
end
end
