function LOOCV(all_sub_diffs)
% perform_LOOCV.m
% Performs Leave-One-Out Cross-Validation (LOOCV) for GA Feature Selection
% Author: Jane (OpenClaw)
% Usage: Run this AFTER you have generated 'all_sub_diffs' in the main script.

% clearvars -except all_sub_diffs config tile_t tile_f; 
clc;

fprintf('=== Starting Leave-One-Out Cross-Validation (LOOCV) ===\n');
tile_t = all_sub_diffs.sub_1.t;
tile_f = all_sub_diffs.sub_1.f;
% 1. Prepare Feature Matrix (Same as Main Script)
% Ensure 'all_sub_diffs' exists in workspace or load it
if ~exist('all_sub_diffs', 'var')
    try
        % Try loading from default path if not in workspace
        % Adjust path if needed
        load("output\ga_feature\sub_tile.mat"); 
        load("output\ga_feature\sub_tile_foot.mat");
        % Re-running feature extraction logic if needed...
        % (Ideally, save 'all_sub_diffs' before running this)
        error('Variable all_sub_diffs not found. Please run main script up to feature extraction.');
    catch
        error('Please load "all_sub_diffs" into the workspace first.');
    end
end

% Extract features (Full Dataset)
[eeg_feature_all, foot_feature_all] = feature_ext_sig(all_sub_diffs);
eeg_feature_U_vs_E = eeg_feature_all(:,401:800);
foot_feature_U_vs_E = foot_feature_all(:,401:800);
% Apply same cleaning/masking as main script
% eeg_feature_all(:,1:400) = 0;
% eeg_feature_all(:,801:1200) = 0;
% foot_feature_all(:,1:400) = 0;
% foot_feature_all(:,801:1200) = 0;

eeg_feature_U_vs_E = remove_unwanted_area(all_sub_diffs, eeg_feature_U_vs_E, tile_f, [5, 40]);
foot_feature_U_vs_E = remove_unwanted_area(all_sub_diffs, foot_feature_U_vs_E, tile_f, [0, 15]);

n_subs = size(eeg_feature_U_vs_E, 1);
predicted_foot_scores = zeros(n_subs, 1);
actual_foot_scores = zeros(n_subs, 1);

fprintf('Total Subjects: %d\n', n_subs);

% 2. LOOCV Loop
for i = 1:n_subs
    fprintf('Fold %d/%d: Leaving out Subject %d...', i, n_subs, i);
    
    % A. Define Training Set (N-1)
    train_idx = setdiff(1:n_subs, i);
    eeg_train = eeg_feature_U_vs_E(train_idx, :);
    eeg_train(:,eeg_train(1,:)==0) = [];
    foot_train = foot_feature_U_vs_E(train_idx, :);
    foot_train(:,foot_train(1,:)==0) = [];
    
    % B. Run GA on Training Set
    % Note: Reduced generations for speed during LOOCV? Or full?
    % Let's use full to be rigorous, but suppress output.
    [~, cfg_ga, best_solution] = ga_iteration(eeg_train, foot_train);
    
    % C. Extract Best Mask from Training
    num_genes_x = cfg_ga.num_genes_x;
    best_mask_x = best_solution(1 : num_genes_x);
    best_mask_y = best_solution(num_genes_x + 1 : end);
    
    % D. Test on Left-Out Subject
    eeg_test = eeg_feature_U_vs_E(i, :); % here
    eeg_test(:,eeg_test(1,:)==0) = [];
    foot_test = foot_feature_U_vs_E(i, :);
    foot_test(:,foot_test(1,:)==0) = [];
    
    % Apply Training Mask to Test Data
    % Score = Mean of selected tiles
    % Handle empty mask case (avoid NaN)
    if sum(best_mask_x) > 0
        brain_score = mean(eeg_test(best_mask_x == 1));
    else
        brain_score = 0; 
    end
    
    if sum(best_mask_y) > 0
        foot_score = mean(foot_test(best_mask_y == 1));
    else
        foot_score = 0;
    end
    
    % Store Results
    predicted_foot_scores(i) = brain_score; % This is the "Neural Prediction"
    actual_foot_scores(i) = foot_score;     % This is the Ground Truth
    
    fprintf(' Done. (Corr: %.2f)\n', cfg_ga.max_fitness_ever);
end

[~,idx]=find(eeg_feature_U_vs_E(1,:)~=0);
idx = idx(logical(best_mask_x));
plot_mask = zeros(size(eeg_feature_U_vs_E(1,:)));
plot_mask(idx) = 1;
[selected_t, selected_f] = mask_visualiser(plot_mask,tile_f,tile_t);
% -------------------
dif=mean(eeg_feature_U_vs_E,1);
dif = reshape(dif, [length(tile_f) length(tile_t)]);

fig_ersp_ave = figure('Color', 'w'); % 'Name' sets the window title
visualisation_ersp(dif, fig_ersp_ave, tile_t, tile_f, [-0.2 0.2]); %title(title_str);
hold on
plot(selected_t, selected_f, 'k.', 'MarkerSize', 5, 'LineWidth', 1.5);
hold off

% -----------------------
% 3. Final Validation Stats
[r_valid, p_valid] = corr(predicted_foot_scores, actual_foot_scores, 'Type', 'Spearman');

fprintf('\n=== CROSS-VALIDATION RESULTS ===\n');
fprintf('Spearman Correlation (LOOCV): r = %.3f\n', r_valid);
fprintf('P-Value: p = %.4f\n', p_valid);

% 4. Plot
figure('Color', 'w');
scatter(predicted_foot_scores, actual_foot_scores, 100, 'filled', 'MarkerFaceColor', 'b');
lsline;
xlabel('Predicted Brain Score (Left-Out Subject)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Actual Foot Score', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('LOOCV Validation (r=%.3f, p=%.3f)', r_valid, p_valid), 'FontSize', 14);
grid on;

% 5. Save
save('output/ga_feature/loocv_results.mat', 'predicted_foot_scores', 'actual_foot_scores', 'r_valid', 'p_valid');
exportgraphics(gcf, 'output/ga_feature/loocv_scatter.png', 'Resolution', 300);