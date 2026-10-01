%% Custom Genetic Algorithm for Feature Selection (No Toolbox Required)
clear; clc; close all;
addpath('functions\for_obtainERD'); % Add the 'functions' folder to MATLAB's search path
addpath('functions\for_ersp'); % Add the 'functions' folder to MATLAB's search path
addpath('functions\for_GA'); % Add the 'functions' folder to MATLAB's search path
ersp_epoch_bootstraped_store = ersp_epoch_bootstraped_store_get();


%%
%visualise ersp and save ersps
config=configure_parameters(); % configure parameters
fields = fieldnames(ersp_epoch_bootstraped_store);           
current_name = fields{1};           % Get the name (string)  
                
[time_range, freq_range] = range_compute(config);
fig_ersp = figure('Color', 'w'); % 'Name' sets the window title
figure(fig_ersp)
tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');   
visualisation_ersp(ersp_epoch_bootstraped_store.(current_name).unexpected_ersp_ave - ersp_epoch_bootstraped_store.(current_name).expected_ersp_ave, fig_ersp, time_range, freq_range); title("subA unexpected - expected ERSP");
visualisation_ersp(ersp_epoch_bootstraped_store.(current_name).right_after_ersp_ave - ersp_epoch_bootstraped_store.(current_name).right_before_ersp_ave, fig_ersp, time_range, freq_range); title("subA right after - right before ERSP");
visualisation_ersp(ersp_epoch_bootstraped_store.(current_name).expected_duplicated_ave - ersp_epoch_bootstraped_store.(current_name).expected_ersp_ave, fig_ersp, time_range, freq_range); title("subA expected - expected ERSP");

fig = gcf; % Get current figure handle
exportgraphics(fig, "output\all_in_one_pic\sub_A_avergae_multiple_bootstrap.png", 'Resolution', 300);

filename = "output\ersp_each_sub_ave\ersp_epoch_bootstraped_store.mat";        
save(filename, 'ersp_epoch_bootstraped_store', '-v7.3');

%%
% make ersp to tiles and save
fig_ersp_tile = figure('Color', 'w'); % 'Name' sets the window title
figure(fig_ersp_tile)
% tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');   

for struct_num = 1:length(fields)
    current_name = fields{struct_num};           % Get the name (string)  
    sub_tile.(current_name).expected = ersp_tile(ersp_epoch_bootstraped_store.(current_name).expected_ersp_ave,time_range, freq_range);
    sub_tile.(current_name).unexpected = ersp_tile(ersp_epoch_bootstraped_store.(current_name).unexpected_ersp_ave,time_range, freq_range);
    sub_tile.(current_name).right_before = ersp_tile(ersp_epoch_bootstraped_store.(current_name).right_before_ersp_ave,time_range, freq_range);
    sub_tile.(current_name).right_after = ersp_tile(ersp_epoch_bootstraped_store.(current_name).right_after_ersp_ave,time_range, freq_range);
    sub_tile.(current_name).duplicated = ersp_tile(ersp_epoch_bootstraped_store.(current_name).expected_duplicated_ave,time_range, freq_range);
    % visualisation_ersp(sub_tile.(current_name).unexpected.vector - sub_tile.(current_name).expected.vector, fig_ersp_tile, sub_tile.(current_name).expected.t, sub_tile.(current_name).expected.f); title("subA unexpected - expected ERSP");

end

filename = "output\ersp_each_sub_ave\sub_tile.mat";        
save(filename, 'sub_tile', '-v7.3');


fig_ersp_tile = figure('Color', 'w'); % 'Name' sets the window title
figure(fig_ersp_tile)
tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');   
visualisation_ersp(sub_tile.(current_name).unexpected.vector - sub_tile.(current_name).expected.vector, fig_ersp_tile, sub_tile.(current_name).expected.t, sub_tile.(current_name).expected.f); title("subA unexpected - expected ERSP");
visualisation_ersp(sub_tile.(current_name).right_after.vector - sub_tile.(current_name).right_before.vector, fig_ersp_tile, sub_tile.(current_name).expected.t, sub_tile.(current_name).expected.f); title("subA right after - right before ERSP");
visualisation_ersp(sub_tile.(current_name).duplicated.vector- sub_tile.(current_name).expected.vector, fig_ersp_tile, sub_tile.(current_name).expected.t, sub_tile.(current_name).expected.f); title("subA expected - expected ERSP");

%%
foot_cleaned_session_combined = load("output\foot\foot_cleaned_session_combined.mat").foot_data_boot_combined;
fields = fieldnames(foot_cleaned_session_combined);           
foot_data_store = struct();    
for struct_num = 1:length(fields)
    current_name = fields{struct_num};
    fields_types = fieldnames(foot_cleaned_session_combined.(current_name)); 
    type_name = fields_types{struct_num};
    for type_num = 1:length(type_name)
        data = foot_cleaned_session_combined.(current_name).(type_name);
        %here think parameters to suggest emobodiment
        foot_data_store.(current_name).(type_name).normal_mean = mean(data,1); 
        foot_data_store.(current_name).(type_name).normal_mean = mean(data,1); 

    end
end

%%

data = load("output\all_in_one_pic\avergae_multiple_bootstrap.mat").ave_plot;
% --- 1. Mock Data Generation (Replace this with your real X and Y) ---
% X: 25 Participants x 40 "Tiles" (Time-Freq bins)
% Y: 25 Participants x 1 CoP Metric
num_subs = 25;
num_features = 40; % e.g., 8 time-bins * 5 freq-bins
X_features = randn(num_subs, num_features); 
Y_cop = randn(num_subs, 1);

% --- 2. Hyperparameters ---
pop_size = 50;          % How many "chromosomes" in the pool
max_gens = 100;         % How many times to evolve
mutation_rate = 0.05;   % 5% chance a gene flips (maintains diversity)
tournament_size = 3;    % For selection (higher = more pressure to pick best)

% --- 3. Initialization ---
% Create random binary masks (0 or 1)
population = randi([0, 1], pop_size, num_features);
best_fitness_history = zeros(max_gens, 1);
best_solution = [];
max_fitness_ever = -inf;

fprintf('Starting Evolution...\n');

% --- 4. Main Evolution Loop ---
for gen = 1:max_gens
    
    % A. Evaluate Fitness for entire population
    fitness_scores = zeros(pop_size, 1);
    for i = 1:pop_size
        mask = population(i, :);
        fitness_scores(i) = calculate_correlation(mask, X_features, Y_cop);
    end
    
    % B. Elitism: Track and Keep the Best
    [current_max, best_idx] = max(fitness_scores);
    if current_max > max_fitness_ever
        max_fitness_ever = current_max;
        best_solution = population(best_idx, :);
    end
    best_fitness_history(gen) = max_fitness_ever;
    
    % C. Create New Generation
    new_population = zeros(size(population));
    
    % Elitism: Always carry over the single best parent unchanged
    new_population(1, :) = best_solution;
    
    % Fill the rest of the new population
    for i = 2:pop_size
        % 1. Selection (Tournament)
        parent1 = tournament_select(population, fitness_scores, tournament_size);
        parent2 = tournament_select(population, fitness_scores, tournament_size);
        
        % 2. Crossover (Single Point)
        child = crossover(parent1, parent2);
        
        % 3. Mutation
        child = mutate(child, mutation_rate);
        
        new_population(i, :) = child;
    end
    
    % Update population for next round
    population = new_population;
    
    % Display progress every 10 generations
    if mod(gen, 10) == 0
        fprintf('Gen %d | Best Correlation: %.4f\n', gen, max_fitness_ever);
    end
end

% --- 5. Visualization & Results ---
figure; plot(best_fitness_history, 'LineWidth', 2);
xlabel('Generation'); ylabel('Correlation (Fitness)');
title('Evolution of Best Correlation');
grid on;

fprintf('\nFinal Best Correlation: %.4f\n', max_fitness_ever);
fprintf('Selected Features (Indices): \n');
disp(find(best_solution == 1));

%% --- HELPER FUNCTIONS (The "Engine" Parts) ---

% 1. Fitness Function: Calculates Correlation
function r = calculate_correlation(mask, X, Y)
    if sum(mask) == 0
        r = 0; return; % Avoid crash if mask is empty
    end
    % Filter X using the mask
    selected_data = X(:, mask == 1);
    % Average the selected features (mean power)
    mean_val = mean(selected_data, 2);
    % Calculate correlation
    corr_mat = corrcoef(mean_val, Y);
    r = abs(corr_mat(1,2)); % Maximise absolute correlation
end

% 2. Tournament Selection: Pick winners to be parents
function winner = tournament_select(pop, fitness, k)
    pop_size = size(pop, 1);
    % Pick 'k' random indices
    contestants_idx = randi(pop_size, k, 1);
    % Find which one has the highest fitness
    [~, best_loc] = max(fitness(contestants_idx));
    % Return that chromosome
    winner = pop(contestants_idx(best_loc), :);
end

% 3. Crossover: Mix two parents
function child = crossover(p1, p2)
    num_genes = length(p1);
    % Pick a random split point
    split_point = randi([1, num_genes-1]);
    % Combine: Head of P1 + Tail of P2
    child = [p1(1:split_point), p2(split_point+1:end)];
end

% 4. Mutation: Randomly flip bits
function child = mutate(child, rate)
    for j = 1:length(child)
        if rand < rate
            child(j) = ~child(j); % Flip 0 to 1, or 1 to 0
        end
    end
end