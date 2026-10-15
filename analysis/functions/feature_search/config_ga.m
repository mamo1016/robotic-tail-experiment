function cfg_ga = config_ga(eeg_feature, foot_feature)

cfg_ga =struct();
% --- 1. Mock Data Generation ---
% (Your data loading code remains the same)
% eeg_feature = eeg_feature;
% foot_feature = foot_feature;

% --- 2. Hyperparameters ---
cfg_ga.pop_size = 20; % Increased from 50 (Strategy 2)
cfg_ga.max_gens = 100; % Increased from 500 (Strategy 2)
cfg_ga.mutation_rate = 0.05;
cfg_ga.tournament_size = 3;

% --- 3. Initialization ---
cfg_ga.num_genes_x = size(eeg_feature, 2); 
cfg_ga.num_genes_y = size(foot_feature, 2);
cfg_ga.total_genes = cfg_ga.num_genes_x + cfg_ga.num_genes_y; % Combine lengths

% Create ONE giant population containing both X and Y genes
cfg_ga.population = randi([0, 1], cfg_ga.pop_size, cfg_ga.total_genes);

% x = cfg_ga.population(:,1:1200);
% x(:,eeg_feature(1,:)==0) = 0;
% cfg_ga.population(:,1:1200) = x;
% 
% x = cfg_ga.population(:,1201:end);
% x(:,foot_feature(1,:)==0) = 0;
% cfg_ga.population(:,1201:end) = x;

cfg_ga.best_fitness_history = zeros(cfg_ga.max_gens, 1);
cfg_ga.best_solution = [];
cfg_ga.max_fitness_ever = -inf;