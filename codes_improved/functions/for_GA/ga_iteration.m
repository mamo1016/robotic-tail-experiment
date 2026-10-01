function [best_fitness_history, cfg_ga, best_solution] = ga_iteration(eeg_feature, foot_feature, method)
    if nargin < 3
        method = 'spearman';
    end

    cfg_ga = config_ga(eeg_feature, foot_feature);
    fprintf('Starting Evolution (%s mode)...\n', method);
    best_fitness_history = zeros(1, cfg_ga.max_gens);
    
    % --- 4. Main Evolution Loop ---
    for gen = 1:cfg_ga.max_gens
        
        % A. Evaluate Fitness
        fitness_scores = zeros(cfg_ga.pop_size, 1);
        for i = 1:cfg_ga.pop_size
            chromosome = cfg_ga.population(i, :);
            
            % SPLIT the chromosome back into X and Y for testing
            mask_x = chromosome(1 : cfg_ga.num_genes_x);
            mask_y = chromosome(cfg_ga.num_genes_x + 1 : end);
            
            % --- STRICT "JUNK DNA" REMOVAL ---
            % Force genes to be 0 if the underlying feature data is all 0.
            % Since 'remove_unwanted_area' sets columns to 0, we can detect them easily.
            % We create a static mask once (outside the loop would be faster, but this is safe).
            if ~exist('static_dead_genes_x', 'var')
                 static_dead_genes_x = sum(abs(eeg_feature)) == 0;
                 static_dead_genes_y = sum(abs(foot_feature)) == 0;
            end
            
            mask_x(static_dead_genes_x) = 0;
            mask_y(static_dead_genes_y) = 0;
            
            % Updated Fitness Function to accept the split masks
            fitness_scores(i) = calculate_correlation(mask_x, mask_y, eeg_feature, foot_feature, method);
        end
        
        % B. Elitism
        [current_max, best_idx] = max(fitness_scores);
        if current_max > cfg_ga.max_fitness_ever
            cfg_ga.max_fitness_ever = current_max;
            best_solution = cfg_ga.population(best_idx, :);
        end
        best_fitness_history(gen) = cfg_ga.max_fitness_ever;
        
        % C. Create New Generation
        new_population = zeros(size(cfg_ga.population));
        
        % Elitism: Keep the best combined pair
        new_population(1, :) = best_solution;
        
        for i = 2:cfg_ga.pop_size
            % 1. Selection (Now picks a whole row, keeping X and Y linked)
            parent1 = tournament_select(cfg_ga.population, fitness_scores, cfg_ga.tournament_size);
            parent2 = tournament_select(cfg_ga.population, fitness_scores, cfg_ga.tournament_size);
            
            % 2. Crossover (Mixes the super-chromosomes)
            child = crossover(parent1, parent2);
            
            % 3. Mutation
            child = mutate(child, cfg_ga.mutation_rate);

            % --- RE-APPLY JUNK DNA MASK TO CHILD ---
            child_x = child(1 : cfg_ga.num_genes_x);
            child_y = child(cfg_ga.num_genes_x + 1 : end);
            
            if ~exist('static_dead_genes_x', 'var')
                 static_dead_genes_x = sum(abs(eeg_feature)) == 0;
                 static_dead_genes_y = sum(abs(foot_feature)) == 0;
            end
            child_x(static_dead_genes_x) = 0;
            child_y(static_dead_genes_y) = 0;
            
            child = [child_x, child_y];
            
            % temp_child = child(1,1:1200);
            % temp_child(:,eeg_feature(1,:)==0) = 0;
            % child(1,1:1200) = temp_child;
            % temp_child = child(1,1201:end);
            % temp_child(:,foot_feature(1,:)==0) = 0;
            % child(1,1201:end) = temp_child;
    
            new_population(i, :) = child;
        end
        
        cfg_ga.population = new_population;
        
        if mod(gen, 100) == 0
            % fprintf('Gen %d | Best Correlation: %.4f\n', gen, cfg_ga.max_fitness_ever);
        end
    end