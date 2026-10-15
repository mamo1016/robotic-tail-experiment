function winner = tournament_select(pop, fitness, k)
    cfg_ga.pop_size = size(pop, 1);
    % Pick 'k' random indices
    contestants_idx = randi(cfg_ga.pop_size, k, 1);
    % Find which one has the highest fitness
    [~, best_loc] = max(fitness(contestants_idx));
    % Return that chromosome
    winner = pop(contestants_idx(best_loc), :);
end