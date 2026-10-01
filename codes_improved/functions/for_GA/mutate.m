% 4. Mutation: Randomly flip bits
function child = mutate(child, rate)
    for j = 1:length(child)
        if rand < rate
            child(j) = ~child(j); % Flip 0 to 1, or 1 to 0
        end
    end
end