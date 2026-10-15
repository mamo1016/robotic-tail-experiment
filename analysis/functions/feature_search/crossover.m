% 3. Crossover: Mix two parents
function child = crossover(p1, p2)
    num_genes = length(p1);
    % Pick a random split point
    split_point = randi([1, num_genes-1]);
    % Combine: Head of P1 + Tail of P2
    child = [p1(1:split_point), p2(split_point+1:end)];
end