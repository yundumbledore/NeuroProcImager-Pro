function avg_local_efficiency = directed_weighted_local_efficiency(A)
% A is the adjacency matrix of a directed weighted graph
n = size(A, 1);
local_efficiency = zeros(n, 1);
for i = 1:n
    % Get the neighbors of node i
    neighbors = find(A(i,:) > 0);
    % Get the subgraph induced by the neighbors of node i
    A_subgraph = A(neighbors, neighbors);
    % Calculate the local efficiency of node i
    if length(neighbors) > 1
        local_efficiency(i) = efficiency_wei(A_subgraph);
    end
end
% Calculate the average local efficiency
avg_local_efficiency = sum(local_efficiency)/n;
end
