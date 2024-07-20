function L = norm_char_path_length(A)
% Computes the normalized characteristic path length for a directed weighted adjacency matrix
% A should be a square matrix with non-negative elements
% Returns L, the normalized characteristic path length

% Compute the shortest path lengths between all pairs of nodes
D = distance_wei(A);

% Compute the average path length of the graph
L = mean(D(~isinf(D)));

% Compute the theoretical minimum path length of the graph
n = size(A, 1);
L_min = log(n) / log(min(sum(A)));

% Compute the normalized characteristic path length as L divided by L_min
L = L / L_min;
end