function SWI = small_world_index(A)
% Computes the small-world index for a directed weighted adjacency matrix
% A should be a square matrix with non-negative elements
% Returns SWI, the small-world index

% Compute the clustering coefficient of the graph
C = mean(clustering_coef_wu(A));

% Generate a random graph with the same degree sequence as A
B = randmio_und(A, 100); % Adjust the number of iterations as needed

% Compute the clustering coefficient of the random graph
C_rand = mean(clustering_coef_wu(B));

% Compute the small-world index as the ratio of C to C_rand
SWI = C / C_rand;
end
