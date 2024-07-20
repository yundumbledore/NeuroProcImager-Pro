function DC = degree_centrality(A)
% Computes the degree centrality of a directed weighted adjacency matrix
% A should be a square matrix with non-negative elements
% Returns DC, a vector of degree centralities for each node

% Compute the total degree of each node (in-degree + out-degree)
D = sum(A, 1) + sum(A, 2).';

% Compute the degree centrality of each node by dividing its degree by the
% maximum possible degree (which is n-1 for a directed graph)
n = size(A, 1);
DC = D / (n-1);
end