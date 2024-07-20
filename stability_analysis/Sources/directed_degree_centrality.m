function deg_centrality = directed_degree_centrality(adj_matrix)
% Compute the degree centrality for a directed weighted graph adjacency matrix
%
% Inputs:
% adj_matrix: n x n directed weighted graph adjacency matrix
%
% Outputs:
% deg_centrality: n x 1 vector of degree centralities, where the i-th entry
% represents the degree centrality of the i-th node

% Compute the total degree (in-degree + out-degree) of each node
total_degree = sum(adj_matrix, 1);

% Compute the maximum total degree in the graph
max_degree = max(total_degree);

% Normalize the total degree by the maximum degree to obtain the degree centrality
deg_centrality = total_degree;

end
