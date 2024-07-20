function [centrality] = eigenvector_centrality(W)
% Input:
%   W: weighted adjacency matrix of a directed graph
% Output:
%   centrality: eigenvector centrality of each node in the graph

% normalize the adjacency matrix
D = diag(sum(W,2)); % degree matrix
P = D\W; % normalized adjacency matrix

% compute the dominant eigenvector
[V, ~] = eigs(P',1); % compute the dominant eigenvector of P'
centrality = abs(V);%/sum(abs(V)); % normalize the eigenvector

end