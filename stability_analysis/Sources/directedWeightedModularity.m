function Q = directedWeightedModularity(A)
% A: directed weighted adjacency matrix

n = size(A, 1); % number of nodes
m = sum(A(:)); % total weight of edges

% Compute community assignments using Louvain algorithm
[~, communities] = community_louvain(A, 1000);

% Compute modularity matrix
ki = sum(A, 2);
kj = sum(A, 1);
pij = ki*kj/m;
B = A - pij;

% Compute modularity
Q = 0;
for i = 1:max(communities)
    nodes = find(communities == i);
    ei = sum(B(nodes, :), 1);
    Q = Q + sum(ei(communities == i)) - sum(ki(nodes))^2/m;
end
Q = Q/m;
