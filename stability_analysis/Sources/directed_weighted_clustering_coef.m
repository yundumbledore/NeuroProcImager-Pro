function c_norm = directed_weighted_clustering_coef(A)
% A normalized clustering coefficient of 1 indicates that every node in 
% the graph is connected to every other node in its neighborhood, 
% which implies a very high degree of local connectivity. 
% On the other hand, a normalized clustering coefficient of 0 indicates 
% that there are no triangles in the graph, and therefore no local connectivity.

% A is the adjacency matrix of a directed weighted graph
n = size(A, 1);
c = zeros(n, 1);
for i = 1:n
    % Get the neighbors of node i
    neighbors = find(A(i,:) > 0);
    % Count the number of connections between neighbors of i
    num_edges = 0;
    for j = 1:length(neighbors)
        for k = j+1:length(neighbors)
            if A(neighbors(j),neighbors(k)) > 0
                num_edges = num_edges + A(neighbors(j),neighbors(k));
            end
        end
    end
    % Calculate the clustering coefficient of node i
    if length(neighbors) > 1
        c(i) = 2*num_edges/(length(neighbors)*(length(neighbors)-1));
    end
end
% Calculate the normalized clustering coefficient
c_norm = sum(c)/n;
end
