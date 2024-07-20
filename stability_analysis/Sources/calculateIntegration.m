function integration = calculateIntegration(adj_matrix)
% Here is a Matlab function that can calculate the integration for a directed 
% weighted graph with an adjacency matrix that can contain elements greater than 1 or even negative values

% In this function, adj_matrix is the adjacency matrix of the graph, which can 
% contain elements greater than 1 or negative values. The function calculates 
% the in-degrees and out-degrees of each node in the graph, and then uses these values 
% to calculate the integration of the graph using the formula I = ΣiΣj(aij*di+dj-)/(Σi(di+^2)+Σj(dj-^2)-ΣiΣj(aijdi+*dj-)), 
% where di+ and dj- are the in-degree and out-degree of node i and j, respectively.

% This formula measures the degree to which a network is functionally integrated, 
% meaning the degree to which information can flow between different regions of 
% the network. It is based on the idea that nodes that are highly connected with 
% nodes that have different input/output patterns (i.e. nodes with high out-degree 
% connected to nodes with high in-degree) are important for integration. 
% The numerator of the formula calculates the total weight of all the edges 
% in the network, weighted by the in-degree and out-degree of the nodes they connect. 
% The denominator is a normalization factor that takes into account the total number 
% of edges in the network and the degree distribution of the nodes.

n = size(adj_matrix, 1);  % Number of nodes in the graph

in_degrees = sum(adj_matrix, 1);  % Calculate in-degrees of nodes
out_degrees = sum(adj_matrix, 2)';  % Calculate out-degrees of nodes

sum_in = sum(in_degrees .^ 2);
sum_out = sum(out_degrees .^ 2);

sum_in_out = 0;
for i = 1:n
    for j = 1:n
        sum_in_out = sum_in_out + (adj_matrix(i,j) * in_degrees(i) * out_degrees(j)) / (in_degrees(i) + out_degrees(j) - adj_matrix(i,j));
    end
end

integration = sum_in_out / (sum_in + sum_out - sum_in_out);

end

