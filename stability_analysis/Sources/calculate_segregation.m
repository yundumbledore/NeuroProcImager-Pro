function [segregation, clustering_coefficient] = calculate_segregation(conn_matrix)
% This function calculates the segregation of a Pearson correlation-based
% functional connectivity matrix using the clustering coefficient.
%
% Inputs:
%   conn_matrix - N x N Pearson correlation-based functional connectivity
%   matrix (N is the number of nodes)
%
% Outputs:
%   segregation - Segregation of the network, defined as the average clustering
%   coefficient across all nodes
%   clustering_coefficient - Clustering coefficient of each node in the network
%
% Author: ChatGPT

% Calculate the clustering coefficient of each node
clustering_coefficient = clustering_coef_wu(conn_matrix);

% Calculate the average clustering coefficient across all nodes
segregation = mean(clustering_coefficient);
end
