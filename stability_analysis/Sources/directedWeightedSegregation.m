function [S, Q] = directedWeightedSegregation(A)
% This MATLAB function calculates the segregation and weighted clustering 
% coefficient of a directed, weighted adjacency matrix A, which represents 
% the connections between nodes in a network. The function first computes 
% the incoming and outgoing weights for each node, as well as the node strengths. 
% It then calculates the weighted clustering coefficient using a nested loop over all node pairs. 
% Next, it computes the segregation coefficient Q, which is the average clustering 
% coefficient across all nodes. Finally, the function computes the within-group 
% fraction of edges S, which quantifies the extent to which edges between nodes 
% are concentrated within specific groups, and normalizes it by the maximum possible value of S.

% A: directed weighted adjacency matrix

n = size(A, 1); % number of nodes
m = sum(A(:)); % total weight of edges

% Compute incoming and outgoing weights for each node
w_in = sum(A, 1)';
w_out = sum(A, 2);

% Compute node strengths
s_in = sum(A, 1);
s_out = sum(A, 2);

% Compute weighted clustering coefficient
C = zeros(n, 1);
for i = 1:n
    w_ii_in = 0;
    w_ii_out = 0;
    for j = 1:n
        if A(i,j) > 0
            w_ij_in = A(j,i);
            w_ii_in = w_ii_in + w_ij_in;
            for k = 1:n
                if A(j,k) > 0
                    w_ik_out = A(j,k);
                    w_ii_out = w_ii_out + w_ij_in*w_ik_out;
                end
            end
        end
    end
    C(i) = w_ii_in*w_ii_out/(s_in(i)*s_out(i));
end

% Compute segregation coefficient
Q = sum(C)/n;

% Compute maximum possible value of S
S_max = 0;
for i = 1:n
    for j = 1:n
        if A(i,j) > 0
            S_max = S_max + (w_in(i)*w_in(j) + w_out(i)*w_out(j));
        end
    end
end
S_max = S_max/(2*m)^2;

% Compute normalized weighted within-group fraction of edges
S = 0;
for i = 1:n
    for j = 1:n
        if A(i,j) > 0
            S = S + A(i,j)*(w_in(i)*w_in(j) + w_out(i)*w_out(j))/(2*m);
        end
    end
end
S = S/(2*m)^2/S_max; % normalize by maximum possible value of S

end

% function [S, Q] = directedWeightedSegregation(A)
% % A: directed weighted adjacency matrix
% 
% n = size(A, 1); % number of nodes
% m = sum(A(:)); % total weight of edges
% 
% % Compute incoming and outgoing weights for each node
% w_in = sum(A, 1)';
% w_out = sum(A, 2);
% 
% % Compute node strengths
% s_in = sum(A, 1);
% s_out = sum(A, 2);
% 
% % Compute weighted clustering coefficient
% C = zeros(n, 1);
% for i = 1:n
%     w_ii_in = 0;
%     w_ii_out = 0;
%     for j = 1:n
%         if A(i,j) > 0
%             w_ij_in = A(j,i);
%             w_ii_in = w_ii_in + w_ij_in;
%             for k = 1:n
%                 if A(j,k) > 0
%                     w_ik_out = A(j,k);
%                     w_ii_out = w_ii_out + w_ij_in*w_ik_out;
%                 end
%             end
%         end
%     end
%     C(i) = w_ii_in*w_ii_out/(s_in(i)*s_out(i));
% end
% 
% % Compute segregation coefficient
% Q = sum(C)/n;
% 
% % Compute weighted within-group fraction of edges
% S = 0;
% for i = 1:n
%     for j = 1:n
%         if A(i,j) > 0
%             S = S + A(i,j)*(w_in(i)*w_in(j) + w_out(i)*w_out(j))/(2*m);
%         end
%     end
% end
% S = S/m;
% end
