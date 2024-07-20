function segregation = calculateSegregation(adjacencyMatrix)
% Calculates the segregation of a directed weighted graph given its
% adjacency matrix using the formula presented in the paper "A Framework for
% Comparing Structural and Functional Network Segregation" by Bassett and
% Bullmore (2006).

% Step 1: Compute the binary adjacency matrix
binaryMatrix = double(adjacencyMatrix > 0);

% Step 2: Compute the in-degree and out-degree of each node
inDegree = sum(binaryMatrix, 1)';
outDegree = sum(binaryMatrix, 2);

% Step 3: Compute the total weight of incoming and outgoing edges for each node
inWeight = sum(adjacencyMatrix, 1)';
outWeight = sum(adjacencyMatrix, 2);

% Step 4: Compute the mean in-degree and out-degree of the graph
meanInDegree = mean(inDegree);
meanOutDegree = mean(outDegree);

% Step 5: Compute the mean weight of incoming and outgoing edges for the graph
meanInWeight = mean(inWeight);
meanOutWeight = mean(outWeight);

% Step 6: Compute the segregation coefficient of each node
segregation = zeros(size(adjacencyMatrix, 1), 1);
for i = 1:size(adjacencyMatrix, 1)
    % Compute the total weight of incoming and outgoing edges for nodes
    % with the same in-degree as node i
    inWeightSameInDegree = sum(adjacencyMatrix(inDegree == inDegree(i), i));
    outWeightSameInDegree = sum(adjacencyMatrix(i, outDegree == outDegree(i)));
    
    % Compute the total weight of incoming and outgoing edges for nodes
    % with the same out-degree as node i
    inWeightSameOutDegree = sum(adjacencyMatrix(outDegree == outDegree(i), i));
    outWeightSameOutDegree = sum(adjacencyMatrix(i, inDegree == inDegree(i)));
    
    % Compute the expected total weight of incoming and outgoing edges for
    % nodes with the same in-degree and out-degree as node i
    expectedInWeight = (inDegree(i) / meanInDegree) * meanInWeight;
    expectedOutWeight = (outDegree(i) / meanOutDegree) * meanOutWeight;
    
    % Compute the segregation coefficient of node i
    numerator = inWeightSameInDegree - expectedInWeight + outWeightSameOutDegree - expectedOutWeight;
    denominator = inWeightSameInDegree + expectedInWeight + outWeightSameOutDegree + expectedOutWeight - adjacencyMatrix(i, i);
    segregation(i) = numerator / denominator;
end

% Step 7: Compute the segregation as the average of the segregation coefficients
segregation = mean(segregation);

end
