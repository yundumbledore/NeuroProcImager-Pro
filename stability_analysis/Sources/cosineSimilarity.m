function sim = cosineSimilarity(v1, v2)
    % Compute the dot product of the two input vectors
    dotProduct = dot(v1, v2);
    
    % Compute the norm (magnitude) of each input vector
    normV1 = norm(v1);
    normV2 = norm(v2);
    
    % Compute the cosine similarity
    sim = dotProduct / (normV1 * normV2);
end