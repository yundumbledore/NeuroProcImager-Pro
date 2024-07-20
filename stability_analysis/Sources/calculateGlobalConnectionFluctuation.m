function [s] = calculateGlobalConnectionFluctuation(a, window_size)
    s = movstd(a, window_size);
end