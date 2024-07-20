function [A] = compute_A(X, Y)
    A = mvregress(X,Y,'algorithm','cwls'); % multivariate regression
end