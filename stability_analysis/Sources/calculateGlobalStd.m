function cosSim = calculateGlobalStd(a, b)
    a = reshape(a, [1, 78*78]);
    b = reshape(b, [1, 78*78]);
    cosSim = abs((a(:).'*b(:))/sqrt(sum(a.^2)*sum(b.^2)));
end