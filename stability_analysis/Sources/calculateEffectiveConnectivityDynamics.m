function ecd = calculateEffectiveConnectivityDynamics(A_sequence, time_window)
    kk = movmean(A_sequence, time_window);
    for i = 1:size(A_sequence,3)
        for j = 1:size(A_sequence,3)
            a = kk(:,:,i);
            b = kk(:,:,j);
            ecd(i,j) = corr2(a(:),b(:));
        end
    end
end