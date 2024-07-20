function draw_centrality()
    centrality = centrality/max(centrality(:,50:1000),[],'all'); % Normalize centrality
    imagesc(log10(centrality(:,50:1000)))
    c = colorbar;
    c.Title.String = 'Centrality';
    time_avg = mean(centrality,1);
    hold on
    plot3(1:951, ones(1,951), time_avg(50:1000), 'magenta')
    region_avg = mean(centrality(:,50:1000),2);
    [M,I] = maxk(region_avg, 5);
    stem3(951*ones(1,78), 1:78, region_avg, 'red', 'Marker', 'none')
    text(951*ones(1,5),I,M,string(I),'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
    hold off
end