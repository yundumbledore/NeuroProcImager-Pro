function show_dynamic_stability(data_file)
%% Load data
load(['./output/Seizure_' num2str(data_file) '_Jacobian_eig.mat'], 'J_eigenvalues_matrix')

%% Histogram real part of eigenvalues of Jacobi
k = real(J_eigenvalues_matrix);
kk = real(J_eigenvalues_matrix);
kk(k>prctile(k, 99.99, 'all')) = prctile(k, 99.99, 'all');
kk(k<prctile(k, 1, 'all')) = prctile(k, 1, 'all');
clear k J_eigenvalues_matrix

min_val = min(kk,[],'all');
max_val = max(kk,[],'all');
min_val = round(min_val, -1);
max_val = round(max_val, -1);
kk(kk>max_val) = max_val;
kk(kk<min_val) = min_val;

edges = min_val:5:max_val;

for i = 1:size(kk, 2)
    h = histogram(kk(:,i), edges);
    y(:,i) = h.Values;
    close all
end
y = flipud(y);

positive_number = find(edges == max_val) - find(edges == 0);
count = sum(y(1:positive_number,:),1); % find the number of critical mode > 0

%% Adjust count, responsiveness level, gas level for display
startidx = 1;
endidx = numel(count);
countt = count(startidx:endidx);
yy = y(:,startidx:endidx);

countt = movmean(countt, 10);
length = numel(countt);

%% Plot
tiledlayout(2,1);
nexttile
imagesc(log(yy(:,1:length)))
xline([60*10 length-10*10],'w-',{'Seizure onset', 'Seizure offset'},'LabelOrientation','aligned','LineWidth',1, 'FontSize',13, 'LabelVerticalAlignment', 'middle')

inverse_edges = sort(edges,'descend');
idx100 = find(inverse_edges == 100);
idx80 = find(inverse_edges == 80);
idx60 = find(inverse_edges == 60);
idx40 = find(inverse_edges == 40);
idx20 = find(inverse_edges == 20);
idx0 = find(inverse_edges == 0);
idxneg20 = find(inverse_edges == -20);
idxneg40 = find(inverse_edges == -40);
idxneg50 = find(inverse_edges == -50);

yticks([idx100 idx80 idx60 idx40 idx20 idx0 idxneg20 idxneg40])
yticklabels({'100', '80', '60', '40', '20', '0', '-20', '-40'})
ylim([idx100, idx0])

c = colorbar;
c.Title.String = 'log(count)';
ylabel('Real part of eigenvalues')
  
set(gca,'FontSize',12)
set(gca, 'box', 'off')
set(gca,'XTick',[])

nexttile
plot(0.1:1/10:length/10,countt(1:length),'LineWidth',1)
xline([60 length/10-10],'k-',{'Seizure onset', 'Seizure offset'},'LabelOrientation','aligned','LineWidth',1, 'FontSize',13, 'LabelVerticalAlignment', 'middle')
xlim([1, length/10])
xlabel('Time (seconds)')
ylabel('Number of unstable eigenmode')
set(gca,'FontSize',12)
set(gca, 'box', 'off')
 
end