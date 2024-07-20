function show_dynamic_chaos(data_file)
%% Load data
load(['./output/Seizure_' num2str(data_file) '_Lyapunov_exponents.mat'])

%% Plot
tiledlayout(2,1);
nexttile
imagesc(time_varying_lyapunov_exponents)
xline([60*10 size(time_varying_lyapunov_exponents,2)-10*10],'k-',{'Seizure onset', 'Seizure offset'},'LabelOrientation','aligned','LineWidth',1, 'FontSize',13, 'LabelVerticalAlignment', 'middle')
c = colorbar;
c.Title.String = 'Lyapunov exponent';
ylabel('Model state')
set(gca,'FontSize',12)
set(gca, 'box', 'off')
set(gca,'XTick',[])

nexttile
l_max = max(time_varying_lyapunov_exponents,[],1);
plot(1/10:1/10:size(time_varying_lyapunov_exponents,2)/10,l_max,'LineWidth',1,'Color','black')
xline([60 size(time_varying_lyapunov_exponents,2)/10-10],'k-',{'Seizure onset', 'Seizure offset'},'LabelOrientation','aligned','LineWidth',1, 'FontSize',13, 'LabelVerticalAlignment', 'middle')
xlim([1, size(time_varying_lyapunov_exponents,2)/10])
xlabel('Time (seconds)')
ylabel('Maximal lyapunov exponent')
set(gca,'FontSize',12)
set(gca, 'box', 'off')
end