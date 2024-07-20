function show_estimates_avg_errorband(xi_hat_list, output_dir, sub_ind, i_file, data)
    if size(xi_hat_list,1) == 13
        y_pred = squeeze(xi_hat_list(1,:,:) + xi_hat_list(7,:,:) + xi_hat_list(9,:,:))/50;

        subplot(4,4,1);
        y = nanmean(data, [2])';
        stderr = nanstd(data,0,2)'/sqrt(16);
        plot_mean_errorband_y(y,stderr);
        title('EEG signal')

    
        subplot(4,4,2);
        y = nanmean(y_pred,[2])';
        stderr = nanstd(y_pred,0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Prediction')


        subplot(4,4,5);
        y = nanmean(xi_hat_list(9,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(9,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('External input estimates')
        set(gca,'YTickLabel',[]);
        
        subplot(4,4,9);
        y = nanmean(xi_hat_list(10,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(10,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Inhibitory to pyramidal connectivity estimates')
        set(gca,'YTickLabel',[]);
    
        subplot(4,4,10);
        y = nanmean(xi_hat_list(11,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(11,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Pyramidal to inhibitory connectivity estimates')
        set(gca,'YTickLabel',[]);
        
        subplot(4,4,11);
        y = nanmean(xi_hat_list(12,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(12,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Pyramidal to excitatory connectivity estimates')
        set(gca,'YTickLabel',[]);
    
        subplot(4,4,12);
        y = nanmean(xi_hat_list(13,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(13,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Excitatory to pyramidal connectivity estimates')
        set(gca,'YTickLabel',[]);

        subplot(4,4,13);
        y = nanmean(xi_hat_list(1,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(1,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Inhibitory to pyramidal potential estimates')
        set(gca,'YTickLabel',[]);
    
        subplot(4,4,14);
        y = nanmean(xi_hat_list(3,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(3,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Pyramidal to inhibitory potential estimates')
        set(gca,'YTickLabel',[]);
        
        subplot(4,4,15);
        y = nanmean(xi_hat_list(5,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(5,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Pyramidal to excitatory potential estimates')
        set(gca,'YTickLabel',[]);
    
        subplot(4,4,16);
        y = nanmean(xi_hat_list(7,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(7,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Excitatory to pyramidal potential estimates')
        set(gca,'YTickLabel',[]);
    else
        y_pred = squeeze(xi_hat_list(1,:,:) + xi_hat_list(7,:,:) + xi_hat_list(9,:,:))/50;

        subplot(4,4,1);
        y = nanmean(data, [2])';
        stderr = nanstd(data,0,2)'/sqrt(16);
        plot_mean_errorband_y(y,stderr);
        title('EEG signal')

    
        subplot(4,4,2);
        y = nanmean(y_pred,[2])';
        stderr = nanstd(y_pred,0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Prediction')


        subplot(4,4,5);
        y = nanmean(xi_hat_list(9,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(9,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('External input estimates')
        set(gca,'YTickLabel',[]);

        subplot(4,4,6);
        y = nanmean(xi_hat_list(10,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(10,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('DC offset estimates')
        set(gca,'YTickLabel',[]);
        
        subplot(4,4,9);
        y = nanmean(xi_hat_list(11,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(11,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Inhibitory to pyramidal connectivity estimates')
        set(gca,'YTickLabel',[]);
    
        subplot(4,4,10);
        y = nanmean(xi_hat_list(12,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(12,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Pyramidal to inhibitory connectivity estimates')
        set(gca,'YTickLabel',[]);
        
        subplot(4,4,11);
        y = nanmean(xi_hat_list(13,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(13,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Pyramidal to excitatory connectivity estimates')
        set(gca,'YTickLabel',[]);
    
        subplot(4,4,12);
        y = nanmean(xi_hat_list(14,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(14,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Excitatory to pyramidal connectivity estimates')
        set(gca,'YTickLabel',[]);

        subplot(4,4,13);
        y = nanmean(xi_hat_list(1,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(1,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Inhibitory to pyramidal potential estimates')
        set(gca,'YTickLabel',[]);
    
        subplot(4,4,14);
        y = nanmean(xi_hat_list(3,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(3,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Pyramidal to inhibitory potential estimates')
        set(gca,'YTickLabel',[]);
        
        subplot(4,4,15);
        y = nanmean(xi_hat_list(5,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(5,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Pyramidal to excitatory potential estimates')
        set(gca,'YTickLabel',[]);
    
        subplot(4,4,16);
        y = nanmean(xi_hat_list(7,:,:),[3]);
        stderr = nanstd(squeeze(xi_hat_list(7,:,:)),0,2)'/sqrt(16);
        plot_mean_errorband(y,stderr);
        title('Excitatory to pyramidal potential estimates')
        set(gca,'YTickLabel',[]);
%     set(gcf,'position',[10,10,1400,1600])
%     print(gcf, [output_dir 'S' num2str(sub_ind) '/Seizure_' num2str(i_file) '_estimates.png'],'-dpng','-r300');
%     close all
    end
end

function plot_mean_errorband(y,stderr)
    x = (1:numel(y))/150;
    curve1 = y + 1.96*stderr;
    curve2 = y - 1.96*stderr;
    
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    
%     fill(x2, inBetween, 'g','LineStyle','none');
%     hold on;
    plot(x, y);
    xlabel('Time (Seconds)')
%     hold off
end

function plot_mean_errorband_y(y,stderr)
    x = (1:numel(y))/400;
    curve1 = y + 1.96*stderr;
    curve2 = y - 1.96*stderr;
    
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    
%     fill(x2, inBetween, 'g','LineStyle','none');
%     hold on;
    plot(x, y);
    xlabel('Time (Seconds)')
%     hold off
end