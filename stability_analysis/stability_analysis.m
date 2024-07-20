function stability_analysis(data_file)
    warning('off')
    
    %% Define model parameters
    v0 = 6;
    varsigma = 3.0339;
    scale = 50;
    ex_tau = 0.010;
    in_tau = 0.020;
    d_tau = 1/33;
    
    fs = 400;
    save_fs = 150; % sampling rate estimates were saved
    
    %% Define RUN Settings
    ntasks = feature('numcores'); % num of CPUs
    
    window_size_inSeconds = 0.1;
    window_length = window_size_inSeconds*save_fs;
    
    data_dir = './output';
    output_dir = data_dir;
    
    %% Run
    i_file = data_file;
    moduler(data_dir, i_file, scale, fs, window_length, output_dir, v0, varsigma, ex_tau, in_tau, d_tau, ntasks)
end

%% Subfunctions
function moduler(data_path, i_file, scale, fs, window_length, output_dir, v0, varsigma, ex_tau, in_tau, d_tau, ntasks)
    %% Load everything before run
    load([data_path '/Seizure_' num2str(i_file) '_estimates.mat'], 'xi_hat_list') % load data
    xi_hat_list(isnan(xi_hat_list))=0; % replace nan with 0
    
    n_channels = size(xi_hat_list,3);
    n_timesteps = size(xi_hat_list,2);
    windows_num = floor(n_timesteps/window_length)-1;

    xi_hat_list = double(xi_hat_list);
    v_pyr = (squeeze(xi_hat_list(1,:,:))' + squeeze(xi_hat_list(7,:,:))' + squeeze(xi_hat_list(9,:,:))')/scale;
    v_ip = squeeze(xi_hat_list(1,:,:))'/scale;
    z_ip = squeeze(xi_hat_list(2,:,:))';
    v_pi = squeeze(xi_hat_list(3,:,:))'/scale;
    z_pi = squeeze(xi_hat_list(4,:,:))';
    v_pe = squeeze(xi_hat_list(5,:,:))'/scale;
    z_pe = squeeze(xi_hat_list(6,:,:))';
    v_ep = squeeze(xi_hat_list(7,:,:))'/scale;
    z_ep = squeeze(xi_hat_list(8,:,:))';
    input = squeeze(xi_hat_list(9,:,:))'/scale;
    alpha_ip = fs*squeeze(xi_hat_list(10,:,:))';
    alpha_pi = fs*squeeze(xi_hat_list(11,:,:))';
    alpha_pe = fs*squeeze(xi_hat_list(12,:,:))';
    alpha_ep = fs*squeeze(xi_hat_list(13,:,:))';
    clear xi_hat_list

    load([output_dir '/Seizure_' num2str(i_file) '_connectivity_matrix.mat'], 'A_sequence')

    %% Find an equilibrium point of the system
    H = [1; 0; 0; 0; 0; 0; 1; 0; 1; 0];
    H_matrix = zeros(n_channels*10, n_channels);
    for iCh = 1:n_channels
        H_matrix((iCh-1)*10+1:(iCh-1)*10+10, iCh) = H;
    end
    
    options = optimoptions('fsolve','MaxIterations',1e5,'MaxFunctionEvaluations',1e5,'UseParallel',false, 'FunctionTolerance',1e-5,'OptimalityTolerance',1e-5, 'Display', 'off'); % define fsolve settings before run
    
    fixpoint = nan(n_channels*10, windows_num);
    poolobj = parpool(ntasks);
    parfor (i = 1:windows_num, ntasks)
%         try
            v_ip_avg = scale*nanmean(v_ip(:,(i-1)*window_length+1:i*window_length),2);
            z_ip_avg = nanmean(z_ip(:,(i-1)*window_length+1:i*window_length),2);
            v_pi_avg = scale*nanmean(v_pi(:,(i-1)*window_length+1:i*window_length),2);
            z_pi_avg = nanmean(z_pi(:,(i-1)*window_length+1:i*window_length),2);
            v_pe_avg = scale*nanmean(v_pe(:,(i-1)*window_length+1:i*window_length),2);
            z_pe_avg = nanmean(z_pe(:,(i-1)*window_length+1:i*window_length),2);
            v_ep_avg = scale*nanmean(v_ep(:,(i-1)*window_length+1:i*window_length),2);
            z_ep_avg = nanmean(z_ep(:,(i-1)*window_length+1:i*window_length),2);
            mu_avg = scale*nanmean(input(:,(i-1)*window_length+1:i*window_length),2);
            alpha_ip_avg = nanmean(alpha_ip(:,(i-1)*window_length+1:i*window_length),2);
            alpha_pi_avg = nanmean(alpha_pi(:,(i-1)*window_length+1:i*window_length),2);
            alpha_pe_avg = nanmean(alpha_pe(:,(i-1)*window_length+1:i*window_length),2);
            alpha_ep_avg = nanmean(alpha_ep(:,(i-1)*window_length+1:i*window_length),2);   

            x0 = zeros(1, n_channels*10);
            for iCh = 1:n_channels
                start_idx = (iCh-1)*10+1;
                end_idx = 10*iCh;
                x0(1, start_idx:end_idx) = [v_ip_avg(iCh) z_ip_avg(iCh) v_pi_avg(iCh) z_pi_avg(iCh) v_pe_avg(iCh) z_pe_avg(iCh) v_ep_avg(iCh) z_ep_avg(iCh) mu_avg(iCh) 0];
            end
            w_matrix = A_sequence(:,:,i);

            eq = find_equilibrium_network(x0,options,n_channels,alpha_ip_avg,alpha_pi_avg,alpha_pe_avg,alpha_ep_avg,w_matrix,ex_tau,in_tau,H_matrix,d_tau,v0,varsigma,scale);
            fixpoint(:,i) = eq';
    end
    delete(poolobj) 
    save([output_dir '/Seizure_' num2str(i_file) '_fix_points.mat'], 'fixpoint', '-v7.3')

    %% Calculate the Jacobian matrix of the system at the equilibrium point
    J_tensor = nan(n_channels*10, n_channels*10, windows_num);
    for i = 1:windows_num
        try
            eq = fixpoint(:,i);
            w_matrix = A_sequence(:,:,i);
            alpha_ip_avg = mean(alpha_ip(:,(i-1)*window_length+1:i*window_length),2);
            alpha_pi_avg = mean(alpha_pi(:,(i-1)*window_length+1:i*window_length),2);
            alpha_pe_avg = mean(alpha_pe(:,(i-1)*window_length+1:i*window_length),2);
            alpha_ep_avg = mean(alpha_ep(:,(i-1)*window_length+1:i*window_length),2);  
            J_tensor(:,:,i) = construct_Jacobian(n_channels, eq, w_matrix, alpha_ip_avg, alpha_pi_avg, alpha_pe_avg, alpha_ep_avg, ex_tau, in_tau, d_tau, varsigma, v0, scale);
        catch
            continue
        end
    end
    save([output_dir '/Seizure_' num2str(i_file) '_Jacobian.mat'], 'J_tensor', '-v7.3')
    
    %% Calculate the eigenvalues of the Jacobian matrix
    J_eigenvalues_matrix = nan(n_channels*10, windows_num);
    for i = 1:windows_num
        try
            J_eigenvalues_matrix(:,i) = eig(J_tensor(:,:,i));
        catch
            continue
        end
    end
    save([output_dir '/Seizure_' num2str(i_file) '_Jacobian_eig.mat'], 'J_eigenvalues_matrix', '-v7.3')
end