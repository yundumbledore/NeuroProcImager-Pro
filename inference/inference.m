function inference(data_file)
    ntasks = feature('numcores'); % number of cpus
    
    %%
    data_path = './data';
    
    output_dir = './output';
    
    fs = 400; % sampling rate data is sent to inference
    save_fs = 150; % sampling rate estimates are saved
    
    %% local inference
    [xi, s_y, A, B, C, Q, R, H, varsigma, v0, N_states, N_syn, N_inputs, N_samples] = forward_simulate_nmm([]);
    Q(10:13, 10:13) = 1.5625*eye(4, 4); 
    
    xi0_empirical = mean(xi(:,N_samples/2:end),2);
    P0_empirical = 10*cov(xi(:,N_samples/2:end)');
    P0_empirical(2*N_syn+1:end,2*N_syn+1:end) = eye(N_syn+N_inputs)*10e-2;
    
    i_file = data_file;
    moduler(ntasks, i_file, data_path, output_dir, s_y, xi0_empirical, P0_empirical, A, B, C, Q, R, H, varsigma, v0, N_syn, N_states, N_inputs, save_fs, fs);

    %% inter-regional connectivity inference
    load([output_dir '/Seizure_' num2str(i_file) '_estimates.mat'], 'xi_hat_list') % load data
    xi_hat_list(isnan(xi_hat_list))=0; % replace nan with 0
    
    n_channels = size(xi_hat_list,3);
    n_timesteps = size(xi_hat_list,2);

    scale = 50;
    window_size_inSeconds = 0.1;
    window_length = window_size_inSeconds*save_fs;

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

    %% Calculate connection matrix A (mu = A*g(vp))
    windows_num = floor(n_timesteps/window_length)-1;
    
    kernel = arrayfun(@h,0:1/save_fs:window_size_inSeconds); % discretized kernel

    A_sequence = zeros(n_channels, n_channels, windows_num);

    for i = 1:windows_num
        mu = input(:,(i-1)*window_length+1:i*window_length)*scale;
        vp = v_pyr(:,(i-1)*window_length+1:i*window_length);
        gvp = g_function(vp,v0,varsigma);
        
        for j = 1:size(gvp,1) % convolution of the gvp and the kernel
            convolved = conv(gvp(j,:),kernel);
            gvp(j,:) = convolved(1:size(gvp,2));
        end
        
        A_sequence(:,:,i) = analytic_multiregress(gvp', mu');
    end

    save([output_dir '/Seizure_' num2str(i_file) '_connectivity_matrix.mat'], 'A_sequence', '-v7.3')
end

%% Subfunction
function moduler(ntasks, i_file, data_path, output_dir, s_y, xi0_empirical, P0_empirical, A, B, C, Q, R, H, varsigma, v0, N_syn, N_states, N_inputs, save_fs, fs)
    [file_number] = file_i2file_number(i_file);
    file_path = [data_path '/Seizure_' file_number '.mat'];
    load(file_path, 'data')

    y = (std(s_y)/max(std(data(1:12000,:))))*data; % scale signal to the same magnitude of the simulation

    prediction_error = nan(size(y,1), 16);

    %% inference starts here
    poolobj = parpool(ntasks);
    parfor (iCh = 1:size(y,2), ntasks)
        measurement = y(:,iCh);      
        [xi_hat, ~] = AKF_quick(measurement, xi0_empirical, P0_empirical, A, B, C, Q, R, H, varsigma, v0, N_syn, N_states, N_inputs);

        prediction_error(:,iCh) = (H*xi_hat(:,2:end) - measurement');
        xi_hat_list(:,:,iCh) = single(resample(xi_hat', save_fs, fs)');
    end
    delete(poolobj)

    %% save
    save([output_dir '/Seizure_' num2str(i_file) '_estimates.mat'], 'xi_hat_list', '-v7.3');
end