function chaos_analysis(data_file)
    %% Settings
    window_length = 0.1*150;
    w = 5e2; %*******************
    scale = 50;
    fs = 400;
    
    % Integration parameters
    T = 1;                 % Total integration time
    dt = 1/400;              % Time step
    num_samples = 1;        %*******************
    initial_perturbation = 1e-5;
    
    % Start parallel
    ntasks = feature('numcores'); % num of CPUs
    poolobj = parpool(ntasks);
    
    i_file = data_file;
    
    connectivity_file_path = ['./output/Seizure_' num2str(i_file) '_connectivity_matrix.mat'];
    
    %% Prepare model
    % System dimensionality (number of variables, parameters)
    n = 160;
    n_theta = 320;
    
    % State vector for the Jansen-Rit model
    load(['./output/Seizure_' num2str(i_file) '_estimates.mat'], 'xi_hat_list') % load data
    load(connectivity_file_path, 'A_sequence')
    n_timesteps = size(xi_hat_list,2);
    
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
    
    windows_num = floor(n_timesteps/window_length)-1;
    time_varying_lyapunov_exponents = zeros(n, windows_num);
    clear xi_hat_list
    
    
    % parfor_progress(windows_num);
    parfor (i = 1:windows_num, ntasks)
    % for i = 1:windows_num
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
    
    x = [v_ip_avg z_ip_avg v_pi_avg z_pi_avg v_pe_avg z_pe_avg v_ep_avg z_ep_avg mu_avg zeros(16,1)];
    x = reshape(x',[160,1]);
    alphas = [alpha_ip_avg alpha_pi_avg alpha_pe_avg alpha_ep_avg];
    alphas = reshape(alphas',[64,1]);
    ws = A_sequence(:,:,i);
    ws = reshape(ws, [256,1]);
    model_parameters = [alphas; ws];
    
    % Parameters for the Jansen-Rit network
    [A, B, C, D, E, v0, varsigma] = set_network_params(model_parameters);
    
    %% Computation
    % Initialize variables to store Lyapunov exponents
    lyapunov_exponents = zeros(n+n_theta, num_samples);
    initial_conditions = zeros(n+n_theta, num_samples);
    
    % Loop over different initial conditions
    
    for sample = 1:num_samples
        % Set random initial conditions
        x_original = [x; model_parameters];
        variations = [-w + 2*w*rand(n,1); zeros(n_theta,1)];
        
        if sample == 1
            initial_conditions(:, sample) = x_original;
        else
            initial_conditions(:, sample) = x_original + variations;
            x_original = x_original + variations;
        end
    
        % Perturb the initial condition in different directions
        perturbations = [eye(n) zeros(n, n_theta); zeros(n_theta, n+n_theta)] * initial_perturbation;
        
        % Initialize the perturbed trajectories
        x_perturbed = repmat(x_original, 1, n+n_theta) + perturbations;
    
        % Numerical integration loop
        for t = 0:dt:T
            % Numerical integration step for both trajectories
            in_phi = g(C*x_original, v0, varsigma);
            ex_phi = g(E*x_original, v0, varsigma);
            x_original = A*x_original + B*x_original.*in_phi + D*ex_phi;
    
            in_phi = g(C*x_perturbed, v0, varsigma);
            ex_phi = g(E*x_perturbed, v0, varsigma);
            x_perturbed = A*x_perturbed + B*x_perturbed.*in_phi + D*ex_phi;
    
            % Calculate the separation in Euclidean space for each perturbation
            separation = sqrt(sum((x_original - x_perturbed).^2, 2)) + 1e-36;
    
            % Calculate the Lyapunov exponents for each direction
            lambda = log(separation / initial_perturbation);
            lambda(n+1:end) = 0;
            lyapunov_exponents(:, sample) = lyapunov_exponents(:,sample) + lambda;
    
            % Renormalization step (optional but recommended)
            norm_perturbation = sqrt(sum((x_perturbed - x_original).^2, 2)) + 1e-36;
            x_perturbed(1:n,1:n+n_theta) = x_original(1:n) + (x_perturbed(1:n,1:n+n_theta) - x_original(1:n)) ./ norm_perturbation(1:n);
        end
    end
    
    % Compute the average Lyapunov exponents over all samples
    average_lyapunov_exponents = nanmean(lyapunov_exponents, 2);
    average_lyapunov_exponents = average_lyapunov_exponents(1:n);
    time_varying_lyapunov_exponents(:,i) = average_lyapunov_exponents;
    end
    
    %% save Lya file
    save(['./output/Seizure_' num2str(i_file) '_Lyapunov_exponents.mat'],'time_varying_lyapunov_exponents','-v7.3')
    delete(poolobj) 
end
