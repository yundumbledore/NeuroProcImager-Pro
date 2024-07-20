function [A, B, C, D, E, v0, varsigma] = set_network_params(model_parameters)
%% define model
input = 300;
input_offset = [];
fs = 400;
TimeOfSim = 10;

scale = 50;                         % this is to get states and derivatives on the same order of magnitude

% set units of membrane potentials (1 for mV, 0 for V)
%
mV = 1;
V2mVfactor = 1e3;

%
dt = 1/fs;                     
N_samples = round(TimeOfSim/dt);


if ~isempty(input_offset)
    N_inputs = 2;
else
    N_inputs = 1; 
end

N_model = 16;
N_syn = 5;                         
N_states = N_model*(N_syn*2) + N_model*(N_syn-1) + N_model*N_model;

% %% define the disturbance covariance matrix
% %
% sigma_all = 5e-8;                               % something small for all states
% sigma_input = 5e-4;                             % for input
% sigma_params = 5e-5;%sigma_all;
% sigma_offset = 5e-6;
sigma_R = 1e-3;
% 
% Q = eye(N_states)*(scale*sqrt(dt)*sigma_all)^2;                             % add a tiny bit of noise to all states (added for numerical stability)
% Q(2*N_syn+1:end,2*N_syn+1:end) = eye(N_syn+N_inputs)*(scale*sqrt(dt)*sigma_params)^2;
% % **** HARDCODED NUMBERS HERE
% Q(2*N_syn+1,2*N_syn+1) = (scale*sqrt(dt)*sigma_input)^2;
% if N_inputs > 1
%     Q(2*N_syn+2,2*N_syn+2) = (scale*sqrt(dt)*sigma_offset)^2;
% end
% 
% % measurement disturbance covariance
% %
R = sigma_R^2;

%% General parameters from J&R
%
% sigmoid bits
%
f_max = 2.5;  % maximum firing rate (spikes/s)
r = 560; 
varsigma = 1.699/r;     % slope of the sigmoid                                                     % (spikes/(Vs))
varsigma_sq = varsigma^2;
v0 = 0.006;    % mean firing threshold                                                             % (V)

% synaptic gains
%
alpha_e = 3.25e-3;                                                          % gain of excitatory synapses (V)
alpha_i = -22e-3;                                                           % gain of inhibitory synapses (V)


% synaptic kernel time constants
%
ex_tau = 0.010;                     % excitatory synaptic time constant (s)
in_tau = 0.020;                     % inhibitory synaptic time constant (s)
d_tau = 1/33;

% input to py population
%
input = input*scale;                   % ! scale here ! the input is a membrane potential input the pyramidal population this is similar to setting it at 270
input = input * alpha_e/ex_tau * ex_tau^2; % transformed input
%       ~~~~~   ~~~~~~~~~~~~~~   ~~~~~~~~~
%       input   synaptic gain    integral of kernel

% measurement DC offset
input_offset = input_offset * scale;        
input_offset = input_offset * alpha_e/ex_tau * ex_tau^2;


if mV == 1   % adjust units
%     Q = V2mVfactor^2 * Q;
    R = V2mVfactor^2 * R;
    
    r = r/V2mVfactor;
    varsigma = 1.699/r;                                                     % (spikes/(Vs))
    varsigma_sq = varsigma^2;
    v0 = v0*V2mVfactor;
    alpha_e = alpha_e*V2mVfactor;                                           % gain of excitatory synapses (V)
    alpha_i = alpha_i*V2mVfactor;                                           % gain of inhibitory synapses (V)
    
    input= input*V2mVfactor;
    input_offset = input_offset*V2mVfactor;
end

% conectivity constants to relate to Jansen and Rit 1995 model
%
% ConnectivityConst = 270;                            % Jansen and Rit connectivity parameters. Either 135, 270 or 675
% C1 = ConnectivityConst;
% C2 = 0.8*ConnectivityConst;
% C3 = 0.25*ConnectivityConst;
% C4 = 0.25*ConnectivityConst;

%% this is the observation function.
%
H = zeros(N_model,N_states);        %Initialize to zeros and later add 1s to states that contribute to EEG

% initialize adjancy matrix
%
Gamma = zeros(N_states, N_states);   %  - plus 1 for input

% specify synapses
%
syn_index = 0;

% syn1, connection from I to P
%
syn_index = syn_index + 1;
tau(syn_index) = in_tau;
% alpha(syn_index) = alpha_i*2*f_max*C4*dt / tau(syn_index);          % note the time constant and time step are in the gains
presyn_inputs = 2;                                                  % the presynaptic population is getting inputs from synapses 2
if ~isempty(presyn_inputs)
    Gamma(2*syn_index,2*presyn_inputs-1) = 1;                       % set the entries of Gamma corresponding to indices of presynaptic inputs to 1
end
H(2*syn_index-1) = 1;

% syn2, connection from P to I
%
syn_index = syn_index + 1;
tau(syn_index) = ex_tau;
% alpha(syn_index) = alpha_e*2*f_max*C3*dt / tau(syn_index);          % note the time constsnt and time step are in the gains
presyn_inputs = [1 4 5];                                            % the presynaptic population is getting inputs from synapses 1, 4, 5
if ~isempty(presyn_inputs)
    Gamma(2*syn_index,2*presyn_inputs-1) = 1;
end
H(2*syn_index-1) = 0;                                               % set to one if it contributes to the EEG (i.e. if the synapse is to Py cells)

% syn3, connection from P to E
%
syn_index = syn_index + 1;
tau(syn_index) = ex_tau;
% alpha(syn_index) = alpha_e*2*f_max*C1*dt / tau(syn_index);
presyn_inputs = [1 4 5];                                         	% the presynaptic population is getting inputs from no other synapses (in the model)
if ~isempty(presyn_inputs)
    Gamma(2*syn_index,2*presyn_inputs-1) = 1;
end
H(2*syn_index-1) = 0;

% syn4, connection from E to P
%
syn_index = syn_index + 1;
tau(syn_index) = ex_tau;
% alpha(syn_index) = alpha_e*2*f_max*C2*dt / tau(syn_index);          % note the time constsnt and time step are in the gains
presyn_inputs = 3;                                                  % the presynaptic population is getting inputs from synapse 3
if ~isempty(presyn_inputs)
    Gamma(2*syn_index,2*presyn_inputs-1) = 1;
end
H(2*syn_index-1) = 1;

% syn5
%
syn_index = syn_index + 1;
tau(syn_index) = d_tau;
% alpha(syn_index) = alpha_e*2*f_max*C2*dt / tau(syn_index);          % note the time constsnt and time step are in the gains
presyn_inputs = 3;                                                  % the presynaptic population is getting inputs from synapse 3
if ~isempty(presyn_inputs)
    Gamma(2*syn_index,2*presyn_inputs-1) = 1;
end
H(2*syn_index-1) = 1;

% for input
%
syn_index = syn_index + 1;
H(2*syn_index-1) = 1;           % the input contributes to the observation function

if N_inputs > 1
    % offset term
    H(2*syn_index) = 1;            % offset contributes to the observation function
end

% rescale H
%
H = H/scale;                                                        % !scale! this is help deal with our numerical issues.

% set_ABC.m

% dean freestone

% this script take the parameters from set_params and creates the system
% matrices for the neural mass model (recurrent NN)

%% define A
% A encodes the dynamics induced by the membrane time constants.
% A is made up of the submatrices Psi in a block diagonal structure.
% There is a Psi for each connection in the model. This is where all the
% synaptic time constants enter the system. Further, the scale paramter
% enters here (and with C (multiplicative factor) and with the H (divisor).

Psi = zeros(2*N_syn,2*N_syn);               % initialise Psi, the component of A for fast states
for n=1:N_syn                               % build block diagonal structure
    index = 2*(n-1)+1;
    Psi(index:index+1,index:index+1) = [0 scale ; -1/(scale*tau(n)^2) -2/(tau(n))];
end

% Psi = [Psi zeros(2*N_syn,1); zeros(1,2*N_syn) -(1/tau_s)];

A = zeros(N_states, N_states);
for i = 1:N_model
    A(1+(i-1)*(2*N_syn):2*N_syn+(i-1)*(2*N_syn), 1+(i-1)*(2*N_syn):2*N_syn+(i-1)*(2*N_syn)) = eye(2*N_syn) + dt*Psi;
end
A(N_model*(2*N_syn)+1:end, N_model*(2*N_syn)+1:end) = eye(N_states-N_model*(2*N_syn));


B = zeros(N_states, N_states);
for i = 1:N_model
    B(2+(2*N_syn)*(i-1), N_model*(2*N_syn)+1+(i-1)*(N_syn-1)) = 1;
    B(4+(2*N_syn)*(i-1), N_model*(2*N_syn)+2+(i-1)*(N_syn-1)) = 1;
    B(6+(2*N_syn)*(i-1), N_model*(2*N_syn)+3+(i-1)*(N_syn-1)) = 1;
    B(8+(2*N_syn)*(i-1), N_model*(2*N_syn)+4+(i-1)*(N_syn-1)) = 1;
end

C = zeros(N_states, N_states);
c = [0	0	0	0	0	0	0	0	0	0;
0	0	1	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0;
1	0	0	0	0	0	1	0	1	0;
0	0	0	0	0	0	0	0	0	0;
1	0	0	0	0	0	1	0	1	0;
0	0	0	0	0	0	0	0	0	0;
0	0	0	0	1	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0];
for i = 1:N_model
    C(1+(i-1)*(2*N_syn):(2*N_syn)+(i-1)*(2*N_syn), 1+(i-1)*(2*N_syn):(2*N_syn)+(i-1)*(2*N_syn)) = c;
end
C = C/scale;

D = zeros(N_states, N_states);
ws = reshape(model_parameters(N_model*4+1:end), [16,16]);
for i = 1:N_model
    for j = 1:N_model
        D((i-1)*(2*N_syn)+10,(j-1)*(2*N_syn)+10) = ws(j,i);
    end
end
D = D/scale;

E = zeros(N_states, N_states);
for i = 1:N_model
    E((i-1)*(2*N_syn)+10,(i-1)*(2*N_syn)+1:i*(2*N_syn)) = [1	0	0	0	0	0	1	0	1	0];
end
E = E/scale;

% H = zeros(N_model, N_states);
% for i = 1:N_model
%     H(i,(i-1)*(2*N_syn+1)+1:i*(2*N_syn+1)) = [1	 0	0	0	0	0	1	0	1	0	1];
% end
% H = H/scale;

%%
% xi = zeros(N_states,N_samples);
% 
% xi(:,1) = [zeros(N_model*(2*N_syn+1),1); zeros(N_model,1); model_parameters];

% if insert_where == 0
%     xi(N_model*(2*N_syn+1)+1:N_model*(2*N_syn+1)+N_model,:) = [zeros(N_model,fs*TimeOfSim/2) intensity*ones(N_model,npoints) zeros(N_model, fs*TimeOfSim/2-npoints)];
% else
%     xi(N_model*(2*N_syn+1)+insert_where,:) = [zeros(1,fs*TimeOfSim/2) intensity*ones(1,npoints) zeros(1, fs*TimeOfSim/2-npoints)];
% end

% for n = 1:N_samples-1
%     in_phi = g(C*xi(:,n), v0, varsigma);
%     ex_phi = g(E*xi(:,n), v0, varsigma);
%     tmp = A*xi(:,n) + B*xi(:,n).*in_phi + D*ex_phi;
%     xi(1:N_model*(2*N_syn+1),n+1) = tmp(1:N_model*(2*N_syn+1));
% end
% 
% v = sqrt(R)*randn(N_model,N_samples); % variance
% y = H*xi;% + v; 
end