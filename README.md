# NeuroProcImager-Pro
#### Yun Zhao, Levin Kuhlmann (Monash University, Australia), Email: yun.zhao@monash.edu, levin.kuhlmann@monash.edu

**NeuroProcImager-Pro** was used in the below manuscripts 

1. *Inference-based time-resolved cortical stability and chaos analysis for focal epileptic seizures*,

2. *Cortical local dynamics, connectivity and stability correlates of global conscious states*.

**NeuroProcImager-Pro**, an extension of **NeuroProcImager** ([Github link](https://github.com/yundumbledore/NeuroProcImager/tree/main), [Neuroimage Paper link](https://www.sciencedirect.com/science/article/pii/S1053811922007078)), explores neurophysiological underpinnings of brain functions via analyzing
1. *dynamic cortical stability*
2. *dynamic cortical chaos*

The below demonstrates the methods of **NeuroProcImager-Pro** and its application to analyze cerebral cortical dynamics for focal epileptic seizures. 


## Methods

### Modelling the brain
![](assets/Framework.png)
The schematic of the multiple regions model fitted to 16-channel iEEG data is shown. The left side shows the device used to collect iEEG recording and an example of iEEG recording. The middle shows a neural mass model (NMM) and the right shows the multiple regions model. Each node in the multiple regions model is a NMM. Purple lines represent connections between NMMs. Each iEEG time series is fitted with a NMM. In this study, iEEG time series are used whereas our framework also applies to EEG and MEG data.

The NMM comprises three neural populations, namely excitatory (e), inhibitory (i), and pyramidal (p). The pyramidal population (in infragranular layers) driven by the external input $\mu$, excites the spiny stellate excitatory population (in granular layer IV) and inhibitory interneurons (in supragranular layers), and is excited by the spiny stellate excitatory population and inhibited by the inhibitory interneurons. Neural populations are characterized by their time varying mean (spatial, not time averaged) membrane potential, $v_n$ , which is the sum of contributing population mean post-synaptic potentials, $v_{mn}$ (pre-synaptic and post-synaptic neural populations are indexed by $m$ and $n$) and connected via synapses in which the parameter, $\alpha_{mn}$ quantifies the population averaged connection strength. $\alpha_{mn}$ is referred as the regional neurophysiological variables. In the figure above, there are four variables in a cortical region.

Coupling of two cortical regions is achieved by connecting the output of the pyramidal population in one region to the input of the pyramidal population in another region via a synapse. The synaptic connection strength from region $a$ to $b$ is referred as the inter-regional connectivity parameter $w_{ab}$. For the multiple regions model, the input to the pyramidal population at one cortical region is formed by the combination of post-synaptic membrane potentials induced by the output of the other regions.

### Parameter estimation
To estimate parameters of the whole-cortex model from data, we first treat each NMM in the multiple regions model independent and apply the semi-analytical Kalman filter (AKF) that we developed in [Neuroimage Paper link](https://www.sciencedirect.com/science/article/pii/S1053811922007078) to estimate parameters of each NMM. The AKF is an unbiased estimator, providing the minimum mean square error estimates for model parameters, under the assumption that the underlying probability distribution of the model state is Gaussian. Briefly, the aim of the estimation is to calculate the posterior distribution of model parameters at time point $t$ given measurements up to $t$. This gives time-varying parameter estimates.

To estimate inter-regional connectition strength $w$, we use the multivariate regression model to relate the input of the pyramidal population $\mu$ in one region to the output of the pyramidal population in other regions. The multivariate regression gives an estimate for $w$ and a bias term $e$ representing the contribution from thalamus input. The estimation is done in each time window, which result in time-varying $w$ and $e$ estimates.

### Dynamic cortical stability analysis
Linear stability analysis is performed in each time window where the averaged parameter estimates are used to define the system. Linear stability analysis is a mathematical technique used to determine the stability of equilibrium points in a dynamical system. The process involves linearizing the system of differential equations around the equilibrium point by computing the Jacobian matrix of the system's partial derivatives. This matrix represents the linear approximation of the system near the equilibrium. By analyzing the eigenvalues of the Jacobian, one can infer the stability of the equilibrium: if all eigenvalues have negative real parts, the equilibrium is stable; if any eigenvalue has a positive real part, the equilibrium is unstable. This method provides insights into the local behavior of the system and helps in understanding how small perturbations evolve over time.

![](figures/Seizure_3_dynamic_stability.png)
The figure above shows a dynamic stability analysis for an example seizure recording. The top panel shows the time-varying Jocobi's eigenvalue spectrum. The color indicates the log scaled number of eigenvalue. The bottom panel shows the number of unstable eigenmode (with positive eigenvalues) as a function of time.

### Dynamic cortical chaos analysis
Chaos analysis is performed in each time window where the averaged parameter estimates are used to define the system. Chaos analysis in a dynamical system involves examining the system's behavior to identify the presence of chaotic dynamics, which are characterized by sensitivity to initial conditions, aperiodicity, and long-term unpredictability. To perform chaos analysis, one typically begins by computing Lyapunov exponents, which measure the rate of separation of infinitesimally close trajectories; a positive Lyapunov exponent indicates chaos. Additionally, tools such as phase space reconstruction and Poincaré sections can be used to visualize the system's trajectory and detect chaotic attractors. The bifurcation diagram is another useful tool, showing how the system's behavior changes as a parameter is varied, often revealing transitions from periodic to chaotic states. Lastly, the correlation dimension and fractal dimension can quantify the complexity of the system's attractors, further confirming chaotic behavior. These methods together provide a comprehensive understanding of the system's dynamics and the presence of chaos.

![](figures/Seizure_3_dynamic_chaos.png)
The figure above shows a dynamic chaos analysis for an example seizure recording. The top panel shows the time-varying Lyapunov spectrum. The color indicates the magnitude of exponent. The bottom panel shows the maximal Lyapunov exponent as a function of time.

## Adaptation to your data
The package can be applied to your dataset. Steps are as follows:
1. Save your multi-channel iEEG/EEG/MEG data into a matlab data file with the format t x N where t is the number of time points and N is the number of channels,
2. Name it for example Seizure_009.mat and put it in the data folder,
3. Open run.m and change the value of 'data_file' to be the number seizure number for example 9,
4. Firstly estimate multi-region model parameters using the keyword 'parameter estimation',
5. Then, choose an analysis from 'stability analysis' and 'chaos analysis'.
