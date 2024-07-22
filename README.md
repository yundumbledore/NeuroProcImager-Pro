# NeuroProcImager-Pro
#### Yun Zhao, Levin Kuhlmann (Monash University, Australia), Email: yun.zhao@monash.edu, levin.kuhlmann@monash.edu

**NeuroProcImager-Pro** is associated with the manuscript *Inference-based time-resolved cortical stability and chaos analysis for focal epileptic seizures*.

**NeuroProcImager-Pro**, an extension of **NeuroProcImager** ([Github link](https://github.com/yundumbledore/NeuroProcImager/tree/main), [Neuroimage Paper link](https://www.sciencedirect.com/science/article/pii/S1053811922007078)), explores neurophysiological underpinnings of brain functions via analyzing
1. *dynamic cortical stability*
2. *dynamic cortical chaos*


## Methods

### Modelling the brain
![](assets/Framework.png)
The schematic of the multiple regions model fitted to 16-channel iEEG data is shown. The left side shows the device used to collect iEEG recording and an example of iEEG recording. The middle shows a neural mass model (NMM) and the right shows the multiple regions model. Each node in the multiple regions model is a NMM. Purple lines represent connections between NMMs. Each iEEG time series is fitted with a NMM. In this study, iEEG time series are used whereas our framework also applies to EEG and MEG data.

Briefly, the NMM comprises three neural populations, namely excitatory (e), inhibitory (i), and pyramidal (p). The pyramidal population (in infragranular layers) driven by the external input $\mu$, excites the spiny stellate excitatory population (in granular layer IV) and inhibitory interneurons (in supragranular layers), and is excited by the spiny stellate excitatory population and inhibited by the inhibitory interneurons. Neural populations are characterized by their time varying mean (spatial, not time averaged) membrane potential, $v_n$ , which is the sum of contributing population mean post-synaptic potentials, $v_{mn}$ (pre-synaptic and post-synaptic neural populations are indexed by $m$ and $n$) and connected via synapses in which the parameter, $\alpha_{mn}$ quantifies the population averaged connection strength. $\alpha_{mn}$ is referred as the regional neurophysiological variables. In the figure above, there are four variables in a cortical region.

Coupling of two cortical regions is achieved by connecting the output of the pyramidal population in one region to the input of the pyramidal population in another region via a synapse. The synaptic connection strength from region $a$ to $b$ is referred as the inter-regional connectivity parameter $w_{ab}$. For the multiple regions model, the input to the pyramidal population at one cortical region is formed by the combination of post-synaptic membrane potentials induced by the output of the other regions.

### Parameter estimation
To estimate parameters of the whole-cortex model from data, we first treat each NMM in the multiple regions model independent and apply the semi-analytical Kalman filter (AKF) that we developed in [Neuroimage Paper link](https://www.sciencedirect.com/science/article/pii/S1053811922007078) to estimate parameters of each NMM. The AKF is an unbiased estimator, providing the minimum mean square error estimates for model parameters, under the assumption that the underlying probability distribution of the model state is Gaussian. Briefly, the aim of the estimation is to calculate the posterior distribution of model parameters at time point $t$ given measurements up to $t$. This gives time-varying parameter estimates.

To estimate inter-regional connectition strength $w$, we use the multivariate regression model to relate the input of the pyramidal population $\mu$ in one region to the output of the pyramidal population in other regions. The multivariate regression gives an estimate for $w$ and a bias term $e$ representing the contribution from thalamus input. The estimation is done in each time window, which result in time-varying $w$ and $e$ estimates.

### Dynamic cortical stability analysis
This analysis is performed in each time window where the averaged parameter estimates are used to define the system. 

Equilibrium point of the system is found in each time window. To study the behaviour of the model around the equilibrium point, we linearize the model and calculate its Jacobian matrix $\mathbf{J}$ at the equilibrium point. We calculate the eigenvalues $\lambda_k$ of $\mathbf{J}$ to see if the system is stable around the equilibrium point. Precisely, we perform the eigenvalue decomposition $\mathbf{J} = U \Lambda U^\top$, where $U$ is orthonormal, $\Lambda$ is the diagonal eigenvalue matrix, and $^\top$ means the transpose conjugate. We call the real part of $\lambda_k$ the criticality index. If a criticality index is less than 0, the corresponding mode is stable. A small perturbation along the eigenvector will decay and the system can return to the equilibrium point. Conversely, if a criticality index is greater than 0, even a trivial perturbation along the eigenvector will diverge the system from the equilibrium point. When the criticality index equals 0, the system is at the neutral point along the corresponding eigenvector. In short, the system is unstable if there is at least one unstable eigenmode, and is stable when all eigenmodes are stable.

To visualize the result, we show the time course of the eigenvalue spectrum and indicate the number of unstable eigenvalues at each time point. Each unstable eigenvalue represents a distinct direction in which the system's state can diverge from the equilibrium point. The more unstable eigenvalues there are, the more directions the system has in which it can become unstable. The magnitude (specifically, the real part) of an unstable eigenvalue indicates how rapidly perturbations in that direction will grow. Systems with larger unstable eigenvalues may be more sensitive to initial conditions and external disturbances.

## Demonstration


## Adaptation to your data


## Common issues and solutions
Issue: *Errors with parallel pool*.

Solution: This type of errors occur when the machine runs out of the memory or loses connections to some cores. As the parallel computation heavily relies on the hardward of computers, it is suggested to reduce the number of parallel threads by defining the number of threads you want to run in parallel in line #16 in /Sources/calculate_dynamic_stability.m. By default, all cores are used in computation.
