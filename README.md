# Turbo Comb
A toolbox for correcting dual comb spectroscopy data

This code is intended to demonstrate the algorithm described in
"Generalized method for the computational phase correction of arbitrary dual comb signals" by Burghoff et al.

To run a basic demonstration of the technique, run Demo_Augmented_Kalman.m.
This does the following:

1. Generates artificial data similar to Fig. 1. (Generate_Data.m)
2. Corrects it by running the main script. (Augmented_Kalman.m)
3. Generates plots assessing the efficacy of the correction (Efficacy_Plots.m)

Demo_Augmented_Kalman_EM is similar, but demonstrates the expectation maximization process:

1. Generates artificial data similar to Fig. 1. (Generate_Data.m)
2. Corrects it by running the main script. (Augmented_Kalman.m), using an arbitrary choice of
	param.Q = diag([1,.01,1,.01].^2);
	param.excess_noise = 1;
   Note that param.Q is expressed in the order of [\phi_0,\phi_r,f_0,f_r].
   	The two phases have units of radians, while the two frequencies have units of (Rep rate)^2/(Rep period).
   	In other words, Q_{f0,f0}=1 means that the offset fluctuates by one f_r per T_r.
   param.excess_noise is a multiplier for the measurement noise \sigma.
	The code takes the top 10% of the signal's frequency range, uses the PSD there as a baseline white
  noise level, and then multiplies that by param.excess_noise. The optimal value is usually larger than 1.
3. Generates plots assessing the efficacy of the initial correction (Efficacy_Plots.m).
4. Runs EM, plotting the convergence of the process.
5. Generates plots assessing the efficacy of the final correction (Efficacy_Plots.m).
