%% Demo: correct synthesized dual-comb data with the Augmented Kalman filter
% 1. Generates artificial data similar to Fig. 1 of Burghoff et al. (Generate_Data.m)
% 2. Corrects it by running the main filter (Augmented_Kalman.m)
% 3. Generates plots assessing the efficacy of the correction (Efficacy_Plots.m)
thisdir = fileparts(mfilename('fullpath')); if isempty(thisdir), thisdir = pwd; end
root = fileparts(thisdir);
addpath(thisdir, root, fullfile(root,'dependencies'));

Generate_Data;

%% Setup for solver
data=struct();
data.sn = Vn;                            % input signal
data.dt = dt;                            % input timestep

% Parameters (defaults come from AK_Params; we override the few that matter)
param = AK_Params();
param.fD_range = fDm*[0.75,1.25];        % Frequency range the rep rate is known to be in.
                                         % If you don't know, get it from abs(fft(Vn)).^2.
param.Q = [0.0419 0.000782 -0.00315 0.000438;0.000782 0.000146 -0.00199 -1.1e-05;-0.00315 -0.00199 0.328 -0.000572;0.000438 -1.1e-05 -0.000572 0.00056];
                                         % Process noise matrix (4x4).
                                         % If you don't know, set to diag([1,.01,1,.01].^2) and run EM.
param.excess_noise = 2.76;               % Excess measurement noise factor.
                                         % If you don't know, set to 1 and run EM.

% Post-regularization: 'lpf' is the rubber-band band-limited filter (Xiao &
% Burghoff 2025), low-pass filtering the phase without edge artifacts. With
% uniform valid-region weighting it puts the residuals on the coherent limit.
% ('ekf' uses the one-step EKF, which is preferable for real data -- see
% Demo_MIR_QCL_DCS / Demo_THz_QCL_DCS.) 'lpf' is the AK_Params default, so this
% line is optional.
param.post_regularization = 'lpf';

% Optional: supply the true rep rate / offset so they are plotted alongside
% the estimate (NOT used by the correction itself).
param.knownfD = fD;
param.knownf0 = f0;

oK = Augmented_Kalman(data,param);       % Run the filter!
Efficacy_Plots
