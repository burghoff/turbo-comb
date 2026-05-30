%% Demo: Expectation Maximization of the process/measurement noise
% Same synthesized data as Demo_Augmented_Kalman, but here we start from an
% arbitrary guess of the noise parameters and let EM learn them:
%   param.Q = diag([1,.01,1,.01].^2);   param.excess_noise = 1;
% Augmented_Kalman_EM runs the filter, re-estimates Q and excess_noise
% each iteration, and plots the convergence.
%
% param.Q is expressed in the order [phi_0, phi_r, f_0, f_r]. The two phases
% have units of radians; the two frequencies have units of (rep rate)^2/(rep
% period). Q(f0,f0)=1 means the offset fluctuates by one f_r per T_r.
% param.excess_noise multiplies the white-noise floor (top 10% of the
% frequency range) to set the measurement noise sigma.
thisdir = fileparts(mfilename('fullpath')); if isempty(thisdir), thisdir = pwd; end
root = fileparts(thisdir);
addpath(thisdir, root, fullfile(root,'dependencies'));

Generate_Data;

%% Setup for solver
data=struct();
data.sn = Vn;                            % input signal
data.dt = dt;                            % input timestep

param = AK_Params();
param.fD_range = fDm*[0.75,1.25];        % Frequency range the rep rate is known to be in.
param.Q = diag([1,.01,1,.01].^2);        % Arbitrary initial process-noise guess.
param.excess_noise = 1;                  % Arbitrary initial measurement-noise factor.
param.post_regularization = 'lpf';       % Rubber-band smoother (see Demo_Augmented_Kalman).
param.knownfD = fD;                      % Plotted alongside the estimate (NOT used by the correction).
param.knownf0 = f0;

%% Run EM (10 iterations). Returns the learned Q and noise factor.
if ~exist('nIters','var'), nIters = 10; end
[Q, R, oK, em] = Augmented_Kalman_EM(data, param, nIters);

fprintf('Learned excess_noise = %.4g\n', R);
disp('Learned Q ='); disp(Q);

%% Assess the final correction
Efficacy_Plots
