%% Demo: correct real THz QCL dual-comb data
% Loads a window of a terahertz quantum-cascade-laser multiheterodyne
% dual-comb dataset, forms the complex signal from its I/Q channels, and
% corrects it with the Augmented Kalman filter, learning the noise
% parameters by EM.
%
% As with Demo_MIR_QCL_DCS, real data has no ground truth; efficacy is judged
% by the comb lines emerging from the noise floor in the corrected PSD.
thisdir = fileparts(mfilename('fullpath')); if isempty(thisdir), thisdir = pwd; end
root = fileparts(thisdir);
addpath(thisdir, root, fullfile(root,'dependencies'), fullfile(root,'data'));

%% Load a window of raw I/Q data
d = load(fullfile(root,'data','THz_QCL_DCS.mat'));
dt = d.Tinterval;
chanI = 'C'; chanQ = 'D';

total_time = 100e-6;             % amount of time to process (s)
start_time = 0e-6;               % where to start processing (s)

N  = round(total_time/dt);
i0 = floor(start_time/dt)+1;
sI = d.(chanI); sQ = d.(chanQ);
sn = sI(i0:i0+N-1,1) + 1i*sQ(i0:i0+N-1,1);
N  = length(sn);

%% Setup for solver
data = struct();
data.sn = sn;
data.dt = dt;

param = AK_Params();
param.fD_range = 35e6*[0.86,1.14];

% Tuned process / measurement noise for this dataset (from a prior EM run).
param.Q = [0.222 0.00148 0.0299 0.00014;0.00148 7.9e-05 0.00104 2.59e-06;0.0299 0.00104 0.0622 0.000293;0.00014 2.59e-06 0.000293 4.91e-06];
param.excess_noise = 3.18;

% One-step EKF post-regularization + EKF-based EM (best for real data).
param.post_regularization = 'ekf';
param.EM = 'ekf';

%% Run the filter with EM
if ~exist('nEM','var'), nEM = 10; end
[Q, R, oK, em] = Augmented_Kalman_EM(data, param, nEM);

fprintf('THz_QCL_DCS: learned excess_noise = %.4g\n', R);
