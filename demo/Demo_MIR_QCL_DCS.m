%% Demo: correct real mid-IR QCL dual-comb spectroscopy data
% Loads a single Picoscope block of mid-infrared quantum-cascade-laser
% dual-comb data (df ~ 2.5 kHz), down-selects the comb band, and corrects it
% with the Augmented Kalman filter, learning the noise parameters by EM.
%
% Real data has no ground truth, so efficacy is judged by how far the comb
% lines rise out of the noise in the corrected PSD (see the figures produced
% by Augmented_Kalman_EM: uncorrected vs corrected).
thisdir = fileparts(mfilename('fullpath')); if isempty(thisdir), thisdir = pwd; end
root = fileparts(thisdir);
addpath(thisdir, root, fullfile(root,'dependencies'), fullfile(root,'data'));

%% Load one block of raw data
d = Picoscope_6404D(fullfile(root,'data','MIR_QCL_DCS'));
sn = d{1}.A;
dt = d{1}.dt;
N  = length(sn);
[fs,Fs] = fftshiftfreqs(sn,dt);

frep_guess = 3.31e6;                         % approximate repetition-rate difference (Hz)

%% Down-select the comb band and resample (saves time/memory)
dsr = [150e6, 1200e6];                        % frequency band containing the comb (Hz)
g   = find(fs>dsr(1) & fs<dsr(2));
dsf = N/length(g);                            % effective downsampling factor
Fs_ds = Fs(g);
N_ds  = length(g);
sn_ds = ifftshiftfreqs(Fs_ds);
dt_ds = dt * dsf;

%% Setup for solver
data = struct();
data.sn = sn_ds;
data.dt = dt_ds;

param = AK_Params();
param.fD_range = frep_guess*[0.75,1.25];

% Tuned process / measurement noise for this dataset (from a prior EM run).
% Start EM from these so it converges in a few iterations.
param.Q = [0.168 -0.00137 0.000213 -2.34e-06;-0.00137 1.24e-05 3.71e-05 -2.53e-07;0.000213 3.71e-05 0.00592 -4.28e-05;-2.34e-06 -2.53e-07 -4.28e-05 3.32e-07];
param.excess_noise = 0.875*2;

% One-step EKF post-regularization + EKF-based EM (best for real data).
param.post_regularization = 'ekf';
param.EM = 'ekf';

%% Run the filter with EM
if ~exist('nEM','var'), nEM = 10; end
[Q, R, oK, em] = Augmented_Kalman_EM(data, param, nEM);

fprintf('MIR_QCL_DCS: learned excess_noise = %.4g\n', R);
