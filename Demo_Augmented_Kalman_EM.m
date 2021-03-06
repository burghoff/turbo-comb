addpath('pathfcns');
Generate_Data;

%% Setup for solver
data=struct();
data.sn = Vn;                            % input signal
data.dt = dt;                            % input timestep

% Required parameters
param=struct();
param.fD_range = fDm*[0.75,1.25];        % Frequency range that the rep rate is known to be in.
                                         % If you don't know, get from abs(fft(Vn)).^2.
param.Q = diag([1,.01,1,.01].^2);        % Process noise matrix (4x4)
                                         % If you don't know, set to diag([1,.01,1,.01].^2) and run EM
param.excess_noise = 1;                  % Excess measurement noise factor
                                         % Code will calculate PSD near the Nyquist frequency and multiply by this to determine the noise.
                                         % If you don't know, set to 1 and run EM

% Optional parameters
param.knownfD = fD;                      % User-supplied repetition rate. When supplied, is plotted alongside estimate (NOT used otherwise).
param.knownf0 = f0;                      % User-supplied offset. When supplied, is plotted alongside estimate (NOT used otherwise).

[sqrt(diag(param.Q)).',param.excess_noise,Inf]
oK=Augmented_Kalman(data,param);         % Run the filter!
Efficacy_Plots

%% Expectation maximization
em = [sqrt(diag(param.Q)).',param.excess_noise,Inf];
for iter=1:10
    param.Q = oK.Q_EM;
    param.excess_noise = max(0,oK.excess_noise_EM);
    [sqrt(diag(param.Q)).',param.excess_noise,oK.err]
    em(end+1,:) = [sqrt(diag(param.Q)).',param.excess_noise,oK.err];
    oK=Augmented_Kalman(data,param);
    dfigure('DName','EM iters');
    for ii=1:5
        subplot(1,2,1); semilogy(abs(em(:,ii))); hold all;
        xyt('Iteration','Value','EM values');
        subplot(1,2,2); semilogy(abs(diff(em(:,ii)))); hold all;
        xyt('Iteration','Value','Diffs');
    end
end
clipboard('copy',['param.Q = ',mat2str(oK.Q_EM,3),';',newline,...
                  'param.excess_noise = ',mat2str(oK.excess_noise_EM,3),';']) % copy result to clipboard              
Efficacy_Plots