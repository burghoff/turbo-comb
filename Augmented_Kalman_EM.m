function [Q, R, oK, em] = Augmented_Kalman_EM(data, param, nIters, extrapN)
%AK_RunEM  Standalone EM loop wrapper for Augmented_Kalman.
%
%   [Q,R,oK,em] = AK_RunEM(data,param,nIters)
%   [Q,R,oK,em] = AK_RunEM(data,param,nIters,extrapN)
%
%   extrapN (optional) – every extrapN iterations, fit a decaying
%       exponential  a + b*exp(-c*x)  to the last floor(extrapN/2) EM
%       estimates and extrapolate to x→∞ (i.e. use 'a').
%       • excess_noise: extrapolated value used directly.
%       • Q: every upper-triangular element extrapolated independently,
%         then the full symmetric matrix is reconstructed and negative
%         eigenvalues are zeroed to guarantee positive semi-definiteness.
%
% Outputs:
%   Q  = final oK.Q_EM
%   R  = final oK.excess_noise_EM (scalar)
%   oK = final Augmented_Kalman output struct (post-EM)
%   em = history matrix: [Q_upper_tri_elements  excess_noise  err]

if nargin < 4
    extrapN = [];          % disabled by default
end

% ---- upper-triangle index setup (used throughout) ----
nQ     = size(param.Q, 1);
triIdx = find(triu(ones(nQ)));   % linear indices into Q for upper triangle
nQelem = numel(triIdx);          % e.g. 10 for a 4×4 Q
[ri, ci] = ind2sub([nQ, nQ], triIdx);
diagCols = find(ri == ci)';      % positions in em of Q diagonal elements
plotCols = [diagCols, nQelem+1]; % 4 diag cols + excess_noise col

% ---- run once before EM (for "before" plot) ----
param0 = param; % keep initial
oK0 = Augmented_Kalman(data, param0);

% ---- plot correction BEFORE EM ----
dfigure;
a1 = subplot(2,1,1); semilogy(oK0.fss/1e6,oK0.Psn);
xyt('Frequency (MHz)','Power (1)','Uncorrected')
a2 = subplot(2,1,2); semilogy(oK0.fssDc/1e6,oK0.Psn_fDc);
xyt('Frequency (MHz)','Power (1)','Corrected (pre-EM)');
linkaxes([a1,a2]);
top5 = find(oK0.fssDc>.95*(max(oK0.fssDc)-min(oK0.fssDc)+min(oK0.fssDc)));
nse0 = mean(oK0.Psn_fDc(top5));
maxv0 = max(oK0.Psn_fDc);
ylim(a1,10.^[log10(nse0)-.1*(log10(maxv0)-log10(nse0)),log10(maxv0)+.1*(log10(maxv0)-log10(nse0))]);
xlim(a1,[min(oK0.fssDc/1e6),max(oK0.fssDc/1e6)]);

if nIters == 0
    [Q, R, ~, em] = deal(0);
    oK = oK0;
    return;
end

% ---- init EM history ----
% Columns: [Q(triIdx).'  excess_noise  err]
% All upper-triangular elements of Q stored (not sqrt) so off-diagonal
% signs are preserved.
em = [param0.Q(triIdx).', param0.excess_noise, Inf];

% ---- EM iters figure ----
fh = dfigure('DName','EM iters');
ax1 = subplot(1,2,1,'Parent',fh); set(ax1,'YScale','log');
xyt(ax1,'Iteration','Value','EM values');
ax2 = subplot(1,2,2,'Parent',fh); set(ax2,'YScale','log');
xyt(ax2,'Iteration','Value','Diffs');
labels = ["p0","pD","f0","fD",'n'];

if ~isfield(param,'EM') || isequal(param.EM,0) || isequal(param.EM,false)
    param.EM = 1;
end

% start EM from the already-run filter output
oK = oK0;
for iter = 1:nIters

    % update params from last EM estimate
    param.Q            = oK.Q_EM;
    param.excess_noise = max(0, oK.excess_noise_EM);

    % record: full upper-triangle of Q + noise + err
    em(end+1,:) = [param.Q(triIdx).', param.excess_noise, oK.err];
    format shortG
    disp([param.Q(triIdx).', param.excess_noise, oK.err])

    % ---- exponential extrapolation every extrapN iterations ----
    if ~isempty(extrapN) && extrapN > 0 && mod(iter, extrapN) == 0
        nFit  = floor(extrapN / 2);
        nHist = size(em, 1);
        if nHist >= nFit && nFit >= 3
            fitIdx = (nHist - nFit + 1) : nHist;
            xFit   = fitIdx(:);

            % --- extrapolate excess_noise (column nQelem+1) ---
            noiseCol = em(fitIdx, nQelem + 1);
            aInf = fit_decay_extrap(xFit, noiseCol);
            if isfinite(aInf) && aInf > 0
                param.excess_noise = aInf;
                fprintf('  [extrap] excess_noise → %.6g\n', aInf);
            end

            % --- extrapolate every upper-triangular Q element ---
            Qvec = zeros(nQelem, 1);
            for qi = 1:nQelem
                col_i = em(fitIdx, qi);
                aQ    = fit_decay_extrap(xFit, col_i);
                if isfinite(aQ)
                    Qvec(qi) = aQ;
                else
                    Qvec(qi) = col_i(end);   % fallback: last
                end
            end

            % reconstruct full symmetric Q from upper-triangle vector
            Qnew         = zeros(nQ);
            Qnew(triIdx) = Qvec;
            Qnew         = Qnew + Qnew' - diag(diag(Qnew));  % mirror to lower

            % ensure positive semi-definiteness: zero negative eigenvalues
            [V, D] = eig(Qnew);
            d      = diag(D);
            d(d < 0) = 0;
            Qnew = V * diag(d) * V';
            Qnew = (Qnew + Qnew') / 2;   % enforce exact symmetry

            param.Q = Qnew;
            fprintf('  [extrap] Q →\n'); disp(Qnew)
        end
    end

    % rerun filter
    oK = Augmented_Kalman(data, param);

    % plot EM traces (diagonal Q elements + excess_noise)
    for kk = 1:5
        semilogy(ax1, abs(em(:, plotCols(kk))));
        text(ax1, ax1.XLim(2)+.02*(diff(ax1.XLim)), abs(em(end,plotCols(kk))), labels(kk), ...
            'VerticalAlignment','middle', 'FontSize',10,'Color',dColor(kk));
        y = abs(diff(em(:, plotCols(kk))));
        semilogy(ax2, y);
        text(ax2, ax2.XLim(2)+.02*(diff(ax2.XLim)), y(end), labels(kk), ...
            'VerticalAlignment','middle', 'FontSize',10,'Color',dColor(kk));
        if kk==1
            hold(ax1,'all'); hold(ax2,'all');
        elseif kk==5
            hold(ax1,'off'); hold(ax2,'off');
        end
    end
    drawnow limitrate nocallbacks

    % ---- copy final params to clipboard ----
    clipboard('copy', ['param.Q = ', mat2str(oK.Q_EM,3), ';', newline, ...
        'param.excess_noise = ', mat2str(oK.excess_noise_EM,3), ';']);
end

% final outputs
Q = oK.Q_EM;
R = oK.excess_noise_EM;

% ---- plot correction AFTER EM ----
dfigure;
b1 = subplot(2,1,1); semilogy(oK.fss/1e6,oK.Psn);
xyt('Frequency (MHz)','Power (1)','Uncorrected');
b2 = subplot(2,1,2); semilogy(oK.fssDc/1e6,oK.Psn_fDc);
xyt('Frequency (MHz)','Power (1)','Corrected (post-EM)');
linkaxes([b1,b2,a1,a2]);
top5 = find(oK.fssDc>.95*(max(oK.fssDc)-min(oK.fssDc)+min(oK.fssDc)));
nse  = min(mean(oK.Psn_fDc(top5)), nse0);
maxv = max(max(oK.Psn_fDc), maxv0);
ylim(b1,10.^[log10(nse)-.1*(log10(maxv)-log10(nse)),log10(maxv)+.1*(log10(maxv)-log10(nse))]);
xlim(b1,[min(oK.fssDc/1e6),max(oK.fssDc/1e6)]);
end


% =====================================================================
function aInf = fit_decay_extrap(x, y)
%FIT_DECAY_EXTRAP  Fit  y = a + b*exp(-c*x)  and return a  (x→∞ limit).
%
%   Uses a robust two-stage approach:
%     1. Linearised least-squares to get a good initial guess.
%     2. Nonlinear least-squares (lsqcurvefit / fminsearch) to refine.
%
%   Returns NaN if the fit fails or produces non-physical results.

x = x(:);  y = y(:);
x = x - mean(x);
aInf = NaN;
n = numel(y);
if n < 3, return; end

% --- Stage 1: linearised initialisation ---
dy = diff(y);
if all(dy >= 0)
    a0 = y(end) + 0.05 * abs(y(end) - y(1));
elseif all(dy <= 0)
    a0 = y(end) - 0.05 * abs(y(end) - y(1));
else
    a0 = y(end);
end

z = y - a0;
sgn = sign(z(1));
if sgn == 0, sgn = 1; end
z = sgn * z;
z(z <= 0) = min(z(z>0))*0.01;
logz = log(z);
P  = [x, ones(n,1)] \ logz;
c0 = -P(1);
b0 = sgn * exp(P(2));
if c0 <= 0, c0 = 0.1; end

% --- Stage 2: nonlinear refinement ---
expModel = @(p,xd) p(1) + p(2)*exp(-p(3)*xd);
p0 = [a0, b0, c0];

try
    if exist('lsqcurvefit','file')
        opts = optimoptions('lsqcurvefit','Display','off', ...
            'MaxIterations',500,'MaxFunctionEvaluations',2000, ...
            'FunctionTolerance',1e-14,'StepTolerance',1e-14);
        lb = [-Inf, -Inf, 1e-8];
        ub = [ Inf,  Inf, Inf ];
        pFit = lsqcurvefit(expModel, p0, x, y, lb, ub, opts);
    else
        costFn = @(p) sum((expModel(p,x) - y).^2);
        opts = optimset('Display','off','MaxIter',1000,'TolFun',1e-14,'TolX',1e-14);
        pFit = fminsearch(costFn, p0, opts);
        if pFit(3) <= 0, return; end
    end
catch
    return;
end

aInf = pFit(1);
end