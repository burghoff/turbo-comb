function [Vf, info] = Rubber_Band_Filter_Inline(ts, Vs, freq, varargin)
% Rubber_Band_Filter
% Low-pass filtering without edge artifacts via optimal padding
% solved by Conjugate Gradient (CG).

% examples:
%   Vf = Rubber_Band_Filter(ts, Vs, freq)
%   [Vf,info] = Rubber_Band_Filter(ts, Vs, freq)
%   Vf = Rubber_Band_Filter(ts, Vs, freq, 'exitmode','fixed')
%   Vf = Rubber_Band_Filter(ts, Vs, freq, 'exitmode','thresh','tol',1e-9)
%   Vf = Rubber_Band_Filter(ts, Vs, freq, 'Niter',10)

%
% Required:
%   ts   : time vector (monotonic, ~uniform)
%   Vs   : signal vector (same length as ts)
%   freq : cutoff frequency
%
% Exit modes:
%   'thresh' (default): break when real(r'*r)/real(b'*b) <= tol  (tol default 1e-5)
%   'fixed'           : run exactly Niter iterations

% Defaults implement right-side mirror padding with padded weight
%   0.03 * tukeywin(L, 0.2), and valid region weight = 1 if sample is valid.

% Optional parameters:
%   'exitmode': fixed or threshold
%   'tol'     : none for fixed mode or tolerance for threshold mode
%   'npad'    : number of extra padding regions (default = 1; right-side mirror)
%   'niter'   : number of CG iterations (default = 30)
%   'ws'      : weights for valid region (default = 1 except 0 where Vs is NaN)


    ts = ts(:); Vs = Vs(:);
    assert(numel(ts) == numel(Vs), 'ts and Vs must have same length.');

    % ---- defaults ----
    exitMode = 'thresh';
    tol      = 1e-5;
    Npad     = 1;
    Niter    = 30;
    ws       = ones(size(Vs));

    % case-by-case varargin scan (key at i, value at i+1)
    for i1 = 1:2:length(varargin)
        key = varargin{i1};
        if ischar(key) || (isstring(key) && isscalar(key))
            key = lower(strtrim(char(key)));
            if i1+1 > length(varargin), error('Option "%s" missing value.', key); end
            val = varargin{i1+1};
            switch key
                case 'exitmode', exitMode = canonical_exitmode(val);
                case 'tol',      tol      = val;
                case 'npad',     Npad     = val;
                case 'niter',    Niter    = val;
                case 'ws',       ws       = val;
            end
        end
    end
    if isempty(ws), ws = ones(size(Vs)); else, ws = ws(:); end

    % ---- handle NaNs: deweight and interpolate for processing ----
    validMask = ~isnan(Vs);
    ws = ws(:) .* double(validMask);
    if any(~validMask)
        Vs(~validMask) = interp1(ts(validMask), Vs(validMask), ts(~validMask), 'linear', 'extrap');
    end

    % ---- inline padding ----
    Vso = Vs;                           % keep original for cropping
    if Npad == 0
        % no padding
        wk = ws;
    elseif Npad == 1
        % right-side mirror with padded weight = 0.03 * tukeywin(L, 0.2)
        windowfcn = tukeywin(length(ws), 0.2);
        wsh = 0.03 .* windowfcn;
        Vs  = [Vso; flipud(Vso)];
        wk  = [ws;  wsh];
    else
        % original extended pattern
        hn = hanning(2*length(Vso)); hn = hn(1:length(Vso));
        Vs = [Vso; flipud(Vso.*hn); zeros(length(Vso)*(Npad-2),1); flipud(Vso).*hn];
        wk = [ws; zeros(length(ws)*Npad,1)];
    end

    % ---- frequency grid & passband (odd length) ----
    dt  = mean(diff(ts));
    N   = length(Vs);
    fss = fftshift(1/(N*dt) * (0:N-1)');        % centered after shift
    zi  = find(fss == 0, 1, 'first');
    if ~isempty(zi), fss(1:zi-1) = fss(1:zi-1) - 1/dt; else, error('Zero frequency not found.'); end

    g    = find(abs(fss) <= freq);
    xtra = 0;
    while mod(length(g), 2) == 0                % ensure odd passband length
        xtra = xtra + eps(freq);
        g    = find(abs(fss) <= (freq + xtra));
    end

    % ---- precompute transforms ----
    W = fftshift(fft(wk));
    P = fftshift(fft(Vs));

    % conv over filtered region (nested; uses zi, fss, g, W)
    function co = convf(xi)
        if length(xi) == length(g), zix = find(fss(g) == 0, 1, 'first');
        else,                       zix = find(fss     == 0, 1, 'first');
        end
        co_full = ifft( fft(xi, length(xi)+length(W)-1) .* fft(W, length(xi)+length(W)-1) );
        ko = zi - 1 + zix;
        co = co_full( ko - (length(g)-1)/2 : ko + (length(g)-1)/2 );
    end

    % ---- CG ----
    b = convf(P);
    x0 = P(g); x = x0;                      %#ok<NASGU> preserve intent
    % rng('default'); 
    x0 = zeros(size(P(g))); x = x0;
    r = b - convf(x);
    p = r;
    b_t=ToV(b);
    bnorm2=std(b_t(1:length(Vso)));

    dv = zeros(1, Niter);                   % normalized residual^2
    for ii = 1:Niter
        Ap = convf(p);
        ak = (r' * r) / (p' * Ap);
        xo = x;
        x  = x + ak * p;
        ro = r;
        r  = r - ak * Ap;
        bk = (r' * r) / (ro' * ro);
        % bk = max(0, r'*(r-ro)/(ro'*ro));
        if sum(wk .* abs((ToV(x) - Vs)).^2) > sum(wk .* abs((ToV(xo) - Vs)).^2)
            bk = 0;
        end
        % bk = 0;
        p = r + bk * p;
        Ax = convf(x); Ax_t=ToV(Ax);
        dv(ii) = max(abs(Ax_t(1:length(Vso))-b_t(1:length(Vso)))) / bnorm2;

        if strcmpi(exitMode, 'thresh')
            if dv(ii) <= tol
                dv = dv(1:ii);
                break;
            end
        end
        % 'fixed' -> never break early
    end
    used_iter = min(ii, Niter);

    % ---- outputs ----
    Vf_pad = ToV(x);
    if isreal(Vs), Vf_pad = real(Vf_pad); end
    Vf = Vf_pad(1:length(Vso));           % crop back to original length
    converged = true;                 % assume OK unless shown otherwise
    final_dv  = dv(end);

    if strcmpi(exitMode,'thresh')
    % If we did NOT break early (i.e., hit Niter) and dv still > tol → warn
         if used_iter == Niter && final_dv > tol
             converged = false;
            %  warning('Rubber_Band_Filter:ThresholdNotReached', ...
            % 'CG hit Niter=%d before meeting tol=%.3g (final dv=%.3g). Consider increasing Niter or relaxing tol.', ...
            % Niter, tol, final_dv);
        end
    end
    if nargout > 1
        info = struct( ...
            'dt', dt, ...
            'N_original', length(Vso), ...
            'N_padded',   length(Vs), ...
            'cutoff_used', freq + xtra, ...
            'passband_len', length(g), ...
            'iterations_requested', Niter, ...
            'iterations_used', used_iter, ...
            'exitMode', exitMode, ...
            'tol', tol, ...
            'dv', dv, ...
            'final_dv', final_dv, ...      % <- added
            'converged', converged, ...    % <- added
            'err', sum(wk .* abs((ToV(x) - Vs)).^2), ...
            'Vf_pad', Vf_pad ...
        );
    end

    % ---- nested: spectrum->time ----
    function Vout = ToV(xin)
        VoutSpec = zeros(size(fss));
        VoutSpec(g) = xin;
        Vout = ifft(ifftshift(VoutSpec));
    end
end

function s = canonical_exitmode(x)
    s = lower(string(x));
    if any(s == ["fixed","niter","maxiter","full"])
        s = "fixed";
    else
        s = "thresh";
    end
    s = char(s);
end