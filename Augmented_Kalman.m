function output = Augmented_Kalman(data,param)
%% AUGMENTED_KALMAN: runs the generalized correction procedure in Burghoff et al., Optics Letters (2019).
% If you found this useful in your research, then please cite our work!
% Requirements:
% 1. Newer versions of MATLAB with implicit expansion
% 2. XSum: sum with error compensation (https://www.mathworks.com/matlabcentral/fileexchange/26800-xsum)
%
% Inputs:
%   data: A structure that contains two fields:
%       sn: a complex signal represented by a column vector
%       dt: the timestep
%   param: A structure that configures the parameters of the correction
%       Required parameters:
%       fD_range:     A two-element row vector that determines the frequency range over which the 
%                     repetition rate is known to occur in. Most easily determined by examining the FFT of |sn|^2.
%                     For now, we do not have an easy robust way of automating this. (Suggestions welcome!)
%       Q:            4x4 process noise matrix, arranged in order (p0,pD,f0,fD)
%                     If you don't know, set to diag([1,.01,1,.01].^2) and run EM.
%       excess_noise: Excess noise factor. (How much larger the noise is than the top 10% of the frequency space.)
%                     If you don't know, set to 1 and run EM.
%
%       Optional parameters:
%       plotme:     Whether or not to plot the current status of the correction (default: 1)
%       knownns:    When set, the code will only use the comb lines you specify (default: fill the Nyquist space)
%       knownA:     When set, amplitudes will be set to these values instead of determined by demodulation
%                   Useful for debugging simulated data. (default: demodulate instead)
%       knownfD:    User-supplied repetition rate. When supplied, is plotted alongside estimate (NOT used otherwise).
%       knownf0:    User-supplied offset. When supplied, is plotted alongside estimate (NOT used otherwise).
%       pulsed_mode_apodization: When data is in pulsed mode, enabling this triggers an apodization procedure
%                                that reduces the computation time. Determines the fraction of the data that is
%                                kept (i.e., 0.1 means that only 10% of the data near each burst is used).
%                                Not fully tested, so use with some caution. 
%       initfrac:   When the processing starts, we do not know the amplitudes and need to bootstrap them. Set aside
%                   this fraction of the data to pre-process, so we can determine an initial estimate (default=0.1).
%       Ninits:     How many times to run the initialization procedure (default=2)
%       global_search_stds:      How large of a range to cover when running the global search, in terms of the std
%                                deviations of the current estimate. (default=6, i.e. six sigma)
%       global_search_maxsize:   Can be used to limit the number of elements in the global search tensor. (default=1e6) 
%       EM:                      Whether or not to perform expectation maximization (default=1) 
%       post_regularization:     Post-regularization procedure to use (default='lpf', finds the band-limited signal
%                                that best fits the data)
%       
%% Get data from structures
% To avoid precision issues, scale time to be in units of the approximate rep rate
param = Initialize_Params(param);
dt= data.dt; sn= data.sn(:,1); N=length(data.sn);

Ts = 1/mean(param.fD_range);
dt = dt/Ts;
param.fD_range = param.fD_range*Ts;
ts=dt*[0:N-1]';

%% Figure out measurement noise standard deviation from white noise level, then multiply by param.excess_noise
fs = 1/(N*dt)*[0:N-1]';
fss = fftshift(fs); fss(1:find(diff(fss)<0))=fss(1:find(diff(fss)<0))-1/dt;     % fftshifted-fs
Fsn = fft(sn);

psd  = fftshift(abs(Fsn)).^2*dt/N;
npsd = mean(psd(abs(fss)>=1/dt/2*0.9));
nstd = sqrt((npsd/2)*(1/dt));                                                   % noise in each quadrature
nstd = nstd*param.excess_noise;
% figure; semilogy(fss,psd);                                  % for spying
% hold all; p=semilogy(fss,fftshift(abs(fft(nstd*(randn(N,1)+i*randn(N,1)))).^2*dt/N)); uistack(p,'bottom');

[~,ml]=max(psd);
maxf = fss(ml); stdf = sqrt(sum(psd.*(fss-maxf).^2)/sum(psd));
cdf = cumsum(psd)/sum(psd); [~,ml]=min(abs(cdf-0.5)); medf=fss(ml);
%% Kalman filter structure and initial values
fDhi = param.fD_range(2); fDlo = param.fD_range(1);

% Initial values
fD0 = (fDlo+fDhi)/2;    
pD0 = 0;                   % phases                         
f00 = medf;0;
p00 = 0;

givenn = ~isnan(param.knownns);
if ~givenn, n = [round((-1/dt-f00)/fDlo/2):round((1/dt-f00)/fDlo/2)];
else,       n = reshape(param.knownns,1,length(param.knownns));
end
Nl = length(n);
num_stds = param.global_search_stds;

% Initial covariances
PfD  = ((fDhi-fDlo)/2)^2;           % rep rate
Pf0  = (stdf)^2;                    % rep rate
Pp0 =  (2*pi)^2*10;                 % defined at t=0
PpD =  (2*pi)^2*10;                 % defined at t=0

%% Divide data into batches
% In pulsed mode, signal is sparse so apodizing can save a lot of time /
% memory. This is somewhat experimental and untested, though.
M0 = ceil(1/((fDlo+fDhi)/2*dt));
if ~isnan(param.pulsed_mode_apodization)
    [pks,locs]=findpeaks(abs(sn).^2,'MinPeakDistance',0.75*M0);
    M = ceil(param.pulsed_mode_apodization/((fDlo+fDhi)/2*dt));
    if mod(M,2)==1, Mr = [-(M-1)/2,(M-1)/2];
    else,           Mr = [-M/2+1,M/2]; end
    starts = locs(2:end-1)'+Mr(1); stops = locs(2:end-1)'+Mr(2);    % starts and stops of each batch processed by filter
    fstops = [starts(2:end),locs(end)'+Mr(1)]-1;                    % stops of the full batches that include data in between
else
    M = M0;
    starts = [1:M:length(sn)]; stops = [M:M:length(sn)]; starts=starts(1:length(stops));
    fstops = stops;
end

Nb = length(starts);
bl = cell2mat(arrayfun(@(x,y)[x:y]',starts,stops,'un',0));   % batch locations                                             
snr = sn(bl); tsf = ts(bl);

Md =[NaN,diff(starts)];                                      % distance from a start to the previous one
lgf = fstops-starts+1;                                       % length of the kth full batch

%% Create and run Kalman filter
% State structure
lp0=1; lpD=2; lf0 = 3; lfD = 4;
x0 = [p00;pD0;f00;fD0];
P0 = diag([Pp0;PpD;Pf0;PfD]);

% Process function
F1 = eye(length(x0));                   % single timestep
F1(lpD,lfD)=F1(lpD,lfD)+2*pi*dt;
F1(lp0,lf0)=F1(lp0,lf0)+2*pi*dt;

Fs = arrayfun(@(x)F1^x,Md,'un',0);
Q = param.Q;
Qs = arrayfun(@(x)Q/M0*x,Md,'un',0);

% Nonlinear measurement function
tk = dt*[0:M-1]';
function [ho,Ho] = hH(x,GIi,GIni)
%   Method 1: Explicit summation, slower for large Nl (takes A as argument)
%     expv = A.'.*exp(i*( x(lp0)+2*pi*x(lf0)*tk + n.*(x(lpD)+2*pi*x(lfD)*tk) ));
%     hc = sum(expv,2); ho = [real(hc);imag(hc)]; 
%     if nargout>1
%         Hc = i*[hc,sum(expv.*n,2),2*pi*sum(expv.*tk,2),2*pi*sum(expv.*tk.*n,2)];
%         Ho = [real(Hc);imag(Hc)];
%     end
%   Method 2: Interpolate over FFT, faster for large Nl
    tp = (x(lpD) + 2*pi*x(lfD)*tk)/(2*pi); tp=tp-floor(tp);
    exp0 = exp(i*(x(lp0)+2*pi*x(lf0)*tk));
    hc = exp0.*GIi(tp);
    ho = [real(hc);imag(hc)]; 
    if nargout>1
        GInv = GIni(tp);
        Hc = i*[hc,exp0.*GInv    ,2*pi*hc.*tk         ,2*pi*exp0.*GInv.*tk];
        Ho = [real(Hc);imag(Hc)];
    end
end
R = eye(M*2)*nstd^2;            % noise is sigma^2 in each quadrature

function [ho,Ho] = hHf(x,GIi,GIni,Min)
    tkf = dt*[0:Min-1]';
    tp = (x(lpD) + 2*pi*x(lfD)*tkf)/(2*pi); tp=tp-floor(tp);
    exp0 = exp(i*(x(lp0)+2*pi*x(lf0)*tkf));
    hc = exp0.*GIi(tp);
    ho = [real(hc);imag(hc)]; 
    if nargout>1
        GInv = GIni(tp);
        Hc = i*[hc,exp0.*GInv    ,2*pi*hc.*tkf         ,2*pi*exp0.*GInv.*tkf];
        Ho = [real(Hc);imag(Hc)];
    end
end


% Setup an interpolant for the interferogram p(t)
ZPF = 5;                                        % zero padding factors
fsi = [-ZPF*max(abs(n)):ZPF*max(abs(n))]';
imfn = find(ismember(fsi,n));
Fsi = zeros(size(fsi)); Fsi(imfn)=zeros(size(n));
VF = length(Fsi)*ifft(ifftshift(Fsi)); VF=[VF;VF(1)];
tf = 1/(length(Fsi)*1)*[0:length(Fsi)]';
[GI,GIn] = deal(griddedInterpolant(tf,VF,'cubic'));    % change Values in loop instead of regenerating, faster than interp1
function Vo = Amplitude_to_IFG(Ai)                     % calculates p(t) from amplitudes using FFT
    Fsi(imfn)=Ai;
    Vo = length(Fsi)*ifft(ifftshift(Fsi)); Vo=[Vo;Vo(1)];
end

% Setup a demodulation function that uses FFT instead of calculating an
% exponential explicitly, which is slow when Nl is large
function yd = FFT_Demod(x,y)
    fd = 1/(ZPF*length(y)*dt)*[0:ZPF*length(y)-1]'; fds = fftshift(fd); zi=find(fds==0); fds(1:zi-1)=fds(1:zi-1)-1/dt;
    GIf = griddedInterpolant(fds,fds,'cubic','none');
    GIf.Values = fftshift(fft(y,ZPF*length(y)));
    Fy = GIf(x(lf0)+n*x(lfD)); Fy(isnan(Fy))=0;
    yd = exp(-i*(x(lp0)+n*x(lpD))) .* Fy / length(y);
end


% tk1 = dt*[0:fstops(1)-starts(1)+1]';
% function expv = expn(x)
%     expv = exp(-i*( x(lp0)+2*pi*x(lf0)*tk1 + n.*(x(lpD)+2*pi*x(lfD)*tk1) ));
% end

% Other initialization steps
pj1 = [0:-2*pi*fDlo*dt:-pi]; pj2 = [0:2*pi*fDlo*dt:pi];
pj0 = [fliplr(pj1),pj2(2:end)];         % densest phi_r grid we might need

givenA  = ~isnan(param.knownA);
givenf0 = ~isnan(param.knownf0);
givenfD = ~isnan(param.knownfD);
if ~givenA
    An = ones(Nl,Nb)*NaN;
%     An(:,1) = mean(expn(x0).*sn(starts(1):fstops(1)),1).';
%     An(:,1) = FFT_Demod(x0,[snr(:,1);zeros(lgf(1)-M,1)]).';
    An(:,1) = FFT_Demod(x0,sn(starts(1):fstops(1))).';
end
initk = floor(param.initfrac*Nb);
initsleft = param.Ninits;


%% Main loop
tic;
[xk,dPkk] = deal(zeros(Nb,length(x0)));
[xkkm1s,Ck,Pkkm1s,Pkks]= deal(cell(Nb,1));
pkk = zeros(M,Nb);
warning('off','MATLAB:nearlySingularMatrix'); 
if param.plotme, dfigure('DPosition',[1,1],'DName','Augmented Kalman'); subplot(2,2,1); end
k=1;
while k<=Nb
    if k==1
        xkkm1 = x0;
        Pkkm1 = P0;
    else
        xkkm1 = Fs{k}*xkk;
        Pkkm1 = Fs{k}*Pkk*Fs{k}' + Qs{k};
        Pkkm1=(Pkkm1+Pkkm1')/2;
    end
        
    if k>1
        Ck{k} = Pkk*Fs{k}'/Pkkm1;
        xkkm1s{k} = xkkm1;
    end
    if param.EM==1
        Pkkm1s{k} = Pkkm1;
    end
    
    meas = snr(:,k);
    fmeas = sn(starts(k):fstops(k));
    if givenA
        Av = param.knownA;
    else
        nn = find(~isnan(An(1,:)));
        Av = sum(An(:,nn).*lgf(nn),2)./sum(lgf(nn));        % weighted sum
        Av(abs(xkkm1(lf0)+n*xkkm1(lfD))>1/dt) = 0;          % Ensure any lines outside the Nyquist band are 0
    end
    
    GI.Values  = Amplitude_to_IFG(Av);                      % Prepare interpolaters for these amplitudes
    GIn.Values = Amplitude_to_IFG(Av.*n');
    
    %% Find candidate local minima of MAP using a grid search
    % Rows: f_0   (offset freq)
    % Cols: phi_D (rep rate phase)
    % 3s:   f_Ds  (rep rate freqs)
    pjs = xkkm1(lpD) + pj0;
    pj =  pjs  (pjs  >=xkkm1(lpD)-num_stds*sqrt(Pkkm1(lpD,lpD)) & pjs  <=xkkm1(lpD)+num_stds*sqrt(Pkkm1(lpD,lpD)));
    
    xkkm1(lfD) = median([fDlo+eps,xkkm1(lfD),fDhi-eps]);
%     fD1 = [1/xkkm1(lfD)/dt:1/fDlo/dt].'; fD2 = [1/xkkm1(lfD)/dt:-1:1/fDhi/dt].';
%     fDs0 = 1./[flipud(fD1(2:end));fD2]./dt; fDs0 = reshape(fDs0,1,1,length(fDs0));
    fD1 = [xkkm1(lfD):-1/M*xkkm1(lfD)/2:fDlo]'; fD2 = [xkkm1(lfD):1/M*xkkm1(lfD)/2:fDhi]';
    fDs0 = [flipud(fD1(2:end));fD2]; fDs0 = reshape(fDs0,1,1,length(fDs0));
    fDs = fDs0(fDs0>=xkkm1(lfD)-num_stds*sqrt(Pkkm1(lfD,lfD)) & fDs0<=xkkm1(lfD)+num_stds*sqrt(Pkkm1(lfD,lfD)));
    
    f0ss = 1/(ZPF*M*dt)*[0:ZPF*M-1]'; f0ss = fftshift(f0ss); zi=find(f0ss==0); f0ss(1:zi-1)=f0ss(1:zi-1)-1/dt;
    f0ss = xkkm1(lf0)-f0ss;       % f0 sample frequencies
    gf0 = find((f0ss>=xkkm1(lf0)-num_stds*sqrt(Pkkm1(lf0,lf0)) & f0ss<=xkkm1(lf0)+num_stds*sqrt(Pkkm1(lf0,lf0))));
    f0s =  f0ss(gf0);
    
    % To prevent huge global searches, reduce frequency search when Flj would be too large
    % May want to fiddle with this for large problems
    numsearch = length(f0s)*length(pj)*length(fDs); nss = num_stds; 
    while numsearch > param.global_search_maxsize          
        nss = nss/sqrt(numsearch/param.global_search_maxsize);
        fDs = fDs0( fDs0>=xkkm1(lfD) -nss*sqrt(Pkkm1(lfD,lfD))  & fDs0<=xkkm1(lfD)+nss*sqrt(Pkkm1(lfD,lfD)));
        gf0 = find((f0ss>=xkkm1(lf0) -nss*sqrt(Pkkm1(lf0,lf0))  & f0ss<=xkkm1(lf0)+nss*sqrt(Pkkm1(lf0,lf0))));
        f0s =  f0ss(gf0);
%         pj =  pjs  (pjs  >=xkkm1(lpD)-nss*sqrt(Pkkm1(lpD,lpD)) & pjs  <=xkkm1(lpD)+nss*sqrt(Pkkm1(lpD,lpD)));
        numsearch = length(f0s)*length(pj)*length(fDs);
    end
    
%     k
      
    tp = (pj+ 2*pi*fDs.*tk)/(2*pi); tp=tp-floor(tp); % numel(tp)
    Pkj = GI(tp);
    Flj = fftshift(fft(conj(meas.*exp(-i*2*pi*xkkm1(lf0)*tk)).*Pkj,ZPF*M),1);
    Flj = Flj(gf0,:,:); 

    C = sum(abs(meas).^2);
    Slj = repmat(sum(abs(Pkj).^2,1),size(Flj,1),1,1);
    
    % Find optimal phi_0 for each point using Newton's method
    Pi = Pkkm1^-1;
    C2 = Pi(lp0,lpD)*(pj-xkkm1(lpD)) + Pi(lp0,lf0)*(f0s-xkkm1(lf0)) + Pi(lp0,lfD)*(xkkm1(lfD)-xkkm1(lfD));
    M00 = Pi(lp0,lp0);
    aFlj = angle(Flj); mFlj = abs(Flj);
    pn = 2*pi*round(1/2/pi*(xkkm1(lp0) - C2/M00 + aFlj)) - aFlj;
    for iN=1:3
        pn = pn - (Pi(lp0,lp0)*(pn-xkkm1(lp0))+C2 + 1/nstd^2*mFlj.*sin(pn+aFlj)) ./ ...
                  (M00 + 1/nstd^2*mFlj.*cos(pn+aFlj));
    end
    
    % Calculate J(x), the MAP function
    eP = (pn-xkkm1(lp0))        .*(Pi(lp0,lp0)*(pn-xkkm1(lp0))+Pi(lp0,lpD)*(pj-xkkm1(lpD)) + Pi(lp0,lf0)*(f0s-xkkm1(lf0)) + Pi(lp0,lfD)*(xkkm1(lfD)-xkkm1(lfD))) + ...
         (pj-xkkm1(lpD))        .*(Pi(lpD,lp0)*(pn-xkkm1(lp0))+Pi(lpD,lpD)*(pj-xkkm1(lpD)) + Pi(lpD,lf0)*(f0s-xkkm1(lf0)) + Pi(lpD,lfD)*(xkkm1(lfD)-xkkm1(lfD))) + ...
         (f0s-xkkm1(lf0))       .*(Pi(lf0,lp0)*(pn-xkkm1(lp0))+Pi(lf0,lpD)*(pj-xkkm1(lpD)) + Pi(lf0,lf0)*(f0s-xkkm1(lf0)) + Pi(lf0,lfD)*(xkkm1(lfD)-xkkm1(lfD))) + ...
         (xkkm1(lfD)-xkkm1(lfD)).*(Pi(lfD,lp0)*(pn-xkkm1(lp0))+Pi(lfD,lpD)*(pj-xkkm1(lpD)) + Pi(lfD,lf0)*(f0s-xkkm1(lf0)) + Pi(lfD,lfD)*(xkkm1(lfD)-xkkm1(lfD)));
    eR = (C + Slj - 2*real(Flj.*exp(i*pn)))/nstd^2;
    elj = eR + eP;

    % Pick a few candidate minima from the grid
    Nls = min(15,numel(elj));
    [~,mls]=mink(elj(:),Nls);
    [mrs,mcs,m3s]=ind2sub(size(elj),mls);
    xls = zeros(4,Nls+1);  
    for il=1:Nls
        mf0 = f0s(mrs(il)); mpD = pj(mcs(il)); mfD = fDs(m3s(il));
        mp0 = pn(mrs(il),mcs(il),m3s(il));
        xls(:,il) = [mp0;mpD;mf0;mfD];
    end
    xls(:,end) = xkkm1;

    %% Amongst all candidates, do iterated EKF
    Pi = Pkkm1^-1; 
    Ps = ones(1,size(xls,2))*Inf;
    [Hs,xs,sks] = deal(cell(1,size(xls,2)));
    
    for ic=1:size(xls,2)
        xl = xls(:,ic);         % candidate minima of MAP
        Pl = Inf; keepgoing = 1;
        while keepgoing
            [hv,Hv]=hH(xl,GI,GIn);                          
            Pn = sum(abs(hv-[real(meas);imag(meas)]).^2)/nstd^2 + (xl-xkkm1)'*(Pi*(xl-xkkm1));
            if Pl-Pn < .01*Pl
                keepgoing = 0;
            end
            Pl = Pn;
            yk = [real(meas);imag(meas)]-hv-Hv*(xkkm1-xl);  
            xl = xkkm1 + Pkkm1*((Hv'*Hv*Pkkm1+nstd^2*eye(4))\(Hv'*yk));
            % The above line is just a speed-optimized version of the usual
            % Kalman filter update, with K=P*H'/(H*P*H'+R). (It makes use
            % of the identity H'/(H*P*H'+R)=I/(H'*H*P+r*I)*H', which
            % is valid when R is diagonal and has diag(R)=r*I.
        end
        sks{ic} = hH(xl,GI,GIn)-[real(meas);imag(meas)];
        Ps(ic) = sum(abs(sks{ic}).^2)/nstd^2 + (xl-xkkm1)'*(Pi*(xl-xkkm1));
        Hs{ic} = Hv; xs{ic}=xl;
    end
    [~,ml]=min(Ps);
    Hv=Hs{ml}; xkk=xs{ml}; sk=sks{ml};
    Pkk = (eye(4) - Pkkm1*((Hv'*Hv*Pkkm1+nstd^2*eye(4))\(Hv'*Hv)))*Pkkm1;   % speed optimized
     
    % Following state update, demodulate batch
    if ~givenA
        An(:,k) = FFT_Demod(xkk,fmeas).';
    end
    
    hv=hH(xkkm1,GI,GIn); pkk(:,k)=real(hv(1:M))+i*(hv(M+1:end)); % predicted value
    
    xk(k,:) = reshape(xkk,1,length(x0));
    dPkk(k,:)=reshape(diag(Pkk),1,length(x0));
    Pkk = (Pkk+Pkk')/2;
    if param.EM==1, Pkks{k} = Pkk; end

    % Plots
    if (mod(k,50)==1 | k==Nb) & param.plotme
        subplot(2,2,3); p=plot(mean(tsf(:,xk(:,lfD)~=0),1),xk(xk(:,lfD)~=0,lfD)/Ts);
        xyt('Batch','f_r','Repetition rate');
        if givenfD, hold all; plot(ts,param.knownfD(1:N)); hold off; uistack(p,'top'); end
        subplot(2,2,1); p=plot(mean(tsf(:,xk(:,lf0)~=0),1),xk(xk(:,lf0)~=0,lf0)/Ts);
        if givenf0, hold all; plot(ts,param.knownf0(1:N)); hold off; uistack(p,'top'); end
        xyt('Batch','f_0','Offset');
        subplot(2,2,4); semilogy(n,abs(Av).^2);
        xyt('n','Power','Current power estimate');
        drawnow;
    end
   
    k=k+1;
    % Initializations
    % When we start, we have no idea what the interferogram looks like, so our demodulated estimates are bad.
    % Therefore, once some data has been processed we start over with the knowledge of the interferogram.
    if k>initk & initsleft>0 & ~givenA
        xkni=cell(k-1,1);
        for ki=k-1:-1:1                 % do a mini-RTS when 
            xkki = reshape(xk(ki,:),length(x0),1);
            if ki==k-1, xkni{ki} = xkki;
            else,       xkni{ki} = xkki + Ck{ki+1}*(xkni{ki+1}-xkkm1s{ki+1});
            end
        end
        k=1; initsleft=initsleft-1;
        x0 = xkni{1};
%         x0 = xkk; P0 = Pkk;
    end
end,toc

residpwr = mean(mean(abs(snr-pkk).^2));
sigpower = mean(mean(abs(snr).^2));
err=mean(mean(abs(snr-pkk).^2));  residpwr/sigpower;                                   % fractional error

%% RTS smoothing
xkn = xk;
[xkns,Pkns,Exkxkp1p] = deal(cell(Nb,1));
for k=Nb:-1:1
    xkk = reshape(xk(k,:),length(x0),1);
    if k==Nb
        xkn(k,:) = xkk.';
        xkns{k} = xkn(k,:).';
        Pkns{k} = Pkks{k};
    else    
        xkn(k,:) = (xkk + Ck{k+1}*(xkn(k+1,:).'-xkkm1s{k+1})).'; 
        if param.EM==1
            xkns{k} = xkn(k,:).';
            Pkns{k} = Pkks{k} + Ck{k+1}*(Pkns{k+1}-Pkkm1s{k+1})*Ck{k+1}'; Pkns{k}=(Pkns{k}+Pkns{k}')/2;
            Exkxkp1p{k}=xkk*xkn(k+1,:) + Ck{k+1}*(Pkns{k+1}+(xkn(k+1,:)'-xkkm1s{k+1})*xkn(k+1,:));  % expected x_k * x_{k+1}'
        end
    end
end

%% Resample onto a 1-timestep grid and perform post-regularization
[ws,p0,pD] = deal(zeros(size(sn)));
for k=1:Nb
    p0(starts(k):fstops(k)) = xkn(k,lp0)+2*pi*(ts(starts(k):fstops(k))-ts(starts(k)))*xkn(k,lf0);
    pD(starts(k):fstops(k)) = xkn(k,lpD)+2*pi*(ts(starts(k):fstops(k))-ts(starts(k)))*xkn(k,lfD);
    ws(starts(k):stops(k))=abs(sn(starts(k):stops(k))).^2;
end
vr = [starts(1):stops(end)]; % valid range
ws = ws(vr); p0 = p0(vr); pD = pD(vr);

switch param.post_regularization
    case 'lpf'
    % Low pass filter under the assumption that noise is bandlimited
        fDm = (pD(end)-pD(1))/(ts(vr(end))-ts(vr(1)))/2/pi;
        p0m = (p0(end)-p0(1))/(ts(vr(end))-ts(vr(1)))*(ts(vr)-ts(vr(1)))+p0(1);
        pDm = (pD(end)-pD(1))/(ts(vr(end))-ts(vr(1)))*(ts(vr)-ts(vr(1)))+pD(1);

        p02 = Weighted_LPF(ts(vr),p0-p0m,fDm/2,ws); p0 = p02+p0m;
        pD2 = Weighted_LPF(ts(vr),pD-pDm,fDm/2,ws); pD = pD2+pDm;
end
f0 = diff(p0)/dt/2/pi; f0=[f0;f0(end)];
fD = diff(pD)/dt/2/pi; fD=[fD;fD(end)];


%% Correct the signal!
% output = Correct_Signal(sn(1:N),dt,p0(1:N),pD(1:N));
output = Correct_Signal(sn(vr),dt,p0,pD);
for ii=2:size(data.sn,2)
    o2 = Correct_Signal(data.sn(vr,ii),dt,p0,pD);
    output.Psn_fDc = [output.Psn_fDc, o2.Psn_fDc];
    output.Psn_f0c = [output.Psn_f0c, o2.Psn_f0c];
    output.Psn     = [output.Psn    , o2.Psn    ];
end
subplot(2,2,2); semilogy(output.fssDc/1e6,output.Psn_fDc);
xlim([-1,1]*1/dt/2/1e6);
xyt('Frequency (MHz)','PSD','Corrected PSD');

%% Expectation maximization
% Sometimes there are precision issues with the traces, so compute
% traces by collecting all terms and using compensated summation.
% Requires XSum on the path
if param.EM==1
    [Qt1,Qt2,Qt3,Qt4,Qt5,Qt6]=deal(zeros(4,4,Nb-1));
    for k=1:Nb-1
        Qt1(:,:,k) = Pkns{k+1}/(Md(k+1)/M0);
        Qt2(:,:,k) = xkns{k+1}*xkns{k+1}'/(Md(k+1)/M0);
        Qt3(:,:,k) = -Fs{k+1}*Exkxkp1p{k}/(Md(k+1)/M0);
        Qt4(:,:,k) = -Exkxkp1p{k}'*Fs{k+1}'/(Md(k+1)/M0);
        Qt5(:,:,k) = Fs{k+1}*Pkns{k}*Fs{k+1}'/(Md(k+1)/M0);
        Qt6(:,:,k) = Fs{k+1}*xkns{k}*xkns{k}'*Fs{k+1}'/(Md(k+1)/M0);
    end                   
%     Qtm = XSum(cat(3,Qt1,Qt2,Qt3,Qt4,Qt5,Qt6),3)/(Nf-1);
    wk = hanning(Nb-1); wk(1)=0;
    wk = reshape(wk,1,1,Nb-1);              % deweight the edges
    Qtm = XSum(cat(3,Qt1.*wk,Qt2.*wk,Qt3.*wk,Qt4.*wk,Qt5.*wk,Qt6.*wk),3)/sum(wk);
    output.Q_EM = Qtm;
%     output.Q_EM = trimmean( XSum(cat(4,Qt1,Qt2,Qt3,Qt4,Qt5,Qt6),4) , 20,3);
%     output.Q_EM = median( XSum(cat(4,Qt1,Qt2,Qt3,Qt4,Qt5,Qt6),4) , 3);

    trms=zeros(1+1+4+4,Nb);
    for k=1:Nb
        [hv,Hv]=hH(xkns{k},GI,GIn);
        Ck = [real(snr(:,k));imag(snr(:,k))]-hv+Hv*xkns{k};
        trms(:,k) = [Ck'*Ck;-Ck'*Hv*xkns{k};-diag(Hv'*Ck*xkns{k}');diag(Hv'*Hv*(Pkns{k}+xkns{k}*xkns{k}'))];
    end  
    wk = hanning(Nb); wk(1)=0;
    nstd2 = XSum(reshape(trms.*wk.',size(trms,1)*Nb,1))/sum(wk)/(2*M);   
    output.excess_noise_EM = sqrt(nstd2)/nstd*param.excess_noise; 
end

%% Set up the output structure, rescaling time axis
output.err = err;
output.npsd = 2*nstd^2*dt;   
output.QM = Q;
output.M = M;
output.p0 = p0;
output.pD = pD;
output.f0 = f0;
output.fD = fD;

output.tsDc = output.tsDc*Ts;
output.dtDc = output.dtDc*Ts;
output.fssDc = output.fssDc/Ts;
output.fss0c = output.fss0c/Ts;
output.fss   = output.fss  /Ts;
output.f0 = output.f0/Ts;
output.fD = output.fD/Ts;
output.Psn_fDc = output.Psn_fDc*Ts;
output.Psn_f0c = output.Psn_f0c*Ts;
output.Psn     = output.Psn    *Ts;
output.npsd    = output.npsd   *Ts;
output.data = data;
end

function Vf = Weighted_LPF(ts,Vs,freq,wk)
% Given cutoff frequency freq and weights wk, finds the x that minimizes
% sum(wk.*abs(xs-Vs).^2) under the constraint that the signal is 0 above
% freq. Basically ends up being similar to a Wiener filter.
%
    dt = mean(diff(ts)); N=length(ts);
    fss = fftshift(1/(N*dt)*[0:N-1]'); zi=find(fss==0); fss(1:zi-1)=fss(1:zi-1)-1/dt;
    Fss = fftshift(fft(Vs));
    Fss(abs(fss)>freq)=0;
    
    % Standard LPF. Doesn't work so well in pulsed mode
    Vf = real(ifft(ifftshift(Fss))); Vf=Vf(1:N);       
    
    % Weighted LPF. Basically does a weighted fit under the constraint that
    % the signal is band-limited.
    g = find(abs(fss)<=freq);       % always odd
    W = fftshift(fft(wk));
    P = fftshift(fft(Vs));
    
    function co=convf(xi)      % conv over filtered region
        if length(xi)==length(g), zix = find(fss(g)==0);       % short version
        else,                     zix = find(fss==0);    end   % long version
        co = ifft(fft(xi,length(xi)+length(W)-1).*fft(W,length(xi)+length(W)-1));   % fft conv
        ko = zi-1+zix;          % location in conv corresponding to f=0
        co = co(ko-(length(g)-1)/2:ko+(length(g)-1)/2);                             % filter region only
    end

    b = convf(P);    
    x0 = P(g); x=x0;
    r = b-convf(x);  p = r;
    for ii=1:5                  % CG iterations
        Ap = convf(p);
        ak = (r'*r)/(p'*Ap);
        x = x + ak*p;
        ro = r;
        r = r - ak*Ap;
        bk = r'*r/(ro'*ro);
        p = r + bk*p;
    end
%     dfigure; semilogy(fss(g),abs(b)  ,fss(g),abs  (convf(x0)),fss(g),abs  (convf(x)));
%     dfigure; plot    (fss(g),angle(b),fss(g),angle(convf(x0)),fss(g),angle(convf(x)));
    Vf = zeros(size(fss)); Vf(g) = x;
    Vf = real(ifft(ifftshift(Vf))); Vf=Vf(1:N);
end

function pout = Initialize_Params(pin)
    pout = pin;
    
    if ~isfield(pout,'plotme'),     pout.plotme = 1;            end
    if ~isfield(pout,'knownns'),    pout.knownns = NaN;         end
    if ~isfield(pout,'knownA'),     pout.knownA =  NaN;         end
    if ~isfield(pout,'knownfD'),    pout.knownfD = NaN;         end
    if ~isfield(pout,'knownf0'),    pout.knownf0 = NaN;         end
    
    if ~isfield(pout,'pulsed_mode_apodization'), pout.pulsed_mode_apodization = NaN;        end
    if ~isfield(pout,'initfrac'),                pout.initfrac = 0.1;                       end
    if ~isfield(pout,'Ninits'),                  pout.Ninits = 2;                           end
    if ~isfield(pout,'global_search_stds'),      pout.global_search_stds= 6;                end
    if ~isfield(pout,'global_search_maxsize'),   pout.global_search_maxsize= 1e6;           end
    if ~isfield(pout,'EM'),                      pout.EM= 1;                                end
    if ~isfield(pout,'post_regularization'),     pout.post_regularization= 'lpf';           end
end