%% Synthesize demonstration dual comb data
rng('default');                 % for reproducibility

dt = 1/(1.25e9);                % time step (s)
Ttotal=50e-6;                   % total time (s)
N = round(Ttotal/dt);           % total number of time steps

ns=[-50:50]';
peak_SNR = 10^6;
SNRs = exp(-ns.'.^2/2/(max(ns)/5).^2)*peak_SNR;
SNRs(1:5:end) = SNRs(1:5:end)/20;

f0m = 0; fDm = 5e6;        % mean offset and rep rate (Hz)
n0_max = fDm*5;    dn0 = fDm*1;
nD_max = fDm*.2;   dnD = 0.05*fDm;        
[nse0,nseD] = deal(zeros(N,1));
M=ceil(1/fDm/dt)+2;
for i1=1:N
    if i1==1
        nse0(i1)=0;-n0_max + 2*n0_max*rand(1);
        nseD(i1)=0;-nD_max + 2*nD_max*rand(1);
    else
        if mod(i1,M)==0
            nse0(i1)=median([-n0_max, nse0(i1-1) + randn(1)*dn0, n0_max]);
            nseD(i1)=median([-nD_max, nseD(i1-1) + randn(1)*dnD, nD_max]);
        else
            nse0(i1)=nse0(i1-1);
            nseD(i1)=nseD(i1-1);
        end
    end
end
f0 = f0m + nse0;
fD = fDm + nseD;

% Ensure noise is bandlimited
fs = 1/N/dt*[0:N-1]';
fss = fftshift(fs); zi=find(fss==0); fss(1:zi-1)=fss(1:zi-1)-1/dt;
Ff0 = fftshift(fft(f0)); Ff0(abs(fss)>fDm/2) = 0;
Ff0(abs(fss)<=fDm/2) = Ff0(abs(fss)<=fDm/2).*0.5.*(cos(fss(abs(fss)<=fDm/2)/(fDm/2)*pi)+1);
FfD = fftshift(fft(fD)); FfD(abs(fss)>fDm/2) = 0;
FfD(abs(fss)<=fDm/2) = FfD(abs(fss)<=fDm/2).*0.5.*(cos(fss(abs(fss)<=fDm/2)/(fDm/2)*pi)+1);
f0 = real(ifft(ifftshift(Ff0)));
fD = real(ifft(ifftshift(FfD)));


wn=1*(randn(N,1)+i*randn(N,1));                                     % white measurement noise
As= sqrt(2/N * 1^2 * SNRs' ) .* exp(1*i*rand(size(ns))*2*pi);       %As=As(randperm(length(ns))); 
nflr_fft = 1^2 * 2 * N; nflr_sig = 2/N*1^2;
% signal FFT = N^2*abs(As).^2, noise  FFT = wn^2 * 2 * N
% ergo signal equal to noise has As^2 = 2/N*wn^2

[Vn,Vc]=deal(wn);
p0 = cumsum(2*pi*f0*dt);
pD = cumsum(2*pi*fD*dt);
p0m = cumsum(2*pi*f0m*ones(N,1))*dt;
pDm = cumsum(2*pi*fDm*ones(N,1))*dt;
Vn = wn + exp(i*(p0 +ns.'.*pD ))*As;
Vc = wn + exp(i*(p0m+ns.'.*pDm))*As;

dfigure('DName','Generated data'); clf;
fs = 1/N/dt*[0:N-1]';
fss = fftshift(fs); zi=find(fss==0); fss(1:zi-1)=fss(1:zi-1)-1/dt;
semilogy(fss/1e6,fftshift(abs(fft(Vn)).^2),fss/1e6,fftshift(abs(fft(Vc)).^2));
xyt('Frequency (MHz)','Power (a.u)','');
legend('Clean','Corrupted');