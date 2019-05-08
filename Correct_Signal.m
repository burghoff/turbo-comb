function output = Correct_Signal(sn,dt,p0,pD)
N = length(sn);

% Correct offset fluctuations (multiplication)
sn_f0c = sn.*exp(-i*p0);

% Correct rep-rate fluctuations (non-uniform time stretch)
% mfD = mean(fD); mf0 = mean(f0);
mfD = (pD(end)-pD(1))/(2*pi*dt*(N-1)); mf0 = (p0(end)-p0(1))/(2*pi*dt*(N-1));

spp = ceil(1/(mfD*dt))*2;                                 % samples per period (samples per 2*pi)
NDc = floor( (max(pD)-min(pD))/(2*pi) )*spp;            % take an integer number of periods to ensure our Fourier frequencies match the corrected frequencies
pD2 = linspace( (max(pD)+min(pD))/2 - floor((max(pD)-min(pD))/(2*pi))*pi, ...
                (max(pD)+min(pD))/2 + floor((max(pD)-min(pD))/(2*pi))*pi,NDc+1)'; pD2 = pD2(1:end-1);

ts = dt*[0:N-1]';
ts2 = interp1(pD,ts,pD2);                                          % times to sample at

US_factor = 10;
Fs = fftshift(fft(sn_f0c));
zi=find(fftshift([0:N-1]')==0);
US_sn = ifft([Fs(zi:end);zeros(N*(US_factor-1),1);Fs(1:zi-1)])*US_factor;
sn_fDc = interp1(dt/US_factor*[0:N*US_factor-1]',US_sn,ts2);

dtDc = 1/(mfD*spp);
tsDc = dtDc*[0:NDc-1]';
fsDc = 1/dtDc/NDc*[0:NDc-1]';

% Compute the FFTs and plot it
fs = 1/(N*dt)*[0:N-1]';
fss   = fftshift(fs);   zi=find(fss==0);   fss(1:zi-1)  =fss(1:zi-1)-1/dt;
fssDc = fftshift(fsDc); zi=find(fssDc==0); fssDc(1:zi-1)=fssDc(1:zi-1)-1/dtDc;
pkdist = floor((max(pD)-min(pD))/(2*pi));
t1 = [zi:-pkdist:1]; t2 = [zi:pkdist:NDc];
pklocs = [fliplr(t1).';t2(2:end)'];

% Psn     = fftshift(abs(fft(sn)).^2)*dt/N;
% Psn_f0c = fftshift(abs(fft(sn_f0c)).^2)*dt/N;
% Psn_fDc = fftshift(abs(fft(sn_fDc)).^2)*dtDc/NDc;

Psn     = fftshift(abs(fft(sn    .*hanning(length(sn    )))).^2)*dt/N     * 1/mean(hanning(length(sn    )))^2;
Psn_f0c = fftshift(abs(fft(sn_f0c.*hanning(length(sn_f0c)))).^2)*dt/N     * 1/mean(hanning(length(sn_f0c)))^2;
Psn_fDc = fftshift(abs(fft(sn_fDc.*hanning(length(sn_fDc)))).^2)*dtDc/NDc * 1/mean(hanning(length(sn_fDc)))^2;

output.tsDc   = tsDc;
output.sn_fDc = sn_fDc;
output.NDc    = NDc;
output.dtDc   = dtDc;

output.fssDc  = fssDc+mf0;
output.Psn_fDc = Psn_fDc;
output.pklocs = pklocs;         % where we expect our comb lines to be

output.fss0c  = fss+mf0;
output.Psn_f0c = Psn_f0c;

output.fss  = fss;
output.Psn  = Psn;
end