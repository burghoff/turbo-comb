function [fs,fss]=fftfreqs(N,dt)
% Generates the frequencies and shifted versions of FFT quantities that are
% N long and have time spacing of dt. Generates column vectors.

fs = 1/(N*dt)*[0:N-1]';
fss = fftshift(fs); zi=find(fss==0); fss(1:zi-1)=fss(1:zi-1)-1/dt;

end