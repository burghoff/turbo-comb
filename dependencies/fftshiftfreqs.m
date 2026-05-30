function [fs,Fs]=fftshiftfreqs(V,dt,varargin)
% A convenience function that combines fft, fftshift, and fftfreqs
% Returns the fftshifted fft of a signal V, along with the corresponding
% shifted frequency axis
% 
% Inputs:
% V:  signal
% dt: corresponding time step
% varargin: any other arguments to be passed to fft
%
% Outputs:
% fs: frequency axis
% Fs: fft'ed and fftshift'ed signal

len = length(V);
if length(varargin)>0
% has Nfft argument
    len = varargin{1};
end
[~,fs]=fftfreqs(len,dt);

if length(varargin)>1
    % FFT has dim argument, shift along that direction
    dim = varargin{2};
    Fs = fftshift(fft(V,varargin{:}),dim);
else
    Fs = fftshift(fft(V,varargin{:}));
end

end