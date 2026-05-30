function V=ifftshiftfreqs(F,varargin)
% A convenience function that undoes fftshiftfreqs
% Returns the original signal from its fft+fftshifted version
% 
% Inputs:
% F:  signal that has been fftshiftfreqs'ed
% varargin: any other arguments to be passed to ifft
%
% Outputs:
% V: original signal

if length(varargin)>1
    % IFFT has dim argument, shift along that direction
    dim = varargin{2};
    V = ifft(ifftshift(F,dim),varargin{:});
else
    V = ifft(ifftshift(F),varargin{:});
end

end