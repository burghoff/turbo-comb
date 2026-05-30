function os = AK_Params(varargin)

% Device params
os.f0_guess  = 0;                   % initial guess of offset (Hz)
os.RTS       = 1;                     % use RTS smoothing? 
os.plotme    = 1;                     % plot results?

os.fD_range = NaN;
os.Q = diag([1,.01,.1,0.001]).^2;
os.excess_noise = 1;

os.init_batches = 50;
os.Ninits = Inf;

os.inittol = .01;
os.maxinits = Inf;

os.global_search_stds = 6;
os.global_search_maxsize = 1e6;

os.post_regularization = 'lpf';
os.EM = 1;

fields = fieldnames(os);
for ii=1:length(varargin)
    for jj=1:length(fields)
        if ischar(varargin{ii}) & isequal(varargin{ii},fields{jj})
            os.(fields{jj}) = varargin{ii+1};
        end
    end     
end

end