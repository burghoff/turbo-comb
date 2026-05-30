function data_out = Picoscope_PSD(filename,blocks_to_load)
% filename: location of .mat file
% blocks_to_load: which blocks to load MATLAB convention

d = load(filename);
if nargin==1, blocks_to_load=[1:d.num_blocks]; end

N = d.samples_per_block;
num_blocks = d.num_blocks;
dt = d.timebase_s;
bin_filename = d.binary_filename; 

[pathstr,name,ext] = fileparts(filename);
if ~isempty(pathstr)
    bin_filename = [pathstr,'/',bin_filename];
end

% Convert range numbers to double vals
cs = d.settings.Channel_settings; crs=[10e-3,20e-3,50e-3,100e-3,200e-3,500e-3,1,2,5,10,20,50,100];
ranges = [cs.Channel_A_Settings.Range,cs.Channel_B_Settings.Range,cs.Channel_C_Settings.Range,cs.Channel_D_Settings.Range];
ranges = crs(ranges+1);

[fid,errmsg] = fopen(bin_filename);
if fid==-1 % file might've been renamed, retry with new name and bin extension
    if isempty(pathstr)
        [pathstr,name,ext] = fileparts(which(filename));
    end
    bin_filename = [pathstr,'/',name,'.bin'];
    [fid,errmsg] = fopen(bin_filename);
end
data_out = cell(length(blocks_to_load),1);
for ib=blocks_to_load
    fseek(fid, N*4*(ib-1)*2, 'bof');            % move to block, int16's are 2 bytes
    [A,cnt]=fread(fid,[N,4],'int16',0,'l');
    % Notes: stored little-endian
    % For each block, each channel stored sequentially (A1, B1, C1, D1, A2, B2,...)
    
    x=struct('dt',dt,'N',N,'settings',d.settings);
    A=A.*repmat(ranges,N,1)./(2^15);
    [x.A,x.B,x.C,x.D] = deal(A(:,1),A(:,2),A(:,3),A(:,4));
    data_out{ib==blocks_to_load}=x;
end
fclose(fid);
% psd = sum_psds/num_blocks;
% 
% 
% fs = 1/(N*dt)*[0:round(N/2)-1];
% psd=2*dt/N*psd(1:round(N/2));
% 
% ts = [0:N-1]'*dt;

end

