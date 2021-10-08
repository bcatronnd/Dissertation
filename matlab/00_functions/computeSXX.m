function [sxx,frequency] = computeSXX(x,varargin)

% TODO
% - Input Error Check
% - Overlap Factors
% - Vector Dimension Input


% Defaults
average = 0;
[blockSize,dim] = max(size(x));
positiveFrequency = 0;
sampleRate = 1;
window = @hann;

% Parse Varagins
for aa=1:2:length(varargin)
    switch lower(varargin{aa})
        case 'average'
            average = varargin{aa+1};
        case 'blocksize'
            blockSize = varargin{aa+1};
        case 'dim'
            dim = varargin{aa+1};
        case {'positivefrequency','positivefreq','positiveonly'}
            positiveFrequency = varargin{aa+1};
        case 'samplerate'
            sampleRate = varargin{aa+1};
        case 'window'
            window = varargin{aa+1};
    end
end

% Input Error Check

% Resize x
blockNumber = floor(size(x,dim)/blockSize);
if mod(size(x,dim),blockSize)~=0
    dimorder = 1:ndims(x);
    dimorder(1) = dim;
    dimorder(dim) = 1;
    x = permute(x,dimorder);
    xSize = size(x);
    xSize(1) = blockSize*blockNumber;
    x = reshape(x,size(x,1),[]);
    x = x(1:blockSize*blockNumber,:);
    x = reshape(x,xSize);
    x = ipermute(x,dimorder);
    clear dimorder xSize;
end
if size(x,dim)~=blockSize
    dimorder = 1:ndims(x);
    dimorder(end) = dim;
    dimorder(dim) = ndims(x);
    x = permute(x,dimorder);
    dimsize = zeros(1,ndims(x)+1);
    for aa=1:ndims(x)-1
        dimsize(aa) = size(x,aa);
    end
    dimsize(end-1) = blockSize;
    dimsize(end) = blockNumber;
    x = reshape(x,dimsize);
    dimorder = [dimorder ndims(x)];
    x = ipermute(x,dimorder);
    clear dimorder dimsize;
end

% Window Check
if isa(window,'function_handle')
    window = window(blockSize);
    dimorder = 1:ndims(x);
    dimorder(1) = dim;
    dimorder(dim) = 1;
    window = permute(window,dimorder);
    clear dimorder;
end
if size(window,dim)~=blockSize
    error('Invalid Window Size');
end

% Main Code
cw = 1/sqrt(sum(window.^2,'all')/blockSize);
sxx = fftshift((abs(cw*fft(x.*window,[],dim))).^2,dim)/blockSize/sampleRate;
if blockNumber>1
    sxx = mean(sxx,ndims(sxx));
end
if average
    sxx = mean(sxx,average,'omitnan');
end

% Frequency
frequency = (-0.5:1/blockSize:0.5-1/blockSize)'*sampleRate;
dimorder = 1:ndims(x);
dimorder(1) = dim;
dimorder(dim) = 1;
frequency = permute(frequency,dimorder);

% Only Positives Out
if positiveFrequency
    % Frequency
    dimorder = 1:ndims(frequency);
    dimorder(1) = dim;
    dimorder(dim) = 1;
    frequency = permute(frequency,dimorder);
    xSize = size(frequency);
    xSize(1) = blockSize/2;
    frequency = reshape(frequency,size(frequency,1),[]);
    frequency = frequency(end/2+1:end,:);
    frequency = reshape(frequency,xSize);
    frequency = ipermute(frequency,dimorder);
    clear dimorder xSize;
    
    % Sxx
    dimorder = 1:ndims(sxx);
    dimorder(1) = dim;
    dimorder(dim) = 1;
    sxx = permute(sxx,dimorder);
    xSize = size(sxx);
    xSize(1) = blockSize/2;
    sxx = reshape(sxx,size(sxx,1),[]);
    sxx = sxx(end/2+1:end,:);
    sxx = reshape(sxx,xSize);
    sxx = ipermute(sxx,dimorder);
    clear dimorder xSize;
end
end

