function [WFdisp,Frequency] = WFdispersion(WF,varargin)
BlockSize = size(WF,3);
SampleRate = ones(1,ndims(WF));
SpatialWindow = ones(size(WF,1,2));
OverlapFactor = 1;
for aa=1:2:nargin-1
    if strcmpi(varargin{aa},'BlockSize')
        BlockSize = varargin{aa+1};
    elseif strcmpi(varargin{aa},'SampleRate')
        SampleRate = varargin{aa+1};
    elseif strcmpi(varargin{aa},'SpatialWindow')
        SpatialWindow = varargin{aa+1};
        SpatialWindow(isnan(SpatialWindow)) = 0;
    elseif strcmpi(varargin{aa},'OverlapFactor')
        OverlapFactor = varargin{aa+1};
    end
end
WF(isnan(WF)) = 0;
nBlocks = floor(size(WF,3)/BlockSize);
if floor((nBlocks+(OverlapFactor-1)/OverlapFactor)*BlockSize)>size(WF,3)
    nBlocks = nBlocks-1;
end
% Pad Size(WF,1,2) to a Power of 2
PadSize = 2.^nextpow2(size(WF,1,2))-size(WF,1,2);
WF = padarray(WF,floor(PadSize/2),0,'pre');
WF = padarray(WF,ceil(PadSize/2),0,'post');
% Windowing
SpatialWindow = padarray(SpatialWindow,floor(PadSize/2),0,'pre');
SpatialWindow = padarray(SpatialWindow,ceil(PadSize/2),0,'post');
TemporalWindow = reshape(hann(BlockSize),1,1,[]);
% Perform FFT
norm = numel(WF(:,:,:,1))*prod(SampleRate);
WFdisp = fftshift(mean(abs(fft(fft(fft(reshape(WF(:,:,1:nBlocks*BlockSize),size(WF,1),size(WF,2),BlockSize,[]).*(SpatialWindow.*TemporalWindow),[],1),[],2),[],3)).^2/norm,4))/OverlapFactor;
for aa=1:OverlapFactor-1
    WFdisp = WFdisp+fftshift(mean(abs(fft(fft(fft(reshape(WF(:,:,(1:BlockSize*(nBlocks-1))+floor(aa*BlockSize/OverlapFactor)),size(WF,1),size(WF,2),BlockSize,[]).*(SpatialWindow.*TemporalWindow),[],1),[],2),[],3)).^2/norm,4))/OverlapFactor;
end
% Calculate Frequency
Frequency.x = (-1/2:1/2.^nextpow2(size(WF,2)):1/2-1/2.^nextpow2(size(WF,2)))*SampleRate(2);
Frequency.y = (-1/2:1/2.^nextpow2(size(WF,1)):1/2-1/2.^nextpow2(size(WF,1)))*SampleRate(1);
Frequency.t = (-1/2:1/BlockSize:1/2-1/BlockSize)*SampleRate(3);
end
