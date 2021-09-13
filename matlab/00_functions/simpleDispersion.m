function [WF,frequency] = simpleDispersion(wf,varargin)
% Default Inputs
BlockSize = size(wf,3);
SampleRate = ones(1,3);
Output = 'log';
% Parse Varargin
for aa=1:2:length(varargin)-1
    switch lower(varargin{aa})
        case 'blocksize'
            BlockSize = varargin{aa+1};
        case 'samplerate'
            SampleRate = varargin{aa+1};
        case 'output'
            Output = varargin{aa+1};
        otherwise
            error('Invalid Input');
    end
end
% 
Mask = wf(:,:,1);
Mask(isnan(Mask)) = 0;
Mask(Mask~=0) = 1;
WindowS = createSpatialWindow(Mask);
WindowT = reshape(hann(BlockSize),1,1,[]);
Window = WindowS.*WindowT;
wf(isnan(wf)) = 0;
% 
wf = reshape(wf,size(wf,1),size(wf,2),BlockSize,size(wf,3)/BlockSize);
[WF,frequency] = simpleSXXn(wf(:,:,:,1),SampleRate,Window);
for aa=2:size(wf,4)
    WF = WF+simpleSXXn(wf(:,:,:,aa),SampleRate,Window);
end
WF = WF/size(wf,4);
if strcmpi(Output,'log')
    WF = log10(WF);
end
end

