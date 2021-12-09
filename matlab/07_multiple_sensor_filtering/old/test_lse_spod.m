close all; clc; clearvars;

% Inputs
Directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
TestPoint = '20210901003';
BlockLength = 2^10;
OverlapFactor = 0;
SensorSelection = 1:16;

% Load Data
[wf,CineInfo,RunLog,WFInfo] = loadWF([Directory TestPoint '_WF.mat'],'Scale',1e6,'ZernikeRemoval',1:3);
load([Directory TestPoint '_CDAQ.mat'],'scanData');

% Basic Setup
StepSize = BlockLength*2^-OverlapFactor;
BlockIndex = 1+(0:StepSize:size(wf,3)-BlockLength);
BlockIndex(2,:) = BlockIndex(1,:)+BlockLength-1;
BlockSize = [size(wf,[1 2]) BlockLength];
BlockNumber = size(BlockIndex,2);
DataLocation = reshape(WFInfo.Mask_WF==1,prod(BlockSize(1:2)),1);
LensletNumber = prod(BlockSize(1:2));
Frequency.yn = (-0.5:1/BlockSize(1):0.5-1/BlockSize(1))';
Frequency.xn = (-0.5:1/BlockSize(2):0.5-1/BlockSize(2));
Frequency.tn = permute((-0.5:1/BlockSize(3):0.5-1/BlockSize(3)),[1 3 2]);
Frequency.y = Frequency.yn*RunLog.samplerate(1);
Frequency.x = Frequency.xn*RunLog.samplerate(2);
Frequency.t = permute(Frequency.tn*RunLog.samplerate(3),[1 3 2]);

% Windowing Functions
Window.wf.s = createSpatialWindow(WFInfo.Mask_WF);
Window.wf.t = reshape(hann(BlockLength),1,1,BlockLength);
Window.wf.w = Window.wf.s.*Window.wf.t;
Window.wf.cw = 1/sqrt(sum(Window.wf.w.^2,'all')/numel(Window.wf.w));
Window.y.w = reshape(hann(BlockLength),1,BlockLength);
Window.y.cw = 1/sqrt(sum(Window.y.w.^2,'all')/numel(Window.y.w));

% Wavefront PSD In Time
wf = wf(:,:,1:BlockIndex(2,end));
wf(isnan(wf)) = 0;
WF = complex(zeros([BlockSize,BlockNumber]));
for aa=1:BlockNumber
    WF(:,:,:,aa) = fftshift(fftn(Window.wf.w.*wf(:,:,BlockIndex(1,aa):BlockIndex(2,aa))));
end
WF = reshape(WF,prod(BlockSize(1:2)),BlockSize(3),BlockNumber);
WF = permute(WF,[2 3 1]);
clear aa wf;



wf = squeeze(mean(conj(WF).*WF,2));





