close all; clc; clearvars; %#ok<*UNRCH>

% Inputs
Directory = '../../data/';
TestPoint = '20210901003';
BlockLength = 2^10;
OverlapFactor = 0;
SensorSelection = 1:16;
fMax = 3e4;

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
WF = permute(WF,[1 3 2]);
clear aa wf;

% Sensor PSD
scanData = scanData(1:BlockIndex(2,end),SensorSelection)';
Y = complex(zeros(size(scanData,1),BlockLength,BlockNumber));
for aa=1:BlockNumber
    Y(:,:,aa) = fftshift(fft(Window.y.w.*scanData(:,BlockIndex(1,aa):BlockIndex(2,aa)),BlockLength,2),2);
end
Y = permute(Y,[1 3 2]);
% Y = reshape(permute(Y,[1 3 2]),length(SensorSelection)*BlockNumber,BlockLength);
% [filter.b,filter.a] = butter(SensorFilterLP(2),SensorFilterLP(1),'low','s');
% gain = reshape(freqs(filter.b,filter.a,reshape(Frequency.t,1,[])),1,[]);
% Y = Y.*gain;
clear aa scanData filter gain;


%%
close all;
wf = complex(zeros(size(WF)));
for aa=1:BlockLength
    if abs(Frequency.t(aa))>fMax
        wf(:,:,aa) = WF(:,:,aa);
    else
        Q = WF(:,:,aa);
        y = Y(:,:,aa);
        [PSI,LAMBDA] = eig(Q'*Q,'vector');
        PHI = Q*PSI;
        L = (PSI*y')/(y*y');
        PSIlse = L*y;
        PSIao = PSI-PSIlse;
        wf(:,:,aa) = (PSIao*PHI')';
    end
end

WF1 = reshape(squeeze(mean(abs(Window.wf.cw*WF).^2,2)/prod(BlockSize)/prod(RunLog.samplerate)),BlockSize);
wf1 = reshape(squeeze(mean(abs(Window.wf.cw*wf).^2,2)/prod(BlockSize)/prod(RunLog.samplerate)),BlockSize);

nSamples = size(WF1);
phase = pi*(2*rand(nSamples)-1);
% phase(nSamples(1)/2+1,nSamples(2)/2+1,nSamples(3)/2+1) = 0;
phase(2:nSamples(1),2:nSamples(2),nSamples(3)/2+2:nSamples(3)) = -flip(flip(flip(phase(2:nSamples(1),2:nSamples(2),2:nSamples(3)/2),1),2),3);



[sxx,frequency] = computeSXX(nanrms(reshape(real(ifftn(ifftshift(sqrt((WF1)*prod(RunLog.samplerate)*numel(WF1)).*exp(1i*phase)))).*WFInfo.Mask_WF,size(WF1,1)*size(WF1,2),size(WF1,3))),'positiveonly',1,'samplerate',RunLog.samplerate(3));
[sxx2] = computeSXX(nanrms(reshape(real(ifftn(ifftshift(sqrt((wf1)*prod(RunLog.samplerate)*numel(wf1)).*exp(1i*phase)))).*WFInfo.Mask_WF,size(wf1,1)*size(wf1,2),size(wf1,3))),'positiveonly',1,'samplerate',RunLog.samplerate(3));


%%
close(findobj('type','figure','number',1));
f1 = figure(1);
colororder(linspecer(5));

semilogx(frequency,10*log10(sxx./sxx2),'linewidth',1.25);
grid on;
xlabel('Frequency (Hz)','interpreter','latex');
ylabel('Reduction in OPD$_{RMS}(t)$ (dB/Hz)','interpreter','latex');
% xlim(10.^[2 5]);
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 5.5 3.25];

saveas(f1,'lse_mspod_freq.eps','epsc');
