close all; clc; clearvars; %#ok<*UNRCH>

% Inputs
Directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
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
wf = reshape(wf(:,:,1:BlockIndex(2,end)),prod(BlockSize(1:2)),BlockIndex(2,end));
wf(isnan(wf)) = 0;
WF = complex(zeros([LensletNumber,BlockLength,BlockNumber]));
for aa=1:BlockNumber
    WF(:,:,aa) = fftshift(fft(permute(Window.wf.t,[1 3 2]).*wf(:,BlockIndex(1,aa):BlockIndex(2,aa)),BlockLength,2),2);
end
WF = permute(WF,[1 3 2]);
% WF = WF(DataLocation,:,:);
clear aa wf;

% Sensor PSD
scanData = scanData(1:BlockIndex(2,end),SensorSelection)';
Y = complex(zeros(size(scanData,1),BlockLength,BlockNumber));
for aa=1:BlockNumber
    Y(:,:,aa) = fftshift(fft(Window.y.w.*scanData(:,BlockIndex(1,aa):BlockIndex(2,aa)),BlockLength,2),2);
end
Y = permute(Y,[1 3 2]);
clear aa scanData;


WFao = complex(zeros(length(DataLocation),BlockNumber,BlockLength));
for aa=1:BlockLength
    if abs(Frequency.t(aa))>fMax
        WFao(DataLocation,:,aa) = WF(DataLocation,:,aa);
    else
        Q = WF(DataLocation,:,aa);
        y = Y(:,:,aa);
        [PSI,LAMBDA] = eig(Q'*Q,'vector');
        PHI = Q*PSI;
        L = (PSI*y')/(y*y');
        PSIlse = L*y;
        PSIao = PSI-PSIlse;
        WFao(DataLocation,:,aa) = (PSIao*PHI')';
    end
end
clear aa PSI PHI PSIao PSIlse Q y;

% Multi-Dimensional Power Spectra - Spatial
WF = reshape(ipermute(WF,[1 3 2]),[BlockSize BlockNumber]);
WF = fftshift(fftshift(fft(fft(Window.wf.s.*WF,BlockSize(1),1),BlockSize(2),2),1),2);
WF = mean(abs(Window.wf.cw*WF).^2/prod(BlockSize)/prod(RunLog.samplerate),4);
WFao = reshape(ipermute(WFao,[1 3 2]),[BlockSize BlockNumber]);
WFao = fftshift(fftshift(fft(fft(Window.wf.s.*WFao,BlockSize(1),1),BlockSize(2),2),1),2);
WFao = mean(abs(Window.wf.cw*WFao).^2/prod(BlockSize)/prod(RunLog.samplerate),4);

disp(sqrt(mean(WF,'all')*prod(RunLog.samplerate)));
disp(sqrt(mean(WFao,'all')*prod(RunLog.samplerate)));

%%
close all;
log_range = -14;
scolor = parula(2);

f1 = figure(1);
subplot(1,2,1);
patch(isocaps(Frequency.x,Frequency.y(end/2+1:end),Frequency.t(end/2+1:end),log10(WF(end/2+1:end,:,end/2+1:end)),log_range,'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(Frequency.x,Frequency.y(end/2+1:end),Frequency.t(end/2+1:end),log10(WF(end/2+1:end,:,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
grid on;
daspect([1 1 50]);
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[-1 1]);
zlim(RunLog.samplerate(3)/2*[0 1]);
xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
zlabel('$f\ (Hz)$','Interpreter','Latex');
title('Unfiltered','Interpreter','Latex');
view([0 0]);
material dull;
camlight;
f1.Children(1).TickLabelInterpreter = 'latex';

subplot(1,2,2);
patch(isocaps(Frequency.x,Frequency.y(end/2+1:end),Frequency.t(end/2+1:end),log10(WFao(end/2+1:end,:,end/2+1:end)),log_range,'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(Frequency.x,Frequency.y(end/2+1:end),Frequency.t(end/2+1:end),log10(WFao(end/2+1:end,:,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
grid on;
daspect([1 1 50]);
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[-1 1]);
zlim(RunLog.samplerate(3)/2*[0 1]);
xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
zlabel('$f\ (Hz)$','Interpreter','Latex');
title('Filtered','Interpreter','Latex');
view([0 0]);
material dull;
camlight;
f1.Children(1).TickLabelInterpreter = 'latex';


f1.Units = 'inches';
f1.Position = [1 1 5.5 3.5];

saveas(f1,'lse_spod.eps','epsc');










