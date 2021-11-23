close all; clc; clearvars; %#ok<*UNRCH>

% Inputs
Directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
TestPoint = '20210901003';
BlockLength = 2^10;
OverlapFactor = 0;
SensorSelection = 1:16;
SensorFilterLP = [1e4 5];

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
WF = permute(WF,[2 3 1]);
clear aa wf;

% Sensor PSD
scanData = scanData(1:BlockIndex(2,end),SensorSelection)';
Y = complex(zeros(size(scanData,1),BlockLength,BlockNumber));
for aa=1:BlockNumber
    Y(:,:,aa) = fftshift(fft(Window.y.w.*scanData(:,BlockIndex(1,aa):BlockIndex(2,aa)),BlockLength,2),2);
end
Y = reshape(permute(Y,[1 3 2]),length(SensorSelection)*BlockNumber,BlockLength);
[filter.b,filter.a] = butter(SensorFilterLP(2),SensorFilterLP(1),'low','s');
gain = reshape(freqs(filter.b,filter.a,reshape(Frequency.t,1,[])),1,[]);
Y = Y.*gain;
clear aa scanData filter gain;

% LSE-SPOD
WFao = complex(zeros(BlockLength,BlockNumber,LensletNumber));
wb = waitbar(0);
for aa=1:LensletNumber
    if DataLocation(aa)
        [PHI,~] = eig(WF(:,:,aa)*WF(:,:,aa)'/BlockNumber,'vector');
        PSI = (WF(:,:,aa)'*PHI);
        PSIphase = angle(PSI);
        L = (PSI*Y')/(Y*Y');
        PSIlse = L*Y;
        PSIao = abs(PSI-PSIlse).*exp(1i*PSIphase);
        WFao(:,:,aa) = (PSIao*PHI')';
    end
    waitbar(aa/LensletNumber,wb);
end
close(wb);

% Multi-Dimensional Power Spectra - Spatial
WF = reshape(ipermute(WF,[2 3 1]),[BlockSize BlockNumber]);
WF = fftshift(fftshift(fft(fft(Window.wf.s.*WF,BlockSize(1),1),BlockSize(2),2),1),2);
WF = mean(abs(Window.wf.cw*WF).^2/prod(BlockSize)/prod(RunLog.samplerate),4);
WFao = reshape(ipermute(WFao,[2 3 1]),[BlockSize BlockNumber]);
WFao = fftshift(fftshift(fft(fft(Window.wf.s.*WFao,BlockSize(1),1),BlockSize(2),2),1),2);
WFao = mean(abs(Window.wf.cw*WFao).^2/prod(BlockSize)/prod(RunLog.samplerate),4);

disp(['OPDrms reduction: ' num2str((1-sqrt(mean(WFao*prod(RunLog.samplerate),'all'))/sqrt(mean(WF*prod(RunLog.samplerate),'all')))*100,'%0.1f') '%']);
save(['lse_spod2_20210901003' replace(num2str(SensorSelection,'_%.0u'),' ','') '.mat'],'WF','WFao','Frequency');

%% Plots
close all;
log_range = -13;
f1 = figure(1);
scolor = parula(2);
subplot(1,2,1);
patch(isocaps(Frequency.x,Frequency.y,Frequency.t(end/2+1:end),log10(WF(:,:,end/2+1:end)),log_range,'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(Frequency.x,Frequency.y,Frequency.t,log10(WF),log_range),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');
grid on;
daspect([1 1 50]);
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[-1 1]);
zlim(RunLog.samplerate(3)/2*[0 1]);
xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
zlabel('$f\ (Hz)$','Interpreter','Latex');
title('Unfiltered','Interpreter','Latex');
view(3);
material dull;
camlight;
view(0,0);
f1.Children(1).TickLabelInterpreter = 'latex';
subplot(1,2,2);
patch(isocaps(Frequency.x,Frequency.y,Frequency.t(end/2+1:end),log10(WFao(:,:,end/2+1:end)),log_range,'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(Frequency.x,Frequency.y,Frequency.t,log10(WFao),log_range),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');
grid on;
daspect([1 1 50]);
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[-1 1]);
zlim(RunLog.samplerate(3)/2*[0 1]);
xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
zlabel('$f\ (Hz)$','Interpreter','Latex');
title('Filtered','Interpreter','Latex');
view(3);
material dull;
camlight;
view(0,0);
f1.Children(1).TickLabelInterpreter = 'latex';
sgtitle({TestPoint;['Block Length = ' num2str(BlockLength)];['Overlap Factor = ' num2str(OverlapFactor)]},'Interpreter','Latex');






% %%
% close all;
% figure(1);
% loglog(squeeze(Frequency.t),mean(abs(squeeze(WF(:,:,aa))).^2,2));
% hold on;
% loglog(squeeze(Frequency.t),mean(abs((PSI*PHI')').^2,2));
% loglog(squeeze(Frequency.t),mean(abs((PSIlse*PHI')').^2,2));
% loglog(squeeze(Frequency.t),mean(abs((PSIao*PHI')').^2,2));
% grid on;
% legend('Orginial','$\Psi\cdot\Phi$','$\Psi_{LSE}\cdot\Phi$','$\Psi_{AO}\cdot\Phi$','interpreter','latex');



