close all; clc; clearvars; %#ok<*UNRCH>

% Only POD in a given frequency range

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
BlockSize = 2^10;
OverlapFactor = 0;
sensors = 1:16;
preFilter = 0;
postFilter = 0;
vFilter = [152 3 0.06 20];

% Load Data
[wf,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6,'ZernikeRemoval',1:3);
load([directory testPoint '_CDAQ.mat'],'scanData');

% Basic Setup
StepSize = BlockSize*2^-OverlapFactor;
BlockIndex = 1+(0:StepSize:size(wf,3)-BlockSize);
BlockIndex(2,:) = BlockIndex(1,:)+BlockSize-1;
M = size(wf,[1 2]);
N = size(BlockIndex,2);
freq.y = (-0.5:1/M(1):0.5-1/M(1))'*RunLog.samplerate(1);
freq.x = (-0.5:1/M(2):0.5-1/M(2))*RunLog.samplerate(2);
freq.t = permute((-0.5:1/BlockSize:0.5-1/BlockSize)*RunLog.samplerate(3),[1 3 2]);
freqNorm.y = (-0.5:1/M(1):0.5-1/M(1))';
freqNorm.x = (-0.5:1/M(2):0.5-1/M(2));
freqNorm.t = permute((-0.5:1/BlockSize:0.5-1/BlockSize),[1 3 2]);


% WF Sxx
windowS = createSpatialWindow(WFInfo.Mask_WF);
windowT = reshape(hann(BlockSize),1,1,BlockSize);
w = windowS.*windowT;
cw = 1/sqrt(sum(w.^2,'all')/numel(w));
Wsize = numel(w);
wf(isnan(wf)) = 0;
WF = zeros(M(1),M(2),BlockSize,N);
for aa=1:N
    WF(:,:,:,aa) = cw*fftshift(fftn(w.*wf(:,:,BlockIndex(1,aa):BlockIndex(2,aa))));
end
if preFilter
    dist = abs(vFilter(1)*RunLog.samplerate(2)/RunLog.samplerate(3)*freqNorm.x+vFilter(2)*RunLog.samplerate(1)/RunLog.samplerate(3)*freqNorm.y-freqNorm.t)/sqrt((vFilter(1)*RunLog.samplerate(2)/RunLog.samplerate(3))^2+(vFilter(2)*RunLog.samplerate(1)/RunLog.samplerate(3))^2+1);
    gain = sqrt(1./(1+(dist/vFilter(3)).^(+2*vFilter(4))));
    WF = WF.*gain;
end
WF = reshape(WF,prod(M),BlockSize,N);
WF = permute(WF,[1 3 2]);
clear aa cw dist gain w wf windowS windowT;

% Sensor Sxx
scanData = scanData';
w = reshape(hann(BlockSize),1,BlockSize);
cw = 1/sqrt(sum(w.^2,'all')/numel(w));
Y = zeros(size(scanData,1),BlockSize,N);
for aa=1:N
    Y(:,:,aa) = cw*fftshift(fft(w.*scanData(:,BlockIndex(1,aa):BlockIndex(2,aa)),BlockSize,2),2);
end
Y = permute(Y,[1 3 2]);
clear aa cw scanData w;

%% LSE-SPOD
WFao = zeros(prod(M),N,BlockSize);
for aa=1:BlockSize
    % Standard POD
    [phi,~] = eig(WF(:,:,aa)*WF(:,:,aa)');
    psi = (WF(:,:,aa)'*WF(:,:,aa))\WF(:,:,aa)'*phi;
    L = (psi'*Y(sensors,:,aa)')/(Y(sensors,:,aa)*Y(sensors,:,aa)');
    psi_lse = (L*Y(sensors,:,aa))';
    psi_ao = psi-psi_lse;
    WFao(:,:,aa) = (psi_ao*phi)';
    
    % Snapshot POD
%     [psi,~] = eig(WF(:,:,aa)'*WF(:,:,aa));
%     % psi = flip(psi,2);
%     phi = WF(:,:,aa)*psi;
%     L = (psi*Y(sensors,:,aa)')/(Y(sensors,:,aa)*Y(sensors,:,aa)');
%     psi_lse = L*Y(sensors,:,aa);
%     psi_ao = psi-psi_lse;
%     WFao(:,:,aa) = (psi_ao*phi')';
    clear L phi psi psi_ao psi_lse;
end
% if postFilter
%     dist = abs(vFilter(1)*RunLog.samplerate(2)/RunLog.samplerate(3)*freqNorm.x+vFilter(2)*RunLog.samplerate(1)/RunLog.samplerate(3)*freqNorm.y-freqNorm.t)/sqrt((vFilter(1)*RunLog.samplerate(2)/RunLog.samplerate(3))^2+(vFilter(2)*RunLog.samplerate(1)/RunLog.samplerate(3))^2+1);
%     gain = sqrt(1./(1+(dist/vFilter(3)).^(+2*vFilter(4))));
%     gain = permute(reshape(gain,prod(M),BlockSize),[1 3 2]);
%     WF = WF.*gain;
%     WFao = WFao.*gain;
% end
WF = (abs(WF)).^2/Wsize/prod(RunLog.samplerate);
WFao = (abs(WFao)).^2/Wsize/prod(RunLog.samplerate);
WF = reshape(squeeze(mean(WF,2)),M(1),M(2),BlockSize);
WFao = reshape(squeeze(mean(WFao,2)),M(1),M(2),BlockSize);
clear aa Y;


%% Plots
close all;
log_range = -14;
f1 = figure(1);
scolor = parula(2);
subplot(1,2,1);
patch(isocaps(freq.x,freq.y,freq.t(end/2+1:end),log10(WF(:,:,end/2+1:end)),log_range,'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq.x,freq.y,freq.t,log10(WF),log_range),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');
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
patch(isocaps(freq.x,freq.y,freq.t(end/2+1:end),log10(WFao(:,:,end/2+1:end)),log_range,'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq.x,freq.y,freq.t,log10(WFao),log_range),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');
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
sgtitle({testPoint;['Block Size = ' num2str(BlockSize)];['Overlap Factor = ' num2str(OverlapFactor)]},'Interpreter','Latex');





