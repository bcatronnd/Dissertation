close all; clc; clearvars; %#ok<*UNRCH>

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
BlockSize = 2^10;
OverlapFactor = 1;
sensors = 1:2;

% Load Data
[wf,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);
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

% WF Sxx
windowS = createSpatialWindow(WFInfo.Mask_WF);
windowT = reshape(hann(BlockSize),1,1,BlockSize);
w = windowS.*windowT;
cw = 1/sqrt(sum(w.^2,'all')/numel(w));
wf(isnan(wf)) = 0;
WF = zeros(M(1),M(2),BlockSize,N);
for aa=1:N
    WF(:,:,:,aa) = (abs(cw*fftshift(fftn(w.*wf(:,:,BlockIndex(1,aa):BlockIndex(2,aa)))))).^2/numel(w)/prod(RunLog.samplerate);
end
WF = reshape(WF,prod(M),BlockSize,N);
WF = permute(WF,[1 3 2]);
clear aa cw w wf windowS windowT;

% Sensor Sxx
scanData = scanData';
w = reshape(hann(BlockSize),1,BlockSize);
cw = 1/sqrt(sum(w.^2,'all')/numel(w));
Y = zeros(size(scanData,1),BlockSize,N);
for aa=1:N
    Y(:,:,aa) = (abs(cw*fftshift(fft(w.*scanData(:,BlockIndex(1,aa):BlockIndex(2,aa)),BlockSize,2),2))).^2/BlockSize/RunLog.samplerate(3);
end
Y = permute(Y,[1 3 2]);
clear aa cw scanData w;

%% LSE-SPOD
WFao = zeros(prod(M),N,BlockSize);
for aa=1:BlockSize
    [psi,~] = eig(WF(:,:,aa)'*WF(:,:,aa));
    psi = flip(psi,2);
    phi = WF(:,:,aa)*psi;
    L = (psi*Y(sensors,:,aa)')/(Y(sensors,:,aa)*Y(sensors,:,aa)');
    psi_lse = L*Y(sensors,:,aa);
    psi_ao = psi-psi_lse;
    WFao(:,:,aa) = (psi_ao*phi')';
    clear L phi psi psi_ao psi_lse;
end
WF = reshape(squeeze(mean(WF,2)),M(1),M(2),BlockSize);
WFao = reshape(squeeze(mean(WFao,2)),M(1),M(2),BlockSize);
clear aa;

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
view(3);
material dull;
camlight;
f1.Children(1).TickLabelInterpreter = 'latex';

WFao(WFao<0) = 1e-21;
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
view(3);
material dull;
camlight;
f1.Children(1).TickLabelInterpreter = 'latex';



% disp(size(WF));
% disp(size(Y));









