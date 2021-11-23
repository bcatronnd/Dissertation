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
DataLocation = WFInfo.Mask_WF==1;
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
WF = mean((abs(Window.wf.cw*WF)).^2/prod(BlockSize)/prod(RunLog.samplerate),4);
clear aa wf;

% Sensor PSD
scanData = scanData(1:BlockIndex(2,end),SensorSelection)';
Y = complex(zeros(size(scanData,1),BlockLength,BlockNumber));
for aa=1:BlockNumber
    Y(:,:,aa) = fftshift(fft(Window.y.w.*scanData(:,BlockIndex(1,aa):BlockIndex(2,aa)),BlockLength,2),2);
end
Y = mean((abs(Window.y.cw*Y)).^2/BlockLength/RunLog.samplerate(3),3);
clear aa scanData;


%% Frequency Bin Plot
close all;
% fBin = [randi(BlockSize(1)) randi(BlockSize(2))];
fBin = [10 27];

log_range = -14;
f1 = figure(1);
subplot(3,2,[1 3 5]);
scolor = parula(2);
patch(isocaps(Frequency.x,Frequency.y,Frequency.t(end/2+1:end),log10(WF(:,:,end/2+1:end)),log_range,'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(Frequency.x,Frequency.y,Frequency.t,log10(WF),log_range),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');
hold on;
plot3(Frequency.x(fBin(2))*[1 1],Frequency.y(fBin(1))*[1 1],RunLog.samplerate(3)/2*[0 1],'m-','linewidth',2);
hold off;
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


error = zeros(1,100);
coeff = zeros(16,100);

% aa=1;

penalty = 1;
X0 = permute(WF(fBin(1),fBin(2),end/2+1:end),[1 3 2]);
Y0 = Y(:,end/2+1:end);
F0 = squeeze(Frequency.t(end/2+1:end));
W0 = 10.^(-log10(1:length(F0))-1/3);

% for aa=1:1
coeff(3,1) = (X0.*W0)/(Y0(3,:).*W0);
% end





% coeff(:,1) = ((X0*Y0')/(Y0*Y0'))';
% coeff(:,1) = ((log(X0)*log(Y0)')/(log(Y0)*log(Y0)'))';
% coeff(:,1) = rand(16,1)-2.5;
R = X0-sum(Y0.*exp(coeff(:,1)),1);
E = R;
E(R<0) = penalty*R(R<0);
error(1) = sum(E.^2);
% coeff(:,2) = coeff(:,1).*((rand(16,1)-0.5)/100+1);
% for aa=2:2
%     R = X0-sum(Y0.*exp(coeff(:,aa)),1);
%     E = R;
%     E(R<0) = penalty*R(R<0);
%     error(aa) = sum(E.^2);
%     coeff(:,aa+1) = coeff(:,aa)-1*error(aa)/(error(aa)-error(aa-1))*(coeff(:,aa)-coeff(:,aa-1));
%     
% end


% [~,index] = min(error);
index = 1;
R = X0-sum(Y0.*exp(coeff(:,index)),1);
    E = R;
    E(R<0) = penalty*R(R<0);

subplot(3,2,2);
loglog(F0,X0,'k-',F0,sum(coeff(:,index).*Y0,1));
grid on;

subplot(3,2,4);
loglog(F0,abs(R),F0,abs(E));
grid on;

subplot(3,2,6);
semilogy(error);
grid on;

f1.Units = 'inches';
f1.Position = [1 1 8 4];

