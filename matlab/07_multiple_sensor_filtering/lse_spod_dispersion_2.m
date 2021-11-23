close all; clc; clearvars; %#ok<*UNRCH>

% Only POD in a given frequency range

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
TestPoint = '20210901003';
BlockLength = 2^10;
OverlapFactor = 0;
PODdimension = 2;
SensorSelection = [1:2 7:16];
% preFilter = 0;
% postFilter = 0;
% vFilter = [152 3 0.06 20];

% Load Data
[wf,CineInfo,RunLog,WFInfo] = loadWF([directory TestPoint '_WF.mat'],'Scale',1e6,'ZernikeRemoval',1:3);
load([directory TestPoint '_CDAQ.mat'],'scanData');

% Basic Setup
StepSize = BlockLength*2^-OverlapFactor;
BlockIndex = 1+(0:StepSize:size(wf,3)-BlockLength);
BlockIndex(2,:) = BlockIndex(1,:)+BlockLength-1;
BlockSize = [size(wf,[1 2]) BlockLength];
BlockNumber = size(BlockIndex,2);
freq.y = (-0.5:1/BlockSize(1):0.5-1/BlockSize(1))'*RunLog.samplerate(1);
freq.x = (-0.5:1/BlockSize(2):0.5-1/BlockSize(2))*RunLog.samplerate(2);
freq.t = permute((-0.5:1/BlockSize(3):0.5-1/BlockSize(3))*RunLog.samplerate(3),[1 3 2]);
freqNorm.y = (-0.5:1/BlockSize(1):0.5-1/BlockSize(1))';
freqNorm.x = (-0.5:1/BlockSize(2):0.5-1/BlockSize(2));
freqNorm.t = permute((-0.5:1/BlockSize(3):0.5-1/BlockSize(3)),[1 3 2]);

% WF Sxx
windowS = createSpatialWindow(WFInfo.Mask_WF);
windowT = reshape(hann(BlockSize(3)),1,1,BlockSize(3));
w = windowS.*windowT;
cw = 1/sqrt(sum(w.^2,'all')/numel(w));
wf(isnan(wf)) = 0;
WF = zeros([BlockSize,BlockNumber]);
for aa=1:BlockNumber
    WF(:,:,:,aa) = cw*fftshift(fftn(w.*wf(:,:,BlockIndex(1,aa):BlockIndex(2,aa))));
end
% if preFilter
%     dist = abs(vFilter(1)*RunLog.samplerate(2)/RunLog.samplerate(3)*freqNorm.x+vFilter(2)*RunLog.samplerate(1)/RunLog.samplerate(3)*freqNorm.y-freqNorm.t)/sqrt((vFilter(1)*RunLog.samplerate(2)/RunLog.samplerate(3))^2+(vFilter(2)*RunLog.samplerate(1)/RunLog.samplerate(3))^2+1);
%     gain = sqrt(1./(1+(dist/vFilter(3)).^(+2*vFilter(4))));
%     WF = WF.*gain;
% end
WF = reshape(WF,prod(BlockSize(1:2)),BlockSize(3),BlockNumber);
WF = permute(WF,[2 3 1]);
clear aa cw dist gain w wf windowS windowT;

% Sensor Sxx
scanData = scanData';
w = reshape(hann(BlockSize(3)),1,BlockSize(3));
cw = 1/sqrt(sum(w.^2,'all')/numel(w));
Y = zeros(size(scanData,1),BlockSize(3),BlockNumber);
for aa=1:BlockNumber
    Y(:,:,aa) = cw*fftshift(fft(w.*scanData(:,BlockIndex(1,aa):BlockIndex(2,aa)),BlockSize(3),2),2);
end
clear aa;
Y = permute(Y,[3 2 1]);
% Y4 = conj(Y)./Y./conj(Y).*Y;
% Y4 = complex(zeros(BlockNumber,BlockNumber,length(SensorSelection)));
% for aa=1:length(SensorSelection)
%     Y4(:,:,aa) = Y(:,:,aa)'/(Y(:,:,aa)*Y(:,:,aa)')*Y(:,:,aa);
% end
clear aa cw scanData w;

%% LSE-SPOD
WFao = complex(zeros(BlockSize(3),BlockNumber,prod(BlockSize(1:2))));
for aa=1%:prod(BlockSize(1:2))
    % Snapshot Method
%     [psi,~] = eig(WF(:,:,aa)'*WF(:,:,aa));
%     % psi = flip(psi,2);
%     phi = WF(:,:,aa)*psi;
%     psi_lse = complex(zeros(BlockNumber,BlockNumber,length(SensorSelection)));
%     for bb=1:length(SensorSelection)
%         psi_lse(:,:,bb) = Y4(:,:,bb)*psi;
%     end
%     psi_ao = psi-sum(psi_lse,3);
%     WFao(:,:,aa) = phi*psi_ao;
%     clear bb L phi psi psi_ao psi_lse;

    % Standard Method
    [phi,~] = eig(WF(:,:,aa)*WF(:,:,aa)');
    phi = flip(phi,2);
    psi = (WF(:,:,aa)'*phi);
    L = sum(psi.*conj(Y(:,:,SensorSelection)),2)./sum(Y(:,:,SensorSelection).*conj(Y(:,:,SensorSelection)),2);
    psi_lse = sum(L.*Y(:,:,SensorSelection),3);
    psi_ao = psi-psi_lse;
    WFao(:,:,aa) = (psi_ao*phi)';
    
    
    loglog(squeeze(freq.t),(abs(psi(1,:))).^2,squeeze(freq.t),(abs(psi_lse(1,:))).^2,squeeze(freq.t),(abs(psi_ao(1,:))).^2)
end
% if postFilter
%     dist = abs(vFilter(1)*RunLog.samplerate(2)/RunLog.samplerate(3)*freqNorm.x+vFilter(2)*RunLog.samplerate(1)/RunLog.samplerate(3)*freqNorm.y-freqNorm.t)/sqrt((vFilter(1)*RunLog.samplerate(2)/RunLog.samplerate(3))^2+(vFilter(2)*RunLog.samplerate(1)/RunLog.samplerate(3))^2+1);
%     gain = sqrt(1./(1+(dist/vFilter(3)).^(+2*vFilter(4))));
%     gain = permute(reshape(gain,prod(M),BlockSize),[1 3 2]);
%     WF = WF.*gain;
%     WFao = WFao.*gain;
% end
% WF = reshape(ipermute(WF,[2 3 1]),[BlockSize BlockNumber]);
% WFao = reshape(ipermute(WFao,[2 3 1]),[BlockSize BlockNumber]);
% WF = mean((abs(WF)).^2/prod(BlockSize)/prod(RunLog.samplerate),4);
% WFao = mean((abs(WFao)).^2/prod(BlockSize)/prod(RunLog.samplerate),4);
% clear aa Y;

%% Plots
% close all;
% log_range = -14;
% f1 = figure(1);
% scolor = parula(2);
% subplot(1,2,1);
% patch(isocaps(freq.x,freq.y,freq.t(end/2+1:end),log10(WF(:,:,end/2+1:end)),log_range,'all'),'facecolor','interp','edgecolor','none','facelighting','none');
% patch(isosurface(freq.x,freq.y,freq.t,log10(WF),log_range),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');
% grid on;
% daspect([1 1 50]);
% xlim(RunLog.samplerate(1)/2*[-1 1]);
% ylim(RunLog.samplerate(2)/2*[-1 1]);
% zlim(RunLog.samplerate(3)/2*[0 1]);
% xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
% ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
% zlabel('$f\ (Hz)$','Interpreter','Latex');
% title('Unfiltered','Interpreter','Latex');
% view(3);
% material dull;
% camlight;
% view(0,0);
% f1.Children(1).TickLabelInterpreter = 'latex';
% subplot(1,2,2);
% patch(isocaps(freq.x,freq.y,freq.t(end/2+1:end),log10(WFao(:,:,end/2+1:end)),log_range,'all'),'facecolor','interp','edgecolor','none','facelighting','none');
% patch(isosurface(freq.x,freq.y,freq.t,log10(WFao),log_range),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');
% grid on;
% daspect([1 1 50]);
% xlim(RunLog.samplerate(1)/2*[-1 1]);
% ylim(RunLog.samplerate(2)/2*[-1 1]);
% zlim(RunLog.samplerate(3)/2*[0 1]);
% xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
% ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
% zlabel('$f\ (Hz)$','Interpreter','Latex');
% title('Filtered','Interpreter','Latex');
% view(3);
% material dull;
% camlight;
% view(0,0);
% f1.Children(1).TickLabelInterpreter = 'latex';
% sgtitle({TestPoint;['Block Length = ' num2str(BlockLength)];['Overlap Factor = ' num2str(OverlapFactor)]},'Interpreter','Latex');

% saveas(f1,['lse_spod_disp_' testPoint '_' num2str(BlockSize) '_' num2str(OverlapFactor) '.png'],'png');





%%
% clim = [-17 -7];
% fPoints = [0 2.5e4];
% f2 = figure(2);
% subplot(3,1,1);
% plot(fPoints,fPoints/RunLog.u,'k-');
% hold on;
% surf(squeeze(freq.t),freq.x,squeeze(log10(WF(end/2+1,:,:))),'linestyle','none','facecolor','interp');
% view(2);
% grid on;
% xlim(RunLog.samplerate(3)/2*[0 1]);
% ylim(125*[-1 1]);
% caxis(clim);
% % colorbar;
% subplot(3,1,2);
% plot(fPoints,fPoints/RunLog.u,'k-');
% hold on;
% surf(squeeze(freq.t),freq.x,squeeze(log10(WFao(end/2+1,:,:))),'linestyle','none','facecolor','interp');
% view(2);
% grid on;
% xlim(RunLog.samplerate(3)/2*[0 1]);
% ylim(125*[-1 1]);
% caxis(clim);
% colorbar;
% subplot(3,1,3);
% plot(fPoints,fPoints/RunLog.u,'k-');
% hold on;
% surf(squeeze(freq.t),freq.x,10*log10(squeeze(WF(end/2+1,:,:))./squeeze(WFao(end/2+1,:,:))),'linestyle','none','facecolor','interp');
% view(2);
% grid on;
% xlim(RunLog.samplerate(3)/2*[0 1]);
% ylim(125*[-1 1]);
% colorbar;
% 
% f2.Children(1).TickLabelInterpreter = 'latex';
% f2.Children(2).TickLabelInterpreter = 'latex';
% f2.Children(3).TickLabelInterpreter = 'latex';
% f2.Children(4).TickLabelInterpreter = 'latex';
% f2.Children(5).TickLabelInterpreter = 'latex';
% f2.Children(2).Layer = 'top';
% f2.Children(4).Layer = 'top';
% f2.Children(5).Layer = 'top';
% 
% f2.Units = 'inches';
% f2.Position = [1 1 5.5 6];
% 
% f2.Children(5).Position(3) = f2.Children(4).Position(3);
% f2.Children(4).Position(1) = f2.Children(4).Position(1);
% f2.Children(3).Position(1) = f2.Children(3).Position(1);
% f2.Children(2).Position(1) = f2.Children(2).Position(1);
% f2.Children(3).Position(4) = f2.Children(5).Position(4)+f2.Children(5).Position(2)-f2.Children(4).Position(2);
% f2.Children(1).Position(1) = f2.Children(3).Position(1);
% f2.Children(1).Position(3) = f2.Children(3).Position(3);
% f2.Children(2).Position(3) = f2.Children(4).Position(3);
% caxis(clim);





% disp(size(WF));
% disp(size(Y));



%%
% figure(3)
% nPlots = 8;
% psiInv = inv(psi);
% for aa=1:nPlots^2
%     subplot(nPlots,nPlots,aa)
% %     surf((abs(reshape(sum(phi*psiInv(:,1:aa),2),M(1),M(2)))).^2,'linestyle','none');
%     surf(reshape((abs(phi(:,aa))).^2,M(1),M(2)),'linestyle','none');
%     view(2);
%     axis image off;
% end




