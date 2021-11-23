close all; clc; clearvars; %#ok<*UNRCH>

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
BlockSize = 2^10;
PartialBlockSize = 2^8;
OverlapFactor = 0;
sensors = [1:2 7:16];
sensorPower = 1;

% Load Data
[wf,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);
load([directory testPoint '_CDAQ.mat'],'scanData');

% Basic Setup
StepSize = BlockSize*2^-OverlapFactor;
BlockIndex = 1+(0:StepSize:size(wf,3)-BlockSize);
BlockIndex(2,:) = BlockIndex(1,:)+BlockSize-1;
PartialBlockIndex = 1+(0:PartialBlockSize:BlockSize-PartialBlockSize);
PartialBlockIndex(2,:) = PartialBlockIndex(1,:)+PartialBlockSize-1;
M = size(wf,[1 2]);
N = size(BlockIndex,2);
freq.y = (-0.5:1/M(1):0.5-1/M(1))'*RunLog.samplerate(1);
freq.x = (-0.5:1/M(2):0.5-1/M(2))*RunLog.samplerate(2);
freq.t = permute((-0.5:1/BlockSize:0.5-1/BlockSize)*RunLog.samplerate(3),[1 3 2]);
freqNorm.y = (-0.5:1/M(1):0.5-1/M(1))';
freqNorm.x = (-0.5:1/M(2):0.5-1/M(2));
freqNorm.t = permute((-0.5:1/BlockSize:0.5-1/BlockSize),[1 3 2]);
dataLocation = reshape(WFInfo.Mask_WF==1,prod(M),1);
window.s = createSpatialWindow(WFInfo.Mask_WF);
window.t = permute(hann(BlockSize),[3 2 1]);
% window.c = window.s.*window.t;
window.cw = 1/sqrt(sum((window.s.*window.t).^2,'all')/(prod(M)*BlockSize));

% Reshape And Trim Data Arrays
wf = reshape(wf,prod(M),size(wf,3));
wf = wf(dataLocation,:);
scanData = permute(scanData,[2 1]);
window.t = permute(window.t,[1 3 2]);
freqIndex = ifftshift(1:BlockSize);

WF = zeros(M(1),M(2),BlockSize);

% Step Through Temporal Frequency
tic;
wb = waitbar(0);
for aa=1:BlockSize/PartialBlockSize
    % Partial Fourier Transform Of Datasets
    WFs = complex(zeros(size(wf,1),PartialBlockSize,N));
    Ys = complex(zeros(size(scanData,1),PartialBlockSize,N));
    for bb=1:N
        WFs(:,:,bb) = pft(wf(:,BlockIndex(1,bb):BlockIndex(2,bb)).*window.t,PartialBlockIndex(1,aa):PartialBlockIndex(2,aa),BlockSize,2);
        Ys(:,:,bb) = pft(scanData(:,BlockIndex(1,bb):BlockIndex(2,bb)).*window.t,PartialBlockIndex(1,aa):PartialBlockIndex(2,aa),BlockSize,2);
    end
    clear bb;
    WFs = permute(WFs,[1 3 2]);
    Ys = permute(Ys,[1 3 2]).^sensorPower;
    for bb=1:PartialBlockSize
        % LSE-SPOD
        [psi,~] = eig(WFs(:,:,bb)'*WFs(:,:,bb));
        phi = WFs(:,:,bb)*psi;
        L = (psi*Ys(sensors,:,bb)')/(Ys(sensors,:,bb)*Ys(sensors,:,bb)');
        WFao = complex(zeros(prod(M),N));
        WFao(dataLocation,:) = ((psi-L*Ys(sensors,:,bb))*phi')';
        WFao = reshape(WFao,M(1),M(2),N);
        % ND-SXX
        WF(:,:,PartialBlockIndex(1,aa)+bb-1) = mean((abs(window.cw*fft(fft(window.s.*WFao,[],1),[],2))).^2/(prod(M)*BlockSize)/prod(RunLog.samplerate),3);
    end
    clear bb;
    waitbar(aa/(BlockSize/PartialBlockSize),wb);
end
close(wb);
WF = fftshift(WF);
toc;

save(['lse_spod_' testPoint '_' num2str(BlockSize) '_' num2str(OverlapFactor) '.mat'],'WF','freq','RunLog');

%% Plots
close all;
log_range = -14;
f1 = figure(1);
scolor = parula(2);
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
% title('Unfiltered','Interpreter','Latex');
view(3);
material dull;
camlight;












% WF-FFT
% wf(isnan(wf)) = 0;
% StepSize = BlockSize*2^-OverlapFactor;
% BlockIndex = 1+(0:StepSize:size(wf,3)-BlockSize);
% BlockIndex(2,:) = BlockIndex(1,:)+BlockSize-1;
% M = size(wf,[1 2]);
% N = size(BlockIndex,2);
% WF = complex(zeros(size(wf,1),size(wf,2),BlockSize,N));
% for aa=1:N
%     WF(:,:,:,aa) = fftshift(fft(wf(:,:,BlockIndex(1,aa):BlockIndex(2,aa)),[],3),3);
% end
% freq = (-0.5:1/BlockSize:0.5-1/BlockSize)*RunLog.samplerate(3);
% WF = permute(WF,[1 2 4 3]);
% WF = reshape(WF,prod(M),N,BlockSize);
% Mask = reshape(WFInfo.Mask_WF,prod(M),1);
% Mask(isnan(Mask)) = 0;
% Mask = logical(Mask);
% WF = WF(Mask,:,:);
% clear aa wf;
% 
% % Sensor-FFT
% scanData = scanData(:,sensors)';
% Y = complex(zeros(size(scanData,1),BlockSize,N));
% for aa=1:N
%     Y(:,:,aa) = fftshift(fft(scanData(:,BlockIndex(1,aa):BlockIndex(2,aa)),[],2),2);
% end
% Y = permute(Y,[1 3 2]);
% clear aa scanData;
% 
% % LSE-SPOD
% WFao = complex(zeros(size(WF)));
% for aa=1:BlockSize
%     [psi,~] = eig(WF(:,:,aa)'*WF(:,:,aa));
%     psi = flip(psi,2);
%     phi = WF(:,:,aa)*psi;
%     L = (psi*Y(:,:,aa)')/(Y(:,:,aa)*Y(:,:,aa)');
%     psi_lse = L*Y(:,:,aa);
%     psi_ao = psi-psi_lse;
%     WFao(:,:,aa) = (psi_ao*phi')';
% end
% 
% % Unmask
% temp = WF;
% WF = complex(zeros(prod(M),N,BlockSize));
% WF(Mask,:,:) = temp;
% temp = WFao;
% WFao = complex(zeros(prod(M),N,BlockSize));
% WFao(Mask,:,:) = temp;
% clear temp;
% 
% % Dispersion
% WFao = reshape(WFao,M(1),M(2),N,BlockSize);
% WFao = squeeze(mean(fftshift(fftshift(fft(fft(WFao,[],1),[],2),1),2),3));
% WFao = (abs(WFao)).^2/numel(WFao)/prod(RunLog.samplerate);
% WF = reshape(WF,M(1),M(2),N,BlockSize);
% WF = squeeze(mean(fftshift(fftshift(fft(fft(WF,[],1),[],2),1),2),3));
% WF = (abs(WF)).^2/numel(WF)/prod(RunLog.samplerate);
% WFsize = size(WFao);
% frequency{1} = (-0.5:1/WFsize(1):0.5-1/WFsize(1))'*RunLog.samplerate(1);
% for aa=2:3
%     frequency{aa} = permute((-0.5:1/WFsize(aa):0.5-1/WFsize(aa))'*RunLog.samplerate(aa),aa:-1:1);
% end
% 
% 
% %%
% figure(1);
% log_range = -15;
% scolor = parula(2);
% patch(isosurface(frequency{2},frequency{1},frequency{3}(end/2+1:end),log10(WFao(:,:,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');
% daspect([1 1 50]);
% xlim(RunLog.samplerate(1)/2*[-1 1]);
% ylim(RunLog.samplerate(2)/2*[-1 1]);
% zlim(RunLog.samplerate(3)/2*[0 1]);
% material dull;
% camlight;
% grid on;
%     
% %%
% figure(2);
% loglog(squeeze(frequency{3}(end/2+1:end)),squeeze(mean(WFao(:,:,end/2+1:end),[1 2])),squeeze(frequency{3}(end/2+1:end)),squeeze(mean(WF(:,:,end/2+1:end),[1 2])));
% grid on;





