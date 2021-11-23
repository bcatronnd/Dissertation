close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
BlockSize = 2^10;
OverlapFactor = 0;
podDim = 3;

% Load Data
[wf,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);
load([directory testPoint '_CDAQ.mat'],'scanData');
% Window
Mask = wf(:,:,1);
Mask(isnan(Mask)) = 0;
Mask(Mask~=0) = 1;
MaskNaN = Mask;
MaskNaN(MaskNaN==0) = NaN;
windowS = createSpatialWindow(Mask);
windowT = reshape(hann(BlockSize),1,1,[]);
w = windowS.*windowT;
cw = 1/sqrt(sum(w.^2,'all')/numel(w));
clear WindowS WindowT;
wf(isnan(wf)) = 0;
% FFTs
StepSize = BlockSize*2^-OverlapFactor;
BlockIndex = 1+(0:StepSize:size(wf,3)-BlockSize);
BlockIndex(2,:) = BlockIndex(1,:)+BlockSize-1;
Q = zeros(size(wf,1),size(wf,2),BlockSize,size(BlockIndex,2));
for aa=1:size(BlockIndex,2)
    Q(:,:,:,aa) = (abs(cw*fftshift(fft(w.*wf(:,:,BlockIndex(1,aa):BlockIndex(2,aa)))))).^2/numel(w)/prod(RunLog.samplerate);
end
clear w wf;
WF = mean(Q,4);
frequency{1} = (-0.5:1/size(WF,1):0.5-1/size(WF,1))'*RunLog.samplerate(1);
frequency{2} = (-0.5:1/size(WF,2):0.5-1/size(WF,2))*RunLog.samplerate(2);
frequency{3} = permute((-0.5:1/size(WF,3):0.5-1/size(WF,3))'*RunLog.samplerate(3),[3 2 1]);
% WFm = mean(WF,[1 2]);

% SPOD
DimOrder = [1:podDim-1 podDim+1:3 4 podDim];
Q = permute(Q,DimOrder);
M = size(Q,[1 2]);
N = size(Q,3);
O = size(Q,4);
Q = reshape(Q,prod(M),N,O);

%%
Phi = zeros(prod(M),prod(M),O);
Psi = zeros(O,prod(M));
wb = waitbar(0);
for aa=1:O
    Qs = Q(:,:,aa)-mean(Q(:,:,aa),'all');
    WFs = reshape(WF(:,:,aa)-mean(WF(:,:,aa),'all'),prod(M),1);
    [phi,~] = eig(Qs*Qs'/(N-1));
    Phi(:,:,aa) = flip(phi,2);
    Psi(aa,:) = (WFs'*phi);
    waitbar(aa/O,wb);
end
close(wb);

%%
figure(1);
loglog(squeeze(frequency{3}(end/2+1:end)),abs(Psi(end/2+1:end,2)));
grid on;













% nModes = 16;
% Phi = complex(zeros(prod(M),N,size(Q,3)));
% Psi = complex(zeros(N,N,size(Q,3)));
% for aa=1:size(Q,3)
% %     [phi,~] = eig(Q(:,:,aa)*Q(:,:,aa)'/(N-1));
% %     Phi(:,1:nModes,aa) = phi(:,1:nModes);
% %     Phi = flip(Phi,2);
% %     Psi = (Q(:,:,aa)'*Q(:,:,aa))\Q(:,:,aa)'*Phi;
% %     Psi = Q(:,:,aa)'*Phi;
%     [Psi(:,:,aa),~] = eig(Q(:,:,aa)'*Q(:,:,aa));
% %     Psi = flip(Psi,2);
%     Phi(:,:,aa) = Q(:,:,aa)*Psi(:,:,aa);
% end
% 
% nPlots = 8;
% close all; figure(1); colormap(redblue);
% for aa=1:nPlots^2
% subplot(nPlots,nPlots,aa); surf(reshape(phi(:,aa),M),'linestyle','none'); view(2); axis image; caxis(max(abs(phi),[],'all')*[-1 1]);
% end


% subplot(1,3,2); surf(imag(reshape(Phi(:,1),M)),'linestyle','none'); view(2); axis image;
% subplot(1,3,3); surf(abs(reshape(Phi(:,1),M)),'linestyle','none'); view(2); axis image;
% subplot(3,3,4); surf(real(reshape(Phi(:,2),M)),'linestyle','none'); view(2); axis image;
% subplot(3,3,5); surf(imag(reshape(Phi(:,2),M)),'linestyle','none'); view(2); axis image;
% subplot(3,3,6); surf(abs(reshape(Phi(:,2),M)),'linestyle','none'); view(2); axis image;
% subplot(3,3,7); surf(real(reshape(Phi(:,3),M)),'linestyle','none'); view(2); axis image;
% subplot(3,3,8); surf(imag(reshape(Phi(:,3),M)),'linestyle','none'); view(2); axis image;
% subplot(3,3,9); surf(abs(reshape(Phi(:,3),M)),'linestyle','none'); view(2); axis image;
% size(Q)














