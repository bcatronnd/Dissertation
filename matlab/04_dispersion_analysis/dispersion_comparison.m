close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[wf,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6,'ZernikeRemoval',1:3);
WF = wf;

blockSize = 2^12;
sampleRate = [277.7778 49020];

% 2-D Dispersion
wf = squeeze(wf(round(size(wf,1)/2),:,:));
wf = wf(:,1:floor(size(wf,2)/blockSize)*blockSize);
wf = reshape(wf,size(wf,1),blockSize,size(wf,2)/blockSize);
window = hann(size(wf,1)).*hann(size(wf,2))';
cw = 1/sqrt(sum(window.^2,'all')/blockSize);
sxx2 = mean(cw*fftshift(fftshift((abs(fft2(wf.*window))).^2,1),2)/prod(size(wf,[1 2]))/prod(sampleRate),3);
for aa=1:2
    freq2{aa} = (-0.5:1/size(wf,aa):0.5-1/size(wf,aa))*sampleRate(aa);
end

% 3-D Dispersion
WF = WF(:,:,1:floor(size(WF,3)/blockSize)*blockSize);
WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));
[WF,freq] = simpleDispersion(WF,'BlockSize',blockSize,'SampleRate',RunLog.samplerate);


blockageRatio = 34.462/36^2;
RunLog.u = RunLog.u*(1+blockageRatio);
% RunLog.u = 0.8*RunLog.u;
lowerRatio = 0.7;

%%
close all;
fPoints = [0 2.5e4];
clim = [-17 -7];
f1 = figure(1);
subplot(3,1,1);
plot(fPoints,fPoints/RunLog.u,'k-',fPoints,fPoints/(RunLog.u+RunLog.c),'k--',fPoints,fPoints/(RunLog.u-RunLog.c),'k-.',fPoints,fPoints/RunLog.u/lowerRatio,'k:');
hold on;
surf(squeeze(freq{3}),squeeze(freq{2}),squeeze(WF(end/2+1,:,:)),'linestyle','none','facecolor','interp');
view(2);
grid on;
colorbar;
xlim([0 2e4]);
ylim(125*[-1 1]);
caxis(clim);
title('3-D Spectrum','interpreter','latex');
colorbar;
legend('$u$','$u+c$','$u-c$','interpreter','latex','location','northwest');
subplot(3,1,2);
plot(fPoints,fPoints/RunLog.u,'k-',fPoints,fPoints/(RunLog.u+RunLog.c),'k--',fPoints,fPoints/(RunLog.u-RunLog.c),'k-.',fPoints,fPoints/RunLog.u/lowerRatio,'k:');
hold on;
surf(squeeze(freq{3}),squeeze(freq{2}),squeeze(log10(trapz(freq{1},10.^WF,1))),'linestyle','none','facecolor','interp');
view(2);
grid on;
xlim([0 2e4]);
ylim(125*[-1 1]);
caxis(clim);
ylabel('Horizontal Spatial Frequency, $\xi_x$ ($m^{-1}$)','interpreter','latex');
title('3-D Spectrum - Integrated','interpreter','latex');
colorbar;
subplot(3,1,3);
plot(fPoints,fPoints/RunLog.u,'k-',fPoints,fPoints/(RunLog.u+RunLog.c),'k--',fPoints,fPoints/(RunLog.u-RunLog.c),'k-.',fPoints,fPoints/RunLog.u/lowerRatio,'k:');
hold on;
surf(freq2{2},freq2{1},log10(sxx2),'linestyle','none','facecolor','interp');
view(2);
grid on;
xlim([0 2e4]);
ylim(125*[-1 1]);
caxis(clim);
xlabel('Temporal Frequency, $f$ ($Hz$)','interpreter','latex');
title('2-D Spectrum','interpreter','latex');

f1.Units = 'inches';
f1.Position = [1 1 5.5 7.5];
f1.Children(1).Position(3) = f1.Children(6).Position(3);
f1.Children(3).Position(3) = f1.Children(6).Position(3);
f1.Children(2).Position(2) = f1.Children(1).Position(2);
f1.Children(2).Position(1) = f1.Children(5).Position(1);
f1.Children(2).Position(4) = f1.Children(3).Position(4)+f1.Children(3).Position(2)-f1.Children(1).Position(2);
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Children(3).TickLabelInterpreter = 'latex';
f1.Children(5).TickLabelInterpreter = 'latex';
f1.Children(6).TickLabelInterpreter = 'latex';
f1.Children(1).Layer = 'top';
f1.Children(3).Layer = 'top';
f1.Children(6).Layer = 'top';
f1.Children(2).Label.String = '$S_{xx}$ ($\mu m^2/Hz/m^{-1}$)';
f1.Children(2).Label.Interpreter = 'latex';
f1.Children(5).Label.String = '$S_{xx}$ ($\mu m^2/Hz/m^{-2}$)';
f1.Children(5).Label.Interpreter = 'latex';
for aa=1:length(f1.Children(2).TickLabels)
    f1.Children(2).TickLabels{aa} = ['$10^{' f1.Children(2).TickLabels{aa} '}$'];
end
for aa=1:length(f1.Children(5).TickLabels)
    f1.Children(5).TickLabels{aa} = ['$10^{' f1.Children(5).TickLabels{aa} '}$'];
end


saveas(f1,'dispersion_comparison.eps','epsc');

