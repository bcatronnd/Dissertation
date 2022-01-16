close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);

BlockSize = 2^12;
WF = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));
[WF,freq] = simpleDispersion(WF,'BlockSize',BlockSize,'SampleRate',RunLog.samplerate);

blockageRatio = 34.462/36^2;
RunLog.u = RunLog.u*(1+blockageRatio);
% RunLog.u = 0.8*RunLog.u;
lowerRatio = 0.7;

% Super Sample Slices
xSlice = squeeze(WF(end/2+1,:,:));
xSlice = repmat(xSlice,3,3);
ySlice = squeeze(WF(:,end/2+1,:));
ySlice = repmat(ySlice,3,3);
tFreq = diff(freq{3}(1:2))*(-size(xSlice,2)/2:size(xSlice,2)/2-1);
xFreq = diff(freq{2}(1:2))*(-size(xSlice,1)/2:size(xSlice,1)/2-1);
yFreq = diff(freq{1}(1:2))*(-size(ySlice,1)/2:size(ySlice,1)/2-1);

%%
close(findobj('type','figure','number',1));
f1 = figure(1);

clim = [-17 -7];
fPoints = [-1.5 1.5];

plot(fPoints,fPoints/RunLog.u*RunLog.samplerate(3)/RunLog.samplerate(2),'k-',fPoints,fPoints/(RunLog.u+RunLog.c)*RunLog.samplerate(3)/RunLog.samplerate(2),'k--',fPoints,fPoints/(RunLog.u-RunLog.c)*RunLog.samplerate(3)/RunLog.samplerate(2),'k-.');
hold on;
surf(tFreq/RunLog.samplerate(3),xFreq/RunLog.samplerate(2),xSlice,'linestyle','none');
view(2);
grid on;
plot([-1.5 1.5 NaN -1.5 1.5 NaN -1.5 1.5 NaN -1.5 1.5],[1.5 1.5 NaN 0.5 0.5 NaN -0.5 -0.5 NaN -1.5 -1.5],'k--','linewidth',1.5);
plot([1.5 1.5 NaN 0.5 0.5 NaN -0.5 -0.5 NaN -1.5 -1.5],[-1.5 1.5 NaN -1.5 1.5 NaN -1.5 1.5 NaN -1.5 1.5],'k--','linewidth',1.5);
plot([0.5 0.5 -0.5 -0.5 0.5],[0.5 -0.5 -0.5 0.5 0.5],'k-','linewidth',2);
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
caxis(clim);
xlabel('Temporal Frequency, $f/f_s$','interpreter','latex');
ylabel('X-Spatial Frequency, $\xi_x/\xi_{x,s}$','interpreter','latex');
title('Horizontal Moving Disturbances','interpreter','latex');
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Children(1).Layer = 'top';
legend('$u$','$u+c$','$u-c$','interpreter','latex','location','northwest');

f1.Units = 'inches';
f1.Position = [1 1 5.5 3.5];

saveas(f1,'dispersion_tiling.eps','epsc');
saveas(f1,'dispersion_tiling.png','png');

