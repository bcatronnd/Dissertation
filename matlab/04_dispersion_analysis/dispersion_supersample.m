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
close all;

clim = [-17 -7];
fPoints = [0 4e4];
f1 = figure(1);
subplot(2,1,1);
plot(fPoints,fPoints/RunLog.u,'k-',fPoints,fPoints/(RunLog.u+RunLog.c),'k--',fPoints,fPoints/(RunLog.u-RunLog.c),'k-.');
hold on;
surf(tFreq,xFreq,xSlice,'edgecolor','none');
view(2);
plot(max(freq{3},[],'all')*[0 1 1 0],max(freq{2},[],'all')*[-1 -1 1 1],'k-','linewidth',2);
colorbar;
grid on;
xlim([0 4e4]);
ylim(200*[-1 1]);
caxis(clim);
ylabel('X-Spatial Frequency, $\xi_x$ ($m^{-1}$)','interpreter','latex');
title('Horizontal Moving Disturbances','interpreter','latex');

subplot(2,1,2);
plot(fPoints,fPoints/(+RunLog.c),'k--',fPoints,fPoints/(-RunLog.c),'k-.');
hold on;
surf(tFreq,yFreq,ySlice,'edgecolor','none');
view(2);
plot(max(freq{3},[],'all')*[0 1 1 0],max(freq{2},[],'all')*[-1 -1 1 1],'k-','linewidth',2);
grid on;
xlim([0 4e4]);
ylim(200*[-1 1]);
caxis(clim);
xlabel('Temporal Frequency, $f$ ($Hz$)','interpreter','latex');
ylabel('Y-Spatial Frequency, $\xi_y$ ($1/m$)','interpreter','latex');
title('Vertical Moving Disturbances','interpreter','latex');

f1.Children(3).Position(3) = f1.Children(3).Position(3);
f1.Children(1).Position(3) = f1.Children(3).Position(3);
f1.Children(2).Position(2) = f1.Children(1).Position(2);
f1.Children(2).Position(4) = f1.Children(3).Position(4)+f1.Children(3).Position(2)-f1.Children(1).Position(2);
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Children(3).TickLabelInterpreter = 'latex';
f1.Children(1).Layer = 'top';
f1.Children(3).Layer = 'top';
f1.Children(2).Label.String = '$S_{xx}$ ($\mu m^2/Hz/m^{-1}/m^{-1}$)';
f1.Children(2).Label.Interpreter = 'latex';
for aa=1:length(f1.Children(2).TickLabels)
    f1.Children(2).TickLabels{aa} = ['$10^{' f1.Children(2).TickLabels{aa} '}$'];
end
f1.Units = 'inches';
f1.Position = [1 1 5.5 5];

% saveas(f1,'dispersion_supersample.eps','epsc');