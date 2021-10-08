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

f1 = figure(1);
fPoints = [0 2.5e4];
clim = [-17 -7];
subplot(2,1,1);
plot(fPoints,fPoints/RunLog.u,'k-',fPoints,fPoints/(RunLog.u+RunLog.c),'k--',fPoints,fPoints/(RunLog.u-RunLog.c),'k-.');
hold on;
surf(squeeze(freq{3}),squeeze(freq{2}),squeeze(WF(end/2+1,:,:)),'linestyle','none','facecolor','interp');
view(2);
colorbar;
xlim([0 2e4]);
ylim(125*[-1 1]);
caxis(clim);
ylabel('X-Spatial Frequency, $\xi_x$ ($1/m$)','interpreter','latex');
legend('$u$','$u+c$','$u-c$','interpreter','latex','location','northwest');
title('Horizontal Moving Disturbances','interpreter','latex');
subplot(2,1,2);
plot(fPoints,fPoints/(+RunLog.c),'k--',fPoints,fPoints/(-RunLog.c),'k-.');
hold on;
surf(squeeze(freq{3}),squeeze(freq{1}),squeeze(WF(:,end/2+1,:)),'linestyle','none','facecolor','interp');
view(2);
xlim([0 2e4]);
ylim(125*[-1 1]);
caxis(clim);
xlabel('Temporal Frequency, $f$ ($Hz$)','interpreter','latex');
ylabel('Y-Spatial Frequency, $\xi_y$ ($1/m$)','interpreter','latex');
legend('$+c$','$-c$','interpreter','latex','location','northwest');
title('Vertical Moving Disturbances','interpreter','latex');
f1.Children(5).Position(3) = f1.Children(5).Position(3);
f1.Children(2).Position(3) = f1.Children(5).Position(3);
f1.Children(4).Position(2) = f1.Children(2).Position(2);
f1.Children(4).Position(4) = f1.Children(5).Position(4)+f1.Children(5).Position(2)-f1.Children(2).Position(2);
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Children(4).TickLabelInterpreter = 'latex';
f1.Children(5).TickLabelInterpreter = 'latex';
f1.Children(2).Layer = 'top';
f1.Children(5).Layer = 'top';
f1.Children(4).Label.String = '$S_{xx}$ ($\mu m^2/Hz/m^{-1}/m^{-1}$)';
f1.Children(4).Label.Interpreter = 'latex';
for aa=1:length(f1.Children(4).TickLabels)
    f1.Children(4).TickLabels{aa} = ['$10^{' f1.Children(4).TickLabels{aa} '}$'];
end
f1.Units = 'inches';
f1.Position = [1 1 5.5 5];

saveas(f1,'dispersion_xy.eps','epsc');








