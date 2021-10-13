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


%%
fLines = [5e3 7.5e3 1e4 1.25e4 1.5e4 1.75e4];
yLines = zeros(size(fLines));

[~,fIndex] = min(abs(squeeze(freq{3})-fLines));
[~,yIndex] = min(abs(squeeze(freq{1})-yLines));
for aa=1:length(fLines)
    fString{aa} = ['f=' num2str(fLines(aa)) ' Hz'];
end
close all;
f1 = figure(1);
plot(fLines(1)./freq{2}/RunLog.u,WF(yIndex(1),:,fIndex(1)));
hold on;
for aa=2:length(fLines)
    plot(fLines(aa)./freq{2}/RunLog.u,WF(yIndex(aa),:,fIndex(aa)));
end
grid on;
xlim([0.6 1.4]);
% ylim([-17 -11]);
% xticks([0.7 0.8 0.9 1.0 1.1 1.2]);
xlabel('$u_{assumed}/u_\infty$','interpreter','latex');
ylabel('$S_{xx}$ ($\mu m^2/Hz/m^{-1}/m^{-1}$)','interpreter','latex');
% legend(fString,'interpreter','latex','location','southoutside','orientation','horizontal');
legend(fString,'interpreter','latex','location','northeast');
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Children(2).XTick = 0.6:0.2:1.4;
for aa=1:length(f1.Children(2).YTickLabels)
    f1.Children(2).YTickLabels{aa} = ['$10^{' f1.Children(2).YTickLabels{aa} '}$'];
end
f1.Units = 'inches';
f1.Position = [1 1 5.5 3.25];

saveas(f1,'dispersion_speed.eps','epsc');