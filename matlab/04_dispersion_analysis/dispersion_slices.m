close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);

BlockSize = 2^12;
WF = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));
[WF,freq] = simpleDispersion(WF,'BlockSize',BlockSize,'SampleRate',RunLog.samplerate);

%%
close all;

fSlices = [0 RunLog.fanfreqency/60*32*1197/60 3*RunLog.fanfreqency/60*32*1197/60 5e3 1e4 2e4];
[~,fIndex] = min(abs(squeeze(freq{3})-fSlices));


f1 = figure(1);
clim = [-17 -7];
for aa=1:length(fSlices)
    subplot(3,2,aa);
    surf(freq{1},freq{2},WF(:,:,fIndex(aa)),'edgecolor','none','facecolor','interp');
    view(2);
    xlim(125*[-1 1]);
    ylim(125*[-1 1]);
    caxis(clim);
    axis equal tight;
    xlabel('$\xi_x\ (m^{-1})$','interpreter','latex');
    ylabel('$\xi_y\ (m^{-1})$','interpreter','latex');
    title(['$f=' num2str(fSlices(aa),'%0.1f') '\ Hz$'],'interpreter','latex');
    f1.Children(1).TickLabelInterpreter = 'latex';
    f1.Children(1).Layer = 'top';
    f1.Children(1).Position(2) = f1.Children(1).Position(2)+0.0375;
end
colorbar('south');
f1.Children(1).Position(1) = f1.Children(3).Position(1)+0.025;
f1.Children(1).Position(2) = f1.Children(1).Position(2)-0.125;
f1.Children(1).Position(3) = f1.Children(2).Position(3)+f1.Children(2).Position(1)-f1.Children(3).Position(1)-0.05;
f1.Children(1).Label.String = '$S_{xx}$ ($\mu m^2/Hz/m^{-1}/m^{-1}$)';
f1.Children(1).Label.Interpreter = 'latex';
f1.Children(1).TickLabelInterpreter = 'latex';
for aa=1:length(f1.Children(1).TickLabels)
    f1.Children(1).TickLabels{aa} = ['$10^{' f1.Children(1).TickLabels{aa} '}$'];
end
f1.Units = 'inches';
f1.Position = [1 1 5.5 7.5];

saveas(f1,'dispersion_slices.eps','epsc');



