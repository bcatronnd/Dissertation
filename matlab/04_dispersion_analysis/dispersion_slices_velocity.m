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

fSlices = [5e3 1e4];
[~,fIndex] = min(abs(squeeze(freq{3})-fSlices));

uFracs = 1:-0.1:0.5;

f1 = figure(1);
clim = [-17 -7];
for aa=1:length(fSlices)
    subplot(1,2,aa);
    surf(freq{1},freq{2},WF(:,:,fIndex(aa)),'edgecolor','none');%,'facecolor','interp');
    view(2);
    xlim(fSlices(aa)./(RunLog.u*[20 0.5]));
    ylim(125*[-1 1]);
    caxis(clim);
%     axis equal tight;
    xix = fSlices(aa)./[RunLog.u+RunLog.c RunLog.u*uFracs];
    xticks(xix);
    xstrings = cell(7,1);
    xstrings(1) = {'u+c'};
    xstrings(1+(1:length(uFracs))) = split(num2str(uFracs,'%0.1f '));
    xticklabels(xstrings);
    xtickangle(60);
    
    xlabel('$u/U$','interpreter','latex');
    ylabel('$\xi_y\ (m^{-1})$','interpreter','latex');
    title(['$f=' num2str(fSlices(aa),'%0.0f') '\ Hz$'],'interpreter','latex');
    f1.Children(1).TickLabelInterpreter = 'latex';
    f1.Children(1).Layer = 'top';
    f1.Children(1).GridAlpha = 0.625;
    f1.Children(1).GridColor = [0 0 0];
    f1.Children(1).Position(2) = f1.Children(1).Position(2)+0.2;
    f1.Children(1).Position(4) = f1.Children(1).Position(4)-0.2;
end
colorbar('south');
f1.Children(1).Position(1) = f1.Children(3).Position(1)+0.025;
f1.Children(1).Position(2) = f1.Children(1).Position(2)-0.3;
f1.Children(1).Position(3) = f1.Children(2).Position(3)+f1.Children(2).Position(1)-f1.Children(3).Position(1)-0.05;
f1.Children(1).Label.String = '$S_{xx}$ ($\mu m^2/Hz/m^{-2}$)';
f1.Children(1).Label.Interpreter = 'latex';
f1.Children(1).TickLabelInterpreter = 'latex';
for aa=1:length(f1.Children(1).TickLabels)
    f1.Children(1).TickLabels{aa} = ['$10^{' f1.Children(1).TickLabels{aa} '}$'];
end
f1.Units = 'inches';
f1.Position = [1 1 5.5 3];
f1.Children(2).LineWidth = 0.8;
f1.Children(3).LineWidth = 0.8;

saveas(f1,'dispersion_slices_velocity.eps','epsc');
saveas(f1,'dispersion_slices_velocity.png','png');


