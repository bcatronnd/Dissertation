close all; clc; clearvars;

directory = '../../data/';
testPoint = {'20210901007' '20210901003' '20210831023' '20210901022'};
angle = [75 90 105 120];

BlockSize = 2^10;
log_range = -14;
for aa=1:length(testPoint)
%     close all;
    [WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint{aa} '_WF.mat'],'Scale',1e6);
    WF = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
    WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
    WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));
    [WF,freq] = simpleDispersion(WF,'BlockSize',BlockSize,'SampleRate',RunLog.samplerate);
    disp(['Total Energy Inside Surface: ' num2str(sum(10.^WF(WF>log_range),'all')/sum(10.^WF,'all')*100,'%0.1f') '%']);
    f1 = figure(aa);
    subplot(1,2,1);
    scolor = parula(2);
    patch(isocaps(freq{1},freq{2},freq{3}(end/2+1:end),WF(:,:,end/2+1:end),log_range(1),'all'),'facecolor','blue','edgecolor','none','facelighting','none');
    patch(isosurface(freq{1},freq{2},freq{3},WF,log_range(1)),'edgecolor','none','facecolor','blue','facelighting','gouraud');%,'specularstrength',0.375);
    grid on;
    daspect([1 1 100]);
    xlim(RunLog.samplerate(1)/2*[-1 1]);
    ylim(RunLog.samplerate(2)/2*[-1 1]);
    zlim(RunLog.samplerate(3)/2*[0 1]);
    xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
    ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
    zlabel('$f\ (Hz)$','Interpreter','Latex');
    view([0 0]);
    material dull;
    camlight;
    f1.Children(1).TickLabelInterpreter = 'latex';
    subplot(1,2,2);
    scolor = parula(2);
    patch(isocaps(freq{1},freq{2},freq{3}(end/2+1:end),WF(:,:,end/2+1:end),log_range(1),'all'),'facecolor','blue','edgecolor','none','facelighting','none');
    patch(isosurface(freq{1},freq{2},freq{3},WF,log_range(1)),'edgecolor','none','facecolor','blue','facelighting','gouraud');%,'specularstrength',0.375);
    grid on;
    daspect([1 1 100]);
    xlim(RunLog.samplerate(1)/2*[-1 1]);
    ylim(RunLog.samplerate(2)/2*[-1 1]);
    zlim(RunLog.samplerate(3)/2*[0 1]);
    xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
    ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
    zlabel('$f\ (Hz)$','Interpreter','Latex');
    view([-90 0]);
    material dull;
    camlight;
    f1.Children(1).TickLabelInterpreter = 'latex';
    sgtitle(['$\gamma$ = ' num2str(angle(aa),'%0.0f') '$^{\circ}$'],'interpreter','latex');
    f1.Units = 'inches';
    f1.Position = [1 1 5.5 2.5];
    
    saveas(f1,['dispersion_angle_' num2str(angle(aa),'%0.0f') '.eps'],'epsc');
    
end














