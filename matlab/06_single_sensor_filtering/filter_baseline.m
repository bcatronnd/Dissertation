close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6,'zernikeremoval',1:3);

BlockSize = 2^10;

WF = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));



[WF,freq] = simpleDispersion(WF,'BlockSize',BlockSize,'SampleRate',RunLog.samplerate);
% Filtering
WF2 = zeros(size(WF));
for aa=1:size(WF,1)
    for bb=1:size(WF,2)
        WF2(aa,bb,:) = ipermute(baseline(permute(WF(aa,bb,:),[3 2 1])),[3 2 1]);
    end
end

%%
close(findobj('type','figure','number',1));
f1 = figure(1);

log_range = -14;
clim = [-14 -7];
scolor = parula(2);
fPlot = [34 43];
PlotColors = linspecer(3);
faceAlpha = 0.2;

subplot(3,2,[1 3]);
patch(isocaps(freq{1}(end/2+1:end),freq{2}(end/2+1:end),freq{3}(end/2+1:end),(WF(end/2+1:end,end/2+1:end,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq{1}(end/2+1:end),freq{2}(end/2+1:end),freq{3}(end/2+1:end),(WF(end/2+1:end,end/2+1:end,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
patch(isocaps(freq{1}(end/2+1:end),freq{2}(1:end/2+1),freq{3}(end/2+1:end),(WF(1:end/2+1,end/2+1:end,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none','facealpha',faceAlpha);
patch(isosurface(freq{1}(end/2+1:end),freq{2}(1:end/2+1),freq{3}(end/2+1:end),(WF(1:end/2+1,end/2+1:end,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud','facealpha',faceAlpha);%,'specularstrength',0.375);
hold on;
plot3(freq{2}(fPlot(2))*[1 1],freq{1}(fPlot(1))*[1 1],RunLog.samplerate(3)/2*[0 1],'linewidth',1.25,'color',PlotColors(3,:));
grid on;
daspect([1 1 50]);
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[-1 1]);
zlim(RunLog.samplerate(3)/2*[0 1]);
caxis(clim);
% xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
zlabel('$f\ (Hz)$','Interpreter','Latex');
title('Unfiltered','Interpreter','Latex');
view(-45,15);
material dull;
camlight;
f1.Children(1).TickLabelInterpreter = 'latex';

subplot(3,2,[2 4]);
patch(isocaps(freq{1}(end/2+1:end),freq{2}(end/2+1:end),freq{3}(end/2+1:end),(WF2(end/2+1:end,end/2+1:end,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq{1}(end/2+1:end),freq{2}(end/2+1:end),freq{3}(end/2+1:end),(WF2(end/2+1:end,end/2+1:end,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
patch(isocaps(freq{1}(end/2+1:end),freq{2}(1:end/2+1),freq{3}(end/2+1:end),(WF2(1:end/2+1,end/2+1:end,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none','facealpha',faceAlpha);
patch(isosurface(freq{1}(end/2+1:end),freq{2}(1:end/2+1),freq{3}(end/2+1:end),(WF2(1:end/2+1,end/2+1:end,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud','facealpha',faceAlpha);%,'specularstrength',0.375);
hold on;
plot3(freq{2}(fPlot(2))*[1 1],freq{1}(fPlot(1))*[1 1],RunLog.samplerate(3)/2*[0 1],'linewidth',1.25,'color',PlotColors(3,:));
grid on;
daspect([1 1 50]);
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[-1 1]);
zlim(RunLog.samplerate(3)/2*[0 1]);
caxis(clim);
xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
% ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
zlabel('$f\ (Hz)$','Interpreter','Latex');
title('Baseline','Interpreter','Latex');
view(-45,15);
material dull;
camlight;
f1.Children(1).TickLabelInterpreter = 'latex';
colorbar('location','south');
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Children(1).Label.String = '$S_{xx}$ ($\mu m^2/Hz/m^{-2}$)';
f1.Children(1).Label.Interpreter = 'latex';



subplot(3,2,[5 6]);
loglog(squeeze(freq{3}),squeeze(10.^WF(fPlot(1),fPlot(2),:)),'linewidth',1.25,'color',PlotColors(1,:));
hold on;
loglog(squeeze(freq{3}),squeeze(10.^WF2(fPlot(1),fPlot(2),:)),'linewidth',1.25,'color',PlotColors(2,:));
grid on;
xlabel('$f\ (Hz)$','Interpreter','Latex');
ylabel('$S_{xx}\ (\mu m^2/Hz/m^{-2})$','Interpreter','Latex');
f1.Children(1).TickLabelInterpreter = 'latex';
legend('Unfiltered','Baseline','Interpreter','Latex');

f1.Units = 'inches';
f1.Position = [1 1 5.5 6];
f1.Children(5).Position(2) = f1.Children(5).Position(2)+0.05;
f1.Children(4).Position(2) = f1.Children(4).Position(2)+0.05;
f1.Children(3).Position(1) = f1.Children(5).Position(1);
f1.Children(3).Position(3) = f1.Children(4).Position(3)+f1.Children(4).Position(1)-f1.Children(5).Position(1);
f1.Children(3).Position(2) = f1.Children(3).Position(2)-0.15;
f1.Children(2).Position(2) = f1.Children(2).Position(2)-0.05;
for aa=1:length(f1.Children(3).Ticks)
    f1.Children(3).TickLabels{aa} = ['$10^{' num2str(f1.Children(3).Ticks(aa)) '}$'];
end
saveas(f1,'filter_baseline.eps','epsc');
saveas(f1,'filter_baseline.png','png');

disp(sum(10.^WF,'all'))
disp(sum(10.^WF2,'all'))
disp(1-sum(10.^WF2,'all')/sum(10.^WF,'all'))
