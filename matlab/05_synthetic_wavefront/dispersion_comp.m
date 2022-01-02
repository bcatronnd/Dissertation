close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/ProcessedWavefronts/';
testPoint = '20180104022';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6,'ZernikeRemoval',1:3);
load('synthetic_wavefront.mat');

blockSize = 2^10;

% 3-D Dispersion
WF = WF(:,:,1:floor(size(WF,3)/blockSize)*blockSize);
WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));
[WF,freq] = simpleDispersion(WF,'BlockSize',blockSize,'SampleRate',RunLog.samplerate);
[WFs,freqs] = simpleDispersion(wf.wf,'BlockSize',blockSize,'SampleRate',wf.sampleRate);







%%
close(findobj('type','figure','number',1));
f1 = figure(1);

log_range = -12.5;
clim = [log_range -7];
scolor = parula(2);
faceAlpha = 0.2;

patch(isocaps(freq{1},freq{2}(end/2+1:end),freq{3}(end/2+1:end),(WF(end/2+1:end,:,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq{1},freq{2}(end/2+1:end),freq{3}(end/2+1:end),(WF(end/2+1:end,:,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
patch(isocaps(freq{1},freq{2}(1:end/2+1),freq{3}(end/2+1:end),(WF(1:end/2+1,:,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none','facealpha',faceAlpha);
patch(isosurface(freq{1},freq{2}(1:end/2+1),freq{3}(end/2+1:end),(WF(1:end/2+1,:,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud','facealpha',faceAlpha);%,'specularstrength',0.375);
grid on;
daspect([1 1 65]);
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[-1 1]);
zlim(RunLog.samplerate(3)/2*[0 1]);
caxis(clim);
xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
zlabel('$f\ (Hz)$','Interpreter','Latex');
view(-45,15);
material dull;
camlight;
f1.Children(1).TickLabelInterpreter = 'latex';
colorbar('location','south');
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Children(1).Label.String = '$S_{xx}$ ($\mu m^2/Hz/m^{-2}$)';
f1.Children(1).Label.Interpreter = 'latex';
f1.Children(2).Position(2) = f1.Children(2).Position(2)+0.075;
f1.Children(1).Position(2) = f1.Children(1).Position(2)-0.175;
f1.Units = 'inches';
f1.Position = [1 1 5.5 6];
for aa=1:length(f1.Children(1).Ticks)
    f1.Children(1).TickLabels{aa} = ['$10^{' num2str(f1.Children(1).Ticks(aa)) '}$'];
end
saveas(f1,'dispersion_comp_real.eps','epsc');


%%
close(findobj('type','figure','number',2));
f2 = figure(2);

log_range = -14.5;
clim = [log_range -8];

% patch(isocaps(freqs{1},freqs{2}(end/2+1:end),freqs{3}(end/2+1:end),WFs(end/2+1:end,:,end/2+1:end),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
% patch(isosurface(freqs{1},freqs{2}(end/2+1:end),freqs{3}(end/2+1:end),WFs(end/2+1:end,:,end/2+1:end),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
patch(isocaps(freqs{1},freqs{2}(end/2+1:end),freqs{3}(end/2+1:end),(WFs(end/2+1:end,:,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freqs{1},freqs{2}(end/2+1:end),freqs{3}(end/2+1:end),(WFs(end/2+1:end,:,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
patch(isocaps(freqs{1},freqs{2}(1:end/2+1),freqs{3}(end/2+1:end),(WFs(1:end/2+1,:,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none','facealpha',faceAlpha);
patch(isosurface(freqs{1},freqs{2}(1:end/2+1),freqs{3}(end/2+1:end),(WFs(1:end/2+1,:,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud','facealpha',faceAlpha);%,'specularstrength',0.375);

grid on;
daspect([1 1 65]);
xlim(wf.sampleRate(1)/2*[-1 1]);
ylim(wf.sampleRate(2)/2*[-1 1]);
zlim(wf.sampleRate(3)/2*[0 1]);
caxis(clim);
xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
zlabel('$f\ (Hz)$','Interpreter','Latex');
view(-45,15);
material dull;
camlight;
f2.Children(1).TickLabelInterpreter = 'latex';
colorbar('location','south');
f2.Children(1).TickLabelInterpreter = 'latex';
f2.Children(1).Label.String = '$S_{xx}$ ($\mu m^2/Hz/m^{-2}$)';
f2.Children(1).Label.Interpreter = 'latex';
f2.Children(2).Position(2) = f2.Children(2).Position(2)+0.075;
f2.Children(1).Position(2) = f2.Children(1).Position(2)-0.175;
f2.Units = 'inches';
f2.Position = [1 1 5.5 6];
for aa=1:length(f2.Children(1).Ticks)
    f2.Children(1).TickLabels{aa} = ['$10^{' num2str(f2.Children(1).Ticks(aa)) '}$'];
end
saveas(f2,'dispersion_comp_synthetic.eps','epsc');

