close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6,'zernikeremoval',1:3);

% WF = WF+WFInfo.WFmean;
BlockSize = 2^10;

WF = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));
[WF,freq] = simpleDispersion(WF,'BlockSize',BlockSize,'SampleRate',RunLog.samplerate);



%%
close all;
f1 = figure(1);
log_range = -14;
clim = [-14 -8];
scolor = parula(2);

patch(isocaps(freq{1},freq{2}(end/2+1:end),freq{3}(end/2+1:end),WF(end/2+1:end,:,end/2+1:end),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq{1},freq{2}(end/2+1:end),freq{3}(end/2+1:end),WF(end/2+1:end,:,end/2+1:end),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
patch(isosurface(freq{1},freq{2}(1:end/2+1),freq{3}(end/2+1:end),WF(1:end/2+1,:,end/2+1:end),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud','facealpha',0.1);%,'specularstrength',0.375);

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
colorbar;
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Children(1).Label.String = '$S_{xx}$ ($\mu m^2/Hz/m^{-2}$)';
f1.Children(1).Label.Interpreter = 'latex';
for aa=1:length(f1.Children(1).TickLabels)
    f1.Children(1).TickLabels{aa} = ['$10^{' f1.Children(1).TickLabels{aa} '}$'];
end
f1.Units = 'inches';
f1.Position = [1 1 5.5 5];

saveas(f1,'dispersion_isosurface.eps','epsc');

