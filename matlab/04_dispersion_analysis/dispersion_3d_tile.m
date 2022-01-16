close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);

BlockSize = 2^10;
WF = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));
[WF,freq] = simpleDispersion(WF,'BlockSize',BlockSize,'SampleRate',RunLog.samplerate);

log_range = -14;
disp(['Total Energy Inside Surface: ' num2str(sum(10.^WF(WF>log_range),'all')/sum(10.^WF,'all')*100,'%0.1f') '%']);

%%
close(findobj('type','figure','number',1));
f1 = figure(1);

scolor = parula(2);
log_range = -14.5;
fRange = [0.4 0.6];

[~,fIndex1] = min(abs(squeeze(freq{3}/RunLog.samplerate(3))-fRange(1)));
[~,fIndex2] = min(abs(squeeze((freq{3}+RunLog.samplerate(3))/RunLog.samplerate(3))-fRange(2)));

patch(isocaps(freq{1}/RunLog.samplerate(1),freq{2}/RunLog.samplerate(2),freq{3}(fIndex1:end)/RunLog.samplerate(3),WF(:,:,fIndex1:end),log_range(1),'all'),'facecolor','blue','edgecolor','none','facelighting','none');
patch(isosurface(freq{1}/RunLog.samplerate(1),freq{2}/RunLog.samplerate(2),freq{3}(fIndex1:end)/RunLog.samplerate(3),WF(:,:,fIndex1:end),log_range(1)),'edgecolor','none','facecolor','blue','facelighting','gouraud');%,'specularstrength',0.375);

patch(isocaps(freq{1}/RunLog.samplerate(1),freq{2}/RunLog.samplerate(2),freq{3}(1:fIndex2)/RunLog.samplerate(3)+1,WF(:,:,1:fIndex2),log_range(1),'all'),'facecolor','red','edgecolor','none','facelighting','none');
patch(isosurface(freq{1}/RunLog.samplerate(1),freq{2}/RunLog.samplerate(2),freq{3}(1:fIndex2)/RunLog.samplerate(3)+1,WF(:,:,1:fIndex2),log_range(1)),'edgecolor','none','facecolor','red','facelighting','gouraud');%,'specularstrength',0.375);

grid on;
daspect([1 1 0.2]);
xlim(1/2*[-1 1]);
ylim(1/2*[-1 1]);
zlim(fRange);
xlabel('$\xi_x/\xi_{s,x}$','Interpreter','Latex');
ylabel('$\xi_y/\xi_{s,y}$','Interpreter','Latex');
zlabel('$f/f_{s}$','Interpreter','Latex');
view(-45,15);
material dull;
camlight;
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 5.5 4.5];
saveas(f1,'dispersion_3d_tile.eps','epsc');

