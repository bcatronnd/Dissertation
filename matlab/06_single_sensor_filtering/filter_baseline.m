close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6,'zernikeremoval',1:3);

BlockSize = 2^10;

WF = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));



[WF,freq] = simpleDispersion(WF,'BlockSize',BlockSize,'SampleRate',RunLog.samplerate);

% Unfiltered
log_range = -14;
clim = [-14 -7];
scolor = parula(2);
f1 = figure(1);
subplot(1,2,1);
patch(isocaps(freq{1},freq{2}(end/2+1:end),freq{3}(end/2+1:end),WF(end/2+1:end,:,end/2+1:end),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq{1},freq{2}(end/2+1:end),freq{3}(end/2+1:end),WF(end/2+1:end,:,end/2+1:end),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
grid on;
daspect([1 1 50]);
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[0 1]);
zlim(RunLog.samplerate(3)/2*[0 1]);
caxis(clim);
xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
zlabel('$f\ (Hz)$','Interpreter','Latex');
title('Unfiltered','Interpreter','Latex');
view(-45,15);
material dull;
camlight;
f1.Children(1).TickLabelInterpreter = 'latex';

% Filtering
for aa=1:size(WF,1)
    for bb=1:size(WF,2)
        WF2(aa,bb,:) = ipermute(baseline(permute(WF(aa,bb,:),[3 2 1])),[3 2 1]);
    end
end
subplot(1,2,2);
patch(isocaps(freq{1},freq{2}(end/2+1:end),freq{3}(end/2+1:end),WF2(end/2+1:end,:,end/2+1:end),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq{1},freq{2}(end/2+1:end),freq{3}(end/2+1:end),WF2(end/2+1:end,:,end/2+1:end),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
grid on;
daspect([1 1 50]);
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[0 1]);
zlim(RunLog.samplerate(3)/2*[0 1]);
caxis(clim);
xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
zlabel('$f\ (Hz)$','Interpreter','Latex');
title('Baseline','Interpreter','Latex');
view(-45,15);
material dull;
camlight;
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 5.5 3.75];

saveas(f1,'filter_baseline.eps','epsc');
saveas(f1,'filter_baseline.png','png');

disp(sum(10.^WF,'all'))
disp(sum(10.^WF2,'all'))
disp(sum(10.^WF2,'all')/sum(10.^WF,'all'))
