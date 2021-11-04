close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);
BlockSize = 2^12;
WF = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));
[WF,freq] = simpleDispersion(WF,'BlockSize',BlockSize,'SampleRate',RunLog.samplerate,'output','log');

offset = mean(WF,'all');
[Phi,a,M,N] = computePOD(WF-offset,3);

%%
close all;
PODmodes = {1:2^4};

for aa=1:length(PODmodes)
    wf = Phi(:,PODmodes{aa})*a(:,PODmodes{aa})';
    wf = reshape(wf,[M N])+offset;
    
    log_range = -14;
    subplot(1,2,1);
    scolor = parula(2);
    patch(isocaps(freq{1},freq{2},freq{3}(end/2+1:end),WF(:,:,end/2+1:end),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
    patch(isosurface(freq{1},freq{2},freq{3},WF,log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
    grid on;
    daspect([1 1 50]);
    xlim(RunLog.samplerate(1)/2*[-1 1]);
    ylim(RunLog.samplerate(2)/2*[-1 1]);
    zlim(RunLog.samplerate(3)/2*[0 1]);
    xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
    ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
    zlabel('$f\ (Hz)$','Interpreter','Latex');
    view(3)
    material dull;
    camlight;
    f1.Children(1).TickLabelInterpreter = 'latex';
    
    subplot(1,2,2);
    scolor = parula(2);
    patch(isocaps(freq{1},freq{2},freq{3}(end/2+1:end),wf(:,:,end/2+1:end),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
    patch(isosurface(freq{1},freq{2},freq{3},wf,log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
    grid on;
    daspect([1 1 50]);
    xlim(RunLog.samplerate(1)/2*[-1 1]);
    ylim(RunLog.samplerate(2)/2*[-1 1]);
    zlim(RunLog.samplerate(3)/2*[0 1]);
    xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
    ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
    zlabel('$f\ (Hz)$','Interpreter','Latex');
    view(3)
    material dull;
    camlight;
    f1.Children(1).TickLabelInterpreter = 'latex';
end


