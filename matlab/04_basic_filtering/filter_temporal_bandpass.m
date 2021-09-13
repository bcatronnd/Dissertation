close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/ProcessedWavefronts/';
filename = '20180104022';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory filename '_WF.mat'],'Permute',[2 1 3],'Scale',1e6,'ZernikeRemoval',1:3);
WFp = WFfilter(WF,'time-bandpass',(532+10*[-1 1])/RunLog.samplerate(3));
WFs = WFfilter(WF,'time-bandstop',(532+10*[-1 1])/RunLog.samplerate(3));

start_frame = 1300;
num_frames = 4;
step_frame = 10;
clim = max(abs(WF(:,:,start_frame+(0:step_frame*num_frames-1))),[],'all','omitnan')*[-1 1];
f1 = figure(1);
colormap(redblue);
for aa=1:num_frames
    subplot(num_frames,2,2*(aa-1)+1);
    surf(WFs(:,:,start_frame+(aa-1)*step_frame),'linestyle','none');
    view(2);
    axis image off tight;
    caxis(clim);
    
    subplot(num_frames,2,2*aa);
    surf(WFp(:,:,start_frame+(aa-1)*step_frame),'linestyle','none');
    view(2);
    axis image off tight;
    caxis(clim);
end
f1.Units = 'inches';
f1.Position = [1 1 5.5 8];

saveas(f1,'filter_temporal_bandpass.eps','epsc');

