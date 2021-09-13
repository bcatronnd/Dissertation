close all; clc; clearvars;

load('synthetic_wavefront.mat');
nframes = 3;
sframe = 1000;
stepFrame = 5;

clim = max(abs(wf.wf),[],'all','omitnan');


f1 = figure(1);
colormap(redblue);
for aa=1:nframes
    subplot(1,nframes,aa);
    surf(wf.wf(:,:,sframe-1+aa*stepFrame),'linestyle','none');
    view(2);
    axis image off tight;
    caxis(clim*[-1 1]);
    sgtitle('Total Signal','interpreter','latex','fontsize',15);
end 
f1.Units = 'inches';
f1.Position = [1+0.5*get(gcf,'number')*[1 1] 6 2];
f2 = figure(2);
colormap(redblue);
for aa=1:nframes
    subplot(1,nframes,aa);
    surf(wf.AO(:,:,sframe-1+aa*stepFrame),'linestyle','none');
    view(2);
    axis image off tight;
    caxis(clim*[-1 1]);
    sgtitle('Aero-Optical Signal','interpreter','latex','fontsize',15);
end
f2.Units = 'inches';
f2.Position = [1+0.5*get(gcf,'number')*[1 1] 6 2];

saveas(f1,'synthetic_frames_total.eps','epsc');
saveas(f2,'synthetic_frames_ao.eps','epsc');




