close all; clc; clearvars;

directory = '../../data/';
testPoint = '20210901003';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6,'zernikeremoval',1:3);

disp(['Original OPDrms: ' num2str(mean(nanrms(reshape(WF,size(WF,1)*size(WF,2),size(WF,3)))),'%0.4f') ' ' char(956) 'm']);

% Options
BlockSize = 2^10;
vFilter = [[0.9 0]*RunLog.u*RunLog.samplerate(1)/RunLog.samplerate(3),0.07,20,0];
faceAlpha = 0.2;

% mask = ones(2.^nextpow2(size(WF,1,2)));
mask = WFInfo.Mask_WF;
mask = cat(1,mask,NaN*zeros(2^nextpow2(size(mask,1))-size(mask,1),size(mask,2)));
mask = cat(2,mask,NaN*zeros(size(mask,1),2^nextpow2(size(mask,2))-size(mask,2)));

WF = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));
[WF,freq] = simpleDispersion(WF,'BlockSize',BlockSize,'SampleRate',RunLog.samplerate,'output','lin');

nSamples = size(WF);
phase = pi*(2*rand(nSamples)-1);
phase(nSamples(1)/2+1,nSamples(2)/2+1,nSamples(3)/2+1) = 0;
phase(2:nSamples(1),2:nSamples(2),nSamples(3)/2+2:nSamples(3)) = -flip(flip(flip(phase(2:nSamples(1),2:nSamples(2),2:nSamples(3)/2),1),2),3);
wf = real(ifftn(ifftshift(sqrt((WF)*prod(RunLog.samplerate)*numel(WF)).*exp(1i*phase)))).*mask;
disp(['Original OPDrms: ' num2str(mean(nanrms(reshape(wf,size(wf,1)*size(wf,2),size(wf,3)))),'%0.4f') ' ' char(956) 'm']);

% Filter - Forward Moving
gain = sign((freq{2}/RunLog.samplerate(2)).*(freq{3}/RunLog.samplerate(3)))/2+0.5;
WF1 = WF.*gain.^2;
wf1 = real(ifftn(ifftshift(sqrt((WF1)*prod(RunLog.samplerate)*numel(WF1)).*exp(1i*phase)))).*mask;
disp(['Forward OPDrms: ' num2str(mean(nanrms(reshape(wf1,size(wf1,1)*size(wf1,2),size(wf1,3)))),'%0.4f') ' ' char(956) 'm']);

% Filter - Velocity
dist = abs(vFilter(1)*((freq{2}/RunLog.samplerate(2))-vFilter(5))+vFilter(2)*(freq{1}/RunLog.samplerate(1))-(freq{3}/RunLog.samplerate(3)))/sqrt((vFilter(1))^2+(vFilter(2))^2+1);
gain = sqrt(1./(1+(dist/vFilter(3)).^(+2*vFilter(4))));
WF2 = WF.*gain.^2;
wf2 = real(ifftn(ifftshift(sqrt((WF2)*prod(RunLog.samplerate)*numel(WF2)).*exp(1i*phase)))).*mask;
disp(['Velocity OPDrms: ' num2str(mean(nanrms(reshape(wf2,size(wf2,1)*size(wf2,2),size(wf2,3)))),'%0.4f') ' ' char(956) 'm']);

% Filter - Baseline
WF3 = zeros(size(WF));
for aa=1:size(WF,1)
    for bb=1:size(WF,2)
        WF3(aa,bb,:) = ipermute(baseline(permute(log10(WF(aa,bb,:)),[3 2 1])),[3 2 1]);
    end
end
WF3 = 10.^WF3;
wf3 = real(ifftn(ifftshift(sqrt((WF3)*prod(RunLog.samplerate)*numel(WF3)).*exp(1i*phase)))).*mask;
disp(['Baseline OPDrms: ' num2str(mean(nanrms(reshape(wf3,size(wf3,1)*size(wf3,2),size(wf3,3)))),'%0.4f') ' ' char(956) 'm']);

% Filter - FVB
dist = abs(vFilter(1)*((freq{2}/RunLog.samplerate(2))-vFilter(5))+vFilter(2)*(freq{1}/RunLog.samplerate(1))-(freq{3}/RunLog.samplerate(3)))/sqrt((vFilter(1))^2+(vFilter(2))^2+1);
gain = sqrt(1./(1+(dist/vFilter(3)).^(+2*vFilter(4))));
WF4 = WF3.*gain.^2;
gain = sign((freq{2}/RunLog.samplerate(2)).*(freq{3}/RunLog.samplerate(3)))/2+0.5;
WF4 = WF4.*gain.^2;
wf4 = real(ifftn(ifftshift(sqrt((WF4)*prod(RunLog.samplerate)*numel(WF4)).*exp(1i*phase)))).*mask;
disp(['FVB OPDrms: ' num2str(mean(nanrms(reshape(wf4,size(wf4,1)*size(wf4,2),size(wf4,3)))),'%0.4f') ' ' char(956) 'm']);

% Filter - VB
dist = abs(vFilter(1)*((freq{2}/RunLog.samplerate(2))-vFilter(5))+vFilter(2)*(freq{1}/RunLog.samplerate(1))-(freq{3}/RunLog.samplerate(3)))/sqrt((vFilter(1))^2+(vFilter(2))^2+1);
gain = sqrt(1./(1+(dist/vFilter(3)).^(+2*vFilter(4))));
WF5 = WF3.*gain.^2;
wf5 = real(ifftn(ifftshift(sqrt((WF5)*prod(RunLog.samplerate)*numel(WF5)).*exp(1i*phase)))).*mask;
disp(['VB OPDrms: ' num2str(mean(nanrms(reshape(wf5,size(wf5,1)*size(wf5,2),size(wf5,3)))),'%0.4f') ' ' char(956) 'm']);


%%
close(findobj('type','figure','number',1));
f1 = figure(1);
log_range = -14;
clim = [-14 -8];
aspect = [1 1 45];
scolor = parula(2);

subplot(2,2,1);
patch(isocaps(freq{1}(end/2+1:end),freq{2}(end/2+1:end),freq{3}(end/2+1:end),log10(WF1(end/2+1:end,end/2+1:end,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq{1}(end/2+1:end),freq{2}(end/2+1:end),freq{3}(end/2+1:end),log10(WF1(end/2+1:end,end/2+1:end,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
patch(isocaps(freq{1}(end/2+1:end),freq{2}(1:end/2+1),freq{3}(end/2+1:end),log10(WF1(1:end/2+1,end/2+1:end,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none','facealpha',faceAlpha);
patch(isosurface(freq{1}(end/2+1:end),freq{2}(1:end/2+1),freq{3}(end/2+1:end),log10(WF1(1:end/2+1,end/2+1:end,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud','facealpha',faceAlpha);%,'specularstrength',0.375);
grid on;
daspect(aspect);
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[-1 1]);
zlim(2e4*[0 1]);
caxis(clim);
xlabel('$\xi_x$','Interpreter','Latex');
ylabel('$\xi_y$','Interpreter','Latex');
zlabel('$f$','Interpreter','Latex');
title('Backward Moving Filter','Interpreter','Latex');
view(-45,15);
material dull;
camlight;

subplot(2,2,2);
patch(isocaps(freq{1},freq{2}(end/2+1:end),freq{3}(end/2+1:end),log10(WF2(end/2+1:end,:,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq{1},freq{2}(end/2+1:end),freq{3}(end/2+1:end),log10(WF2(end/2+1:end,:,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
patch(isocaps(freq{1},freq{2}(1:end/2+1),freq{3}(end/2+1:end),log10(WF2(1:end/2+1,:,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none','facealpha',faceAlpha);
patch(isosurface(freq{1},freq{2}(1:end/2+1),freq{3}(end/2+1:end),log10(WF2(1:end/2+1,:,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud','facealpha',faceAlpha);%,'specularstrength',0.375);
% patch(isosurface(freq{1},freq{2},freq{3},log10(WF2),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
% dist = abs(vFilter(1)*((freq{2}/RunLog.samplerate(2))-vFilter(5))+vFilter(2)*(freq{1}/RunLog.samplerate(1))-(freq{3}/RunLog.samplerate(3)))/sqrt((vFilter(1))^2+(vFilter(2))^2+1);
% gain = sqrt(1./(1+(dist/vFilter(3)).^(+2*vFilter(4))));
% patch(isosurface(freq{1},freq{2},freq{3},gain.^2,0.99),'edgecolor','none','facecolor','red','facelighting','none','facealpha',faceAlpha);
grid on;
daspect(aspect);
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[-1 1]);
zlim(2e4*[0 1]);
caxis(clim);
xlabel('$\xi_x$','Interpreter','Latex');
ylabel('$\xi_y$','Interpreter','Latex');
zlabel('$f$','Interpreter','Latex');
title('Velocity Filter','Interpreter','Latex');
view(-45,15);
material dull;
camlight;

subplot(2,2,3);
patch(isocaps(freq{1},freq{2}(end/2+1:end),freq{3}(end/2+1:end),log10(WF3(end/2+1:end,:,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq{1},freq{2}(end/2+1:end),freq{3}(end/2+1:end),log10(WF3(end/2+1:end,:,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
patch(isocaps(freq{1},freq{2}(1:end/2+1),freq{3}(end/2+1:end),log10(WF3(1:end/2+1,:,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none','facealpha',faceAlpha);
patch(isosurface(freq{1},freq{2}(1:end/2+1),freq{3}(end/2+1:end),log10(WF3(1:end/2+1,:,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud','facealpha',faceAlpha);%,'specularstrength',0.375);
grid on;
daspect(aspect);
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[-1 1]);
zlim(2e4*[0 1]);
caxis(clim);
xlabel('$\xi_x$','Interpreter','Latex');
ylabel('$\xi_y$','Interpreter','Latex');
zlabel('$f$','Interpreter','Latex');
title('Baseline Filter','Interpreter','Latex');
view(-45,15);
material dull;
camlight;

subplot(2,2,4);
patch(isocaps(freq{1}(end/2+1:end),freq{2}(end/2+1:end),freq{3}(end/2+1:end),log10(WF4(end/2+1:end,end/2+1:end,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq{1}(end/2+1:end),freq{2}(end/2+1:end),freq{3}(end/2+1:end),log10(WF4(end/2+1:end,end/2+1:end,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
patch(isocaps(freq{1}(end/2+1:end),freq{2}(1:end/2+1),freq{3}(end/2+1:end),log10(WF4(1:end/2+1,end/2+1:end,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none','facealpha',faceAlpha);
patch(isosurface(freq{1}(end/2+1:end),freq{2}(1:end/2+1),freq{3}(end/2+1:end),log10(WF4(1:end/2+1,end/2+1:end,end/2+1:end)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud','facealpha',faceAlpha);%,'specularstrength',0.375);
grid on;
daspect(aspect);
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[-1 1]);
zlim(2e4*[0 1]);
caxis(clim);
xlabel('$\xi_x$','Interpreter','Latex');
ylabel('$\xi_y$','Interpreter','Latex');
zlabel('$f$','Interpreter','Latex');
title('BM+V+B Filter','Interpreter','Latex');
view(-45,15);
material dull;
camlight;

colorbar('location','south');
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Children(3).TickLabelInterpreter = 'latex';
f1.Children(4).TickLabelInterpreter = 'latex';
f1.Children(5).TickLabelInterpreter = 'latex';
f1.Children(1).Label.String = '$S_{xx}$ ($\mu m^2/Hz/m^{-2}$)';
f1.Children(1).Label.Interpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 5.5 7];
f1.Children(1).Position(1) = f1.Children(3).Position(1);
f1.Children(1).Position(3) = f1.Children(4).Position(3)+f1.Children(4).Position(1)-f1.Children(1).Position(1);
f1.Children(1).Position(2) = f1.Children(1).Position(2)-0.1;
f1.Children(2).Position(2) = f1.Children(2).Position(2)+0.06;
f1.Children(3).Position(2) = f1.Children(3).Position(2)+0.06;
for aa=1:length(f1.Children(1).Ticks)
    f1.Children(1).TickLabels{aa} = ['$10^{' num2str(f1.Children(1).Ticks(aa)) '}$'];
end

% saveas(f1,'dispersion_filters.eps','epsc');

%%
close(findobj('type','figure','number',2));
f2 = figure(2);
colororder(linspecer(5));

[sxx,frequency] = computeSXX(nanrms(reshape(wf,size(wf,1)*size(wf,2),size(wf,3))),'positiveonly',1,'samplerate',RunLog.samplerate(3));
[sxx1] = computeSXX(nanrms(reshape(wf1,size(wf1,1)*size(wf1,2),size(wf1,3))),'positiveonly',1,'samplerate',RunLog.samplerate(3));
[sxx2] = computeSXX(nanrms(reshape(wf2,size(wf2,1)*size(wf2,2),size(wf2,3))),'positiveonly',1,'samplerate',RunLog.samplerate(3));
[sxx3] = computeSXX(nanrms(reshape(wf3,size(wf3,1)*size(wf3,2),size(wf3,3))),'positiveonly',1,'samplerate',RunLog.samplerate(3));
[sxx4] = computeSXX(nanrms(reshape(wf4,size(wf4,1)*size(wf4,2),size(wf4,3))),'positiveonly',1,'samplerate',RunLog.samplerate(3));
% [sxx5] = computeSXX(nanrms(reshape(wf5,size(wf5,1)*size(wf5,2),size(wf5,3))),'positiveonly',1,'samplerate',RunLog.samplerate(3));

loglog(frequency,sxx,frequency,sxx1,frequency,sxx2,frequency,sxx3,frequency,sxx4,'linewidth',1.25);
grid on;
xlabel('Frequency (Hz)','interpreter','latex');
ylabel('$S_{xx}\ (\mu m^2/Hz)$','interpreter','latex');
% title('$OPD_{RMS}(t)$','interpreter','latex');
legend('No Filter','Backward Moving Filter','Velocity Filter','Baseline Filter','BM+V+B Filter','interpreter','latex');
f2.Children(2).TickLabelInterpreter = 'latex';
f2.Units = 'inches';
f2.Position = [1 1 5.5 3.5];

% saveas(f2,'dispersion_filters_sxx.eps','epsc');

%% 
close(findobj('type','figure','number',3));
f3 = figure(3);
fRange = [0 3e3];
xiRange = 40;

subplot(3,2,1);
surf(squeeze(freq{3}),freq{2},log10(squeeze(WF(end/2+1,:,:))),'linestyle','none');
view(2);
grid on;
xlabel('$f\ (Hz)$','interpreter','latex');
ylabel('$\xi_x\ (m^{-1})$','interpreter','latex');
title('No Filter','interpreter','latex');
xlim(fRange);
ylim(xiRange*[-1 1]);
caxis(clim);

subplot(3,2,2);
surf(squeeze(freq{3}),freq{2},log10(squeeze(WF1(end/2+1,:,:))),'linestyle','none');
view(2);
grid on;
xlabel('$f\ (Hz)$','interpreter','latex');
ylabel('$\xi_x\ (m^{-1})$','interpreter','latex');
title('Backward Moving Filter','interpreter','latex');
xlim(fRange);
ylim(xiRange*[-1 1]);
caxis(clim);

subplot(3,2,3);
surf(squeeze(freq{3}),freq{2},log10(squeeze(WF2(end/2+1,:,:))),'linestyle','none');
view(2);
grid on;
xlabel('$f\ (Hz)$','interpreter','latex');
ylabel('$\xi_x\ (m^{-1})$','interpreter','latex');
title('Velocity Filter','interpreter','latex');
xlim(fRange);
ylim(xiRange*[-1 1]);
caxis(clim);

subplot(3,2,4);
surf(squeeze(freq{3}),freq{2},log10(squeeze(WF3(end/2+1,:,:))),'linestyle','none');
view(2);
grid on;
xlabel('$f\ (Hz)$','interpreter','latex');
ylabel('$\xi_x\ (m^{-1})$','interpreter','latex');
title('Baseline Filter','interpreter','latex');
xlim(fRange);
ylim(xiRange*[-1 1]);
caxis(clim);

subplot(3,2,5);
surf(squeeze(freq{3}),freq{2},log10(squeeze(WF5(end/2+1,:,:))),'linestyle','none');
view(2);
grid on;
xlabel('$f\ (Hz)$','interpreter','latex');
ylabel('$\xi_x\ (m^{-1})$','interpreter','latex');
title('V+B Filter','interpreter','latex');
xlim(fRange);
ylim(xiRange*[-1 1]);
caxis(clim);

subplot(3,2,6);
surf(squeeze(freq{3}),freq{2},log10(squeeze(WF4(end/2+1,:,:))),'linestyle','none');
view(2);
grid on;
xlabel('$f\ (Hz)$','interpreter','latex');
ylabel('$\xi_x\ (m^{-1})$','interpreter','latex');
title('BM+V+B Filter','interpreter','latex');
xlim(fRange);
ylim(xiRange*[-1 1]);
caxis(clim);

colorbar('location','east');
f3.Units = 'inches';
f3.Position = [1 1 5.5 4.5];
f3.Children(1).TickLabelInterpreter = 'latex';
spPosition = zeros(7,4);
spPosition(1,:) = f3.Children(1).Position;
for aa=2:7
    f3.Children(aa).TickLabelInterpreter = 'latex';
    f3.Children(aa).Layer = 'top';
	spPosition(aa,:) = f3.Children(aa).Position;
end
f3.Children(1).Label.String = '$S_{xx}$ ($\mu m^2/Hz/m^{-2}$)';
f3.Children(1).Label.Interpreter = 'latex';
spPosition([2 4 6],[1 3]) = spPosition([2 4 6],[1 3])-0.05;
spPosition([3 5 7],3) = spPosition([3 5 7],3)-0.05;
spPosition(2:7,2) = spPosition(2:7,2)+0.015;
spPosition(2:7,4) = spPosition(2:7,4)-0.03;
spPosition(1,2) = spPosition(2,2);
spPosition(1,4) = spPosition(6,4)+spPosition(6,2)-spPosition(1,2);
spPosition(1,1) = spPosition(1,1)+0.05;
for aa=1:7
    f3.Children(aa).Position = spPosition(aa,:);
end
for aa=1:length(f3.Children(1).TickLabels)
    f3.Children(1).TickLabels{aa} = ['$10^{' f3.Children(1).TickLabels{aa} '}$'];
end

saveas(f3,'dispersion_filters_slices.eps','epsc');




%%
% fileID = fopen('dispersion_filters.txt','w');
% fprintf(fileID,'\\begin{tabular}{c c}\n');
% fprintf(fileID,'  Filter & $\\opdrms\\ (\\mu m)$ \\\\ \n');
% fprintf(fileID,'  \\hline \\hline \n');
% fprintf(fileID,['  None & ' num2str(mean(nanrms(reshape(wf,size(wf,1)*size(wf,2),size(wf,3)))),'%0.4f') ' \\\\ \n']);
% fprintf(fileID,['  Forward & ' num2str(mean(nanrms(reshape(wf1,size(wf,1)*size(wf,2),size(wf,3)))),'%0.4f') ' \\\\ \n']);
% fprintf(fileID,['  Velocity & ' num2str(mean(nanrms(reshape(wf2,size(wf,1)*size(wf,2),size(wf,3)))),'%0.4f') ' \\\\ \n']);
% fprintf(fileID,['  Baseline & ' num2str(mean(nanrms(reshape(wf3,size(wf,1)*size(wf,2),size(wf,3)))),'%0.4f') ' \\\\ \n']);
% fprintf(fileID,['  V+B & ' num2str(mean(nanrms(reshape(wf5,size(wf,1)*size(wf,2),size(wf,3)))),'%0.4f') ' \\\\ \n']);
% fprintf(fileID,['  F+V+B & ' num2str(mean(nanrms(reshape(wf4,size(wf,1)*size(wf,2),size(wf,3)))),'%0.4f') ' \\\\ \n']);
% fprintf(fileID,'\\end{tabular}\n');


