close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/ProcessedWavefronts/';
testPoint = {'20210417004' '20210417007' '20210417008'};
plotFrame = [1000 1014 1005];
plotFreq = [9000 14000 18000];
c0 = 340;

load([directory '20210417009' '_RL.mat']);
MicHeight = str2double(extract(RunLog.notes{3},digitsPattern))/1e3;
BeamHeight = str2double(extract(RunLog.notes{4},digitsPattern))/1e3;
R = mean(BeamHeight);
Ap = diff(BeamHeight);
% Lambda = Ap*[1 0.75 0.5];
Lambda = c0./plotFreq;

% Measured
f1 = figure(1);
colormap(redblue);
for aa=1:length(testPoint)
    [WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint{aa} '_WF.mat'],'Scale',1e6);%,'ZernikeRemoval',1:3);
    %     load([directory testPoint{aa} '_CDAQ.mat']);
    %     load('spherical_sample_win.mat');
    runStats = str2double(extract(RunLog.notes{2},digitsPattern));
    WF = WFfilter(WF,'time-bandpass',(runStats(1)+10*[-1 1])/RunLog.samplerate(3));
    subplot(length(testPoint),1,aa);
    surf(WF(:,:,plotFrame(aa)),'linestyle','none');
    view(2);
    axis image off tight;
    caxis(max(abs(WF(:,:,plotFrame(aa))),[],'all')*[-1 1]);
    title([num2str(runStats(1)) ' Hz'],'interpreter','latex');
end
sgtitle('Measured Wavefronts','interpreter','latex');
f1.Units = 'inches';
f1.Position = [1 1 2 5];

%%
plotPhase = [1.75 5.1 -0.2];
% Simulated
nLenslets = 32;
zMax = 5;
dz = 0.005;

kgd = 2.27e-4;
[x,y,z] = meshgrid(-Ap/2*(1-1/nLenslets):Ap/nLenslets:Ap/2*(1-1/nLenslets),-Ap/2*(1-1/nLenslets):Ap/nLenslets:Ap/2*(1-1/nLenslets),-zMax:dz:zMax);
r = sqrt(x(:,:,1).^2+y(:,:,1).^2);
t = atan2(y(:,:,1),x(:,:,1));
mask = ones(nLenslets);
mask(r>0.5*Ap) = NaN;
R2 = sqrt(x.^2+(y+R).^2+z.^2);
win = reshape(tukeywin(size(z,3)),1,1,[]);
k = 2*pi./Lambda;
f2 = figure(2);
colormap(redblue);
for aa=1:length(testPoint)
    
    OPD = mask.*(kgd/c0^2*sum(win./R2.*exp(1i*k(aa)*R2),3)*dz*exp(-1i*plotPhase(aa)));
    OPD = ZernikeRemoval(1:3,OPD,r.*mask,t.*mask,'noll');
    subplot(length(testPoint),1,aa);
    surf(real(OPD),'linestyle','none');
    view(2);
    axis image off tight;
    caxis(max(abs(real(OPD)),[],'all')*[-1 1]);
    title([num2str(plotFreq(aa)) ' Hz'],'interpreter','latex');
end
sgtitle('Simulated Wavefronts','interpreter','latex');
f2.Units = 'inches';
f2.Position = [1 1 2 5];


saveas(f1,'spherical_plot_measured.eps','epsc');
saveas(f2,'spherical_plot_simulated.eps','epsc');





