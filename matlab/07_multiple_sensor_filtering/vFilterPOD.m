close all; clc; clearvars;

u = 152;
v = 3;
d_cut = 0.06;
n = 50;


directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);
load([directory testPoint '_CDAQ.mat'],'scanData');


% WF = WFfilter(WF,'velocity-lowpass',[u*RunLog.samplerate(2)/RunLog.samplerate(3),v*RunLog.samplerate(1)/RunLog.samplerate(3),0.02,20]);
% Dispersion
BlockSize = 2^12;
WF = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));
[WF,freq] = simpleDispersion(WF,'BlockSize',BlockSize,'SampleRate',RunLog.samplerate,'Output','lin');

% Velocity Filter
BlockSize = size(WF);
Frequency.y = reshape(-1/2:1/BlockSize(1):1/2-1/BlockSize(1),[],1);
Frequency.x = reshape(-1/2:1/BlockSize(2):1/2-1/BlockSize(2),1,[]);
Frequency.t = reshape(-1/2:1/BlockSize(3):1/2-1/BlockSize(3),1,1,[]);
dist = abs(u*RunLog.samplerate(2)/RunLog.samplerate(3)*Frequency.x+v*RunLog.samplerate(1)/RunLog.samplerate(3)*Frequency.y-Frequency.t)/sqrt((u*RunLog.samplerate(2)/RunLog.samplerate(3))^2+(v*RunLog.samplerate(1)/RunLog.samplerate(3))^2+1);
gain = sqrt(1./(1+(dist/d_cut).^(+2*n)));
WF = log10(gain.^2.*WF);
% POD
[Phi,a,M,N,DataLocation] = computePOD(WF-mean(WF,'all'),3);
%%
podMode = 1:32;



WFpod = inversePOD(Phi(:,podMode),a(:,podMode),M,N,DataLocation,3)+mean(WF,'all');








%%

close all;
%%%%% Plot
log_range = -14;
f1 = figure(1);
subplot(1,2,1);
scolor = parula(2);
patch(isocaps(freq{1},freq{2},freq{3}(end/2+1:end),(WF(:,:,end/2+1:end)),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq{1},freq{2},freq{3},(WF),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
grid on;
daspect([1 1 50]);
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[-1 1]);
zlim(RunLog.samplerate(3)/2*[0 1]);
xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
zlabel('$f\ (Hz)$','Interpreter','Latex');
material dull;
camlight;
f1.Children(1).TickLabelInterpreter = 'latex';
view(3);

subplot(1,2,2);
scolor = parula(2);
patch(isocaps(freq{1},freq{2},freq{3}(end/2+1:end),real((WFpod(:,:,end/2+1:end))),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq{1},freq{2},freq{3},real((WFpod)),log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
grid on;
daspect([1 1 50]);
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[-1 1]);
zlim(RunLog.samplerate(3)/2*[0 1]);
xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
zlabel('$f\ (Hz)$','Interpreter','Latex');
material dull;
camlight;
f1.Children(1).TickLabelInterpreter = 'latex';
view(3);




