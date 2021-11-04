close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);
BlockSize = 2^12;
WF = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));


BlockSize = [size(WF,[1 2]) BlockSize];
Frequency.y = reshape(-1/2:1/BlockSize(1):1/2-1/BlockSize(1),[],1);
Frequency.x = reshape(-1/2:1/BlockSize(2):1/2-1/BlockSize(2),1,[]);
Frequency.t = reshape(-1/2:1/BlockSize(3):1/2-1/BlockSize(3),1,1,[]);
Frequency.rho = sqrt(Frequency.x.^2+Frequency.y.^2);

width = 0.05;
order = 2;
u = 152;
v = 3;

U = u*RunLog.samplerate(2)/RunLog.samplerate(3);
V = v*RunLog.samplerate(1)/RunLog.samplerate(3);
dist = abs(U*Frequency.x+V*Frequency.y-Frequency.t)/sqrt(U^2+V^2+1);
gain = sqrt(1./(1+(dist/width).^(+2*order))).*sqrt(1./(1+(Frequency.rho/0.1).^(-2*2)));

[Phi,a,M,N,DataLocation] = computePOD(WF,3);
[WF,freq] = simpleDispersion(WF,'BlockSize',BlockSize(3),'SampleRate',RunLog.samplerate);

blMode = zeros(1,size(a,2));
blCut = 0.9;
tic;
for aa=1:size(a,2)
    wf = simpleDispersion(inversePOD(Phi(:,aa),a(:,aa),M,N,DataLocation,3),'BlockSize',BlockSize(3));
    if sum(wf.*gain.^2,'all')/sum(wf,'all')>=blCut
        blMode(aa) = aa;
    end
    clear wf;
end
toc;
blMode = (randi(2,1,size(a,2))-1).*(1:size(a,2));

blMode(blMode==0) = [];
disp(['Number of Modes for BL: ' num2str(length(blMode))]);
[wf] = simpleDispersion(inversePOD(Phi(:,blMode),a(:,blMode),M,N,DataLocation,3),'BlockSize',BlockSize(3),'SampleRate',RunLog.samplerate);

log_range = -14;
f1 = figure(1);
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
view(3);
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
view(3);
material dull;
camlight;
f1.Children(1).TickLabelInterpreter = 'latex';


%%
% close all;
% 
% 
% log_range = -14;

% 
% f1 = figure(1);
%     subplot(2,length(PODmodes),aa);
%     scolor = parula(2);
%     patch(isocaps(freq{1},freq{2},freq{3}(end/2+1:end),wf(:,:,end/2+1:end),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
%     patch(isosurface(freq{1},freq{2},freq{3},wf,log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
%     grid on;
%     daspect([1 1 50]);
%     xlim(RunLog.samplerate(1)/2*[-1 1]);
%     ylim(RunLog.samplerate(2)/2*[-1 1]);
%     zlim(RunLog.samplerate(3)/2*[0 1]);
%     xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
%     ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
%     zlabel('$f\ (Hz)$','Interpreter','Latex');
%     title({['Modes: ' num2str(min(PODmodes{aa})) '-' num2str(max(PODmodes{aa}))];' '},'Interpreter','Latex');
%     view(0,0);
%     material dull;
%     camlight;
%     f1.Children(1).TickLabelInterpreter = 'latex';
%     
%     
%     subplot(2,length(PODmodes),aa+length(PODmodes));
%     scolor = parula(2);
%     patch(isocaps(freq{1},freq{2},freq{3}(end/2+1:end),wf(:,:,end/2+1:end),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
%     patch(isosurface(freq{1},freq{2},freq{3},wf,log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
%     grid on;
%     daspect([1 1 50]);
%     xlim(RunLog.samplerate(1)/2*[-1 1]);
%     ylim(RunLog.samplerate(2)/2*[-1 1]);
%     zlim(RunLog.samplerate(3)/2*[0 1]);
%     xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
%     ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
%     zlabel('$f\ (Hz)$','Interpreter','Latex');
%     view(90,0);
%     material dull;
%     camlight;
%     f1.Children(1).TickLabelInterpreter = 'latex';
% end
% f1.Units = 'inches';
% f1.Position = [1 1 7.5 5.5];
% saveas(f1,'filter_pod.eps','epsc');


