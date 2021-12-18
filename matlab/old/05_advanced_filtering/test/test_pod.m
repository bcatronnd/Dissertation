close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[wf,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);
load([directory testPoint '_CDAQ.mat'],'scanData','timeStamp');
sampleRate = mean(1./diff(timeStamp));

range = [1000 3e4];
[wf_ao] = filterLSE_POD(wf,scanData(1:size(wf,3),1:2)',range/sampleRate);
% wf_ao = wf;






BlockSize = 2^10;
wf_ao = wf_ao(:,:,1:floor(size(wf_ao,3)/BlockSize)*BlockSize);
wf_ao = cat(1,wf_ao,NaN*zeros(2^nextpow2(size(wf_ao,1))-size(wf_ao,1),size(wf_ao,2),size(wf_ao,3)));
wf_ao = cat(2,wf_ao,NaN*zeros(size(wf_ao,1),2^nextpow2(size(wf_ao,2))-size(wf_ao,2),size(wf_ao,3)));
[wf_ao,freq] = simpleDispersion(wf_ao,'BlockSize',BlockSize,'SampleRate',RunLog.samplerate);


log_range = -15;
disp(['Total Energy Inside Surface: ' num2str(sum(10.^wf_ao(wf_ao>log_range),'all')/sum(10.^wf_ao,'all')*100,'%0.1f') '%']);

%%%%% Plot
f1 = figure(1);
scolor = parula(2);
patch(isocaps(freq{1},freq{2},freq{3}(end/2+1:end),wf_ao(:,:,end/2+1:end),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq{1},freq{2},freq{3},wf_ao,log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
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
