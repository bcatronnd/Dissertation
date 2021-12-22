close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);

dbOffset = 1/20;
BlockSize = 2^10;

PlotColors = linspecer(2);
window = permute(hann(BlockSize),[3 2 1]);
cw = 1/sqrt(sum(window.^2,'all')/BlockSize);
freq = (-0.5:1/BlockSize:0.5-1/BlockSize)*RunLog.samplerate(3);

WFtemp = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
WFtemp = reshape(WFtemp,size(WFtemp,1),size(WFtemp,2),BlockSize,size(WFtemp,3)/BlockSize);
wf = log10(permute(mean(abs(cw*fftshift(fft(WFtemp.*window,[],3))).^2/BlockSize/RunLog.samplerate(3),[1 2 4],'omitnan'),[1 3 2]));

% wf = log10(wf');
wfb = baseline(wf');



%% 
close all;
figure(1);
semilogx(freq(end/2+1:end),wf(end/2+1:end),'linewidth',1.25,'color',PlotColors(1,:));
hold on;
semilogx(freq(end/2+1:end),dbOffset+wfb(end/2+1:end),'linewidth',1.25,'color',PlotColors(2,:));
grid on;



