close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);
load([directory testPoint '_CDAQ.mat'],'scanData');

BlockSizeTF = 2^8;
BlockSize = 2^12;

tf = estimateTF(computeSXX(scanData,'blocksize',BlockSizeTF,'dim',1,'window',@hann),squeeze(computeSXX(WF,'blocksize',BlockSizeTF,'dim',3,'window',@hann,'average',[1 2])),BlockSize);

[WF,freq] = computeSXX(WF,'blocksize',BlockSize,'dim',3,'window',@hann,'average',[1 2],'samplerate',RunLog.samplerate(3),'positiveonly',1);
scanData = computeSXX(scanData,'blocksize',BlockSize,'dim',1,'window',@hann,'positiveonly',1);




% figure(1);
% loglog(tf(end/2+1:end,:));
% grid on;
% title('Transfer Functions');

figure(2);
loglog(squeeze(freq),squeeze(WF),'k-',squeeze(freq),scanData.*tf(end/2+1:end,:));
grid on;







