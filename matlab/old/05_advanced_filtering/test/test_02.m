close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[wf,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);
load([directory testPoint '_CDAQ.mat'],'scanData','DAQcopy');

[SXX_WF] = computeSXX(wf,'blockSize',2^12,'sampleRate',DAQcopy.Rate,'positiveOnly',1,'average',[1 2]);
[SXX_DAQ,freq] = computeSXX(scanData','blockSize',2^12,'sampleRate',DAQcopy.Rate,'positiveOnly',1);

f1 = figure(1);
subplot(3,1,1);
loglog(freq,reshape(SXX_WF,[],2^11));
grid on;
xlim(10.^[2 log10(30000)]);
title('Wavefront','interpreter','latex');
ylabel('$S_{XX}$ ($\mu$m$^2$/Hz)','interpreter','latex');
subplot(3,1,2);
loglog(freq,SXX_DAQ(1:6,:));
grid on;
xlim(10.^[2 log10(30000)]);
title('Microphones','interpreter','latex');
ylabel('$S_{XX}$ (V$^2$/Hz)','interpreter','latex');
subplot(3,1,3);
loglog(freq,SXX_DAQ(7:16,:));
grid on;
xlim(10.^[2 log10(30000)]);
title('Accelerometers','interpreter','latex');
xlabel('Frequency (Hz)','interpreter','latex');
ylabel('$S_{XX}$ (V$^2$/Hz)','interpreter','latex');

f1.Children(1).TickLabelInterpreter = 'latex';
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Children(3).TickLabelInterpreter = 'latex';
