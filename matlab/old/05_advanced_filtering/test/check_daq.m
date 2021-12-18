close all; clc; clearvars;

nPtns = 2^12;
nChannels = 16;
directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';


load([directory 'RunLog.mat']);

SXX = NaN*zeros(nChannels,nPtns/2,length(RunLog.TestPoint));
freq = NaN*zeros(length(RunLog.TestPoint),nPtns/2);
sampleRate = NaN*zeros(length(RunLog.TestPoint),2);

for aa=1:length(RunLog.TestPoint)
    if RunLog.DaqData(aa)
        load([directory num2str(RunLog.TestPoint(aa)) '_CDAQ.mat'],'scanData','timeStamp');
        sampleRate(aa,:) = [mean(1./diff(timeStamp)) std(1./diff(timeStamp))];
        [SXX(:,:,aa),freq(aa,:)] = computeSXX(scanData','dim',2,'samplerate',sampleRate(aa,1),'blocksize',nPtns,'positiveonly',1);
        clear scanData DAQcopy;
    end
end
clear aa;

[~,iBPF] = min(abs(freq-RunLog.BladePassingFrequency'),[],2);
RossiterMode = 1;
RossiterFrequency = RossiterMode*RunLog.VelocityFreeStream./(RunLog.Length*0.0254)./(RunLog.MachCalc+1+coth(RossiterMode*pi*RunLog.Depth./RunLog.Length));
[~,iRF] = min(abs(freq-RossiterFrequency'),[],2);

%%
close all;
for aa=1:nChannels
    figure(aa);
    subplot(2,1,1)
    loglog(freq',squeeze(SXX(aa,:,:)),RunLog.BladePassingFrequency,diag(squeeze(SXX(aa,iBPF,:))),'ko',RossiterFrequency,diag(squeeze(SXX(aa,iRF,:))),'k*');
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Sxx (V^2/Hz)');
    subplot(2,1,2);
    loglog(freq'./RunLog.VelocityFreeStream,squeeze(SXX(aa,:,:)),RunLog.BladePassingFrequency./RunLog.VelocityFreeStream,diag(squeeze(SXX(aa,iBPF,:))),'ko',RossiterFrequency./RunLog.VelocityFreeStream,diag(squeeze(SXX(aa,iRF,:))),'k*');
    grid on;
    xlabel('St/l (1/m)');
    ylabel('Sxx (V^2/Hz)');
    
    sgtitle(['Channel ' num2str(aa)]);
    
    
    
    
    
end




% testPoints = num2str([20210831001:20210831007 20210831009:20210831020 20210901001:20210901023]');
% M = [0	0.3	0.4	0.5	0.3	0.4	0.5	0.3	0.4	0.5	0.3	0.4	0.5	0.3	0.4	0.5	0.3	0.4	0.5	0.3	0.4	0.5	0	0.3	0.4	0.5	0.3	0.4	0.5	0.3	0.4	0.5	0.3	0.4	0.5	0.3	0.4	0.5	0.3	0.4	0.5	0]';
% T = [23.9	26.5	27.4	30.1	30.1	31	25.4	25.6	26.5	25.3	25.8	27.1	25.9	26.4	27.5	28.1	28.5	29.8	29.2	29.5	30.6	24.3	24.9	26.1	26.4	27	28	27.9	28.4	29.6	29.5	29.9	30.8	28.5	28.7	29.8	29.7	30	31	29.9	30	31;
%     25	27.4	30.2	30.1	31	32.9	25.6	26.5	28.4	25.8	27.1	28.8	26.4	27.5	29.5	28.5	29.8	31.6	29.5	30.6	32.4	24.9	26.1	28.2	27	28	29.8	28.4	29.6	31.4	29.9	30.8	32.4	28.7	29.8	31.7	30	31	32.7	30	31	32.6];
% T = mean(T,1)'+273;
% 
% 
% c = sqrt(1.4*287*T);
% U = c.*M;
% 
% 
% l = 0.125*0.0254;
% 
% nPtns = 2^12;
% nChannels = 16;
% 
% SXX = zeros(nChannels,nPtns/2,size(testPoints,1));
% freq = zeros(size(testPoints,1),nPtns/2);
% sampleRate = zeros(size(testPoints,1),2);
% 
% for aa=1:size(testPoints,1)
%     load([directory testPoints(aa,:) '_CDAQ.mat'],'scanData','timeStamp');
%     sampleRate(aa,:) = [mean(1./diff(timeStamp)) std(1./diff(timeStamp))];
%     [SXX(:,:,aa),freq(aa,:)] = computeSXX(scanData','dim',2,'samplerate',sampleRate(aa,1),'blocksize',nPtns,'positiveonly',1);
%     clear scanData DAQcopy;
% end
% clear aa;
% 
% St_l = freq./U;
% close all;
% for aa=1:nChannels
%     figure(aa);
%     subplot(2,1,1)
%     loglog(freq',squeeze(SXX(aa,:,:)));
%     grid on;
%     xlabel('Frequency (Hz)');
%     ylabel('Sxx (Units^2/Hz)');
%     subplot(2,1,2);
%     loglog(St_l',squeeze(SXX(aa,:,:)));
%     grid on;
%     xlabel('St/l (1/m)');
%     ylabel('Sxx (Units^2/Hz)');
%     
%     sgtitle(['Channel ' num2str(aa)]);
% end


