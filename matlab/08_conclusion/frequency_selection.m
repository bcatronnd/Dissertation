close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003'; % Flush Window
% testPoint = '20210901010'; % 12.75" Long Cavity
% testPoint = '20210901013'; % 9" Long Cavity

% Load Data
load([directory 'RunLog.mat']);
runlog = RunLog;
index = find(str2double(testPoint)==runlog.TestPoint);
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6,'zernikeremoval',1:3);
load([directory testPoint '_CDAQ.mat'],'scanData');

% Options
BlockSize = 2^10;
n = 1:4;




% Code
Uc = runlog.VelocityFreeStream(index)./(1+coth(n*pi*runlog.Depth(index)/runlog.Length(index)));
fRos = n./(runlog.Length(index)*0.0254/runlog.SpeedOfSound(index)+runlog.Length(index)*0.0254./Uc);


if isempty(RunLog)
    RunLog.samplerate = [2.777777777777778e+02,2.777777777777778e+02,49000];
end
[wf,wffreq] = computeSXX(WF,'samplerate',RunLog.samplerate(3),'blocksize',BlockSize,'average',[1 2]);
wf = permute(wf,[1 3 2]);
wffreq = permute(wffreq,[1 3 2]);
wfb = ipermute(10.^baseline(permute(log10(wf),[2 1])),[2 1]);
% WF = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
% WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
% WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));
% [WF,WFfreq] = simpleDispersion(WF,'BlockSize',BlockSize,'SampleRate',RunLog.samplerate,'output','lin');
[Y,Yfreq] = computeSXX(scanData','samplerate',RunLog.samplerate(3),'blocksize',BlockSize);
Ymic = mean(Y(3:6,:),1);
Ybmic = ipermute(10.^baseline(permute(log10(Ymic),[2 1])),[2 1]);
Yaccel = mean(Y(8:13,:),1);
Ybaccel = ipermute(10.^baseline(permute(log10(Yaccel),[2 1])),[2 1]);
% Yb = zeros(size(Y));
% for aa=1:size(Y,1)
%     Yb(aa,:) = ipermute(10.^baseline(permute(log10(Y(aa,:)),[2 1])),[2 1]);
% end

%%
close(findobj('type','figure','number',1));
f1 = figure(1);
colororder(linspecer(5));
sensors = [3 10];

subplot(4,1,1);
loglog(wffreq(end/2+1:end),wf(end/2+1:end),wffreq(end/2+1:end),wfb(end/2+1:end),'linewidth',1.25);
grid on;

subplot(4,1,2);
loglog(Yfreq(end/2+1:end),Ymic(end/2+1:end),Yfreq(end/2+1:end),Ybmic(end/2+1:end),'linewidth',1.25);
grid on;

subplot(4,1,3);
loglog(Yfreq(end/2+1:end),Yaccel(end/2+1:end),Yfreq(end/2+1:end),Ybaccel(end/2+1:end),'linewidth',1.25);grid on;

subplot(4,1,4);
semilogx(wffreq(end/2+1:end),log10(wf(end/2+1:end))-log10(wfb(end/2+1:end)),Yfreq(end/2+1:end),log10(Ymic(end/2+1:end))-log10(Ybmic(end/2+1:end)),Yfreq(end/2+1:end),log10(Yaccel(end/2+1:end))-log10(Ybaccel(end/2+1:end)),'linewidth',1.25);
grid on;



% dbOffset = 0;
% BlockSize = 2^10;
% 
% PlotColors = linspecer(2);
% window = permute(hann(BlockSize),[3 2 1]);
% cw = 1/sqrt(sum(window.^2,'all')/BlockSize);
% freq = (-0.5:1/BlockSize:0.5-1/BlockSize)*RunLog.samplerate(3);
% 
% WFtemp = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
% WFtemp = reshape(WFtemp,size(WFtemp,1),size(WFtemp,2),BlockSize,size(WFtemp,3)/BlockSize);
% wf = log10(permute(mean(abs(cw*fftshift(fft(WFtemp.*window,[],3))).^2/BlockSize/RunLog.samplerate(3),[1 2 4],'omitnan'),[1 3 2]));
% 
% % wf = log10(wf');
% wfb = baseline(wf');
% 
% 
% 
% %% 
% close all;
% figure(1);
% semilogx(freq(end/2+1:end),wf(end/2+1:end),'linewidth',1.25,'color',PlotColors(1,:));
% hold on;
% semilogx(freq(end/2+1:end),dbOffset+wfb(end/2+1:end),'linewidth',1.25,'color',PlotColors(2,:));
% grid on;



