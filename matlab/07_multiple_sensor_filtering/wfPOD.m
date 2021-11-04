close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);
load([directory testPoint '_CDAQ.mat'],'scanData');

BlockSize = 2^12;
% Wavefront Dispersion and POD
WF = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));
[WF,freq] = simpleDispersion(WF,'BlockSize',BlockSize,'SampleRate',RunLog.samplerate,'output','lin');
[Phi,a,M,N] = computePOD(WF,3);
aSign = sign(a);
% DAQ PSD
sxx = computeSXX(scanData,'blocksize',BlockSize,'dim',1,'window',@hann);
sxx = sxx./(1+(squeeze(freq{3})/2e3).^(2*5))./(1+(squeeze(freq{3})/1e2).^(-2*5));

%%
% LSE
channels = 1:16;
fRange = [1e2 1e4];
fPoints = ones(BlockSize,1);
fPoints(squeeze(abs(freq{3}))<fRange(1)) = 0;
fPoints(squeeze(abs(freq{3}))>fRange(2)) = 0;
fPoints = fPoints.*(1:BlockSize)';
fPoints(fPoints==0) = [];

% L = ((a(fPoints,:)')*sxx(fPoints,channels))/(sxx(fPoints,channels)'*sxx(fPoints,channels));
% L = a\sxx(:,channels);
L = (a'*sxx(:,channels))/(sxx(:,channels)'*sxx(:,channels));
aLSE = (L*sxx(:,channels)')';
aAO = aSign.*(abs(a)-abs(aLSE));


close all;

pMode = 1;
figure(1);
loglog(squeeze(freq{3}(end/2+1:end)),abs(a(end/2+1:end,pMode)),squeeze(freq{3}(end/2+1:end)),abs(aLSE(end/2+1:end,pMode)),squeeze(freq{3}(end/2+1:end)),abs(aAO(end/2+1:end,pMode)));
grid on;
legend('a','alse','aao');

