close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901002';
[wf,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);
load([directory testPoint '_CDAQ.mat'],'scanData','timeStamp');

range = [0 2000];
sensors = 1:16;
nPtns = 2^12;

% 01 - Compute SXX
sampleRate = mean(1./diff(timeStamp));
WF = computeSXX(wf,'dim',3,'samplerate',sampleRate,'blocksize',nPtns);
[Y,freq] = computeSXX(scanData(:,sensors)','dim',2,'samplerate',sampleRate,'blocksize',nPtns);

% 02 - Spectral POD of WF
M = size(WF,1,2);
N = size(WF,3);
WF = reshape(WF,prod(M),N);
DataLocation = ~isnan(WF);
WF = reshape(WF(DataLocation),[],N);
[Phi,~] = eig(1/N*(WF*WF'),'vector');
Phi = flip(Phi,2);
a = (WF'*Phi);

% 03 - LSE Over Range
lPoints = and(abs(freq)<=range(2),abs(freq)>=range(1));
L = (a(lPoints',:)'*Y(:,lPoints)')/(Y(:,lPoints)*Y(:,lPoints)');
a_lse = (L*Y(:,lPoints))';
a_ao = a;
a_ao(lPoints',:) = a(lPoints',:)-a_lse;

f1 = figure(1);
loglog(freq,a(:,1),freq(lPoints),a_lse(:,1),freq,a_ao(:,1));
grid on;
xlim([0 inf]);
