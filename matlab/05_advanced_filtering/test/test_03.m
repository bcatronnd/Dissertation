close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901002';
[wf,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);
load([directory testPoint '_CDAQ.mat'],'scanData','DAQcopy');


sensors = 1:16;
nPtns = 2^12;

% [b,a] = butter(1,5000/DAQcopy.Rate*2);
% scanData = filter(b,a,scanData,[],1);


% 01 - Compute SXX
WF = computeSXX(wf,'dim',3,'samplerate',DAQcopy.Rate,'blocksize',nPtns);
[Y,freq] = computeSXX(scanData(:,sensors)','dim',2,'samplerate',DAQcopy.Rate,'blocksize',nPtns);

f1 = figure(1);
subplot(3,1,1);
loglog(freq(end/2+1:end),squeeze(mean(WF(:,:,end/2+1:end),[1 2],'omitnan')));
grid on;
title('Wavefront','interpreter','latex');
ylabel('$S_{XX}\ (\mu m^2/Hz)$','interpreter','latex');
subplot(3,1,2);
loglog(freq(end/2+1:end),Y(1:6,end/2+1:end));
grid on;
title('Microphones','interpreter','latex');
ylabel('$S_{XX}\ (V^2/Hz)$','interpreter','latex');
legend('Mic 1','Mic 2','Mic 3','Mic 4','Mic 5','Mic 6','interpreter','latex');
subplot(3,1,3);
loglog(freq(end/2+1:end),Y(7:16,end/2+1:end));
grid on;
title('Accelerometers','interpreter','latex');
ylabel('$S_{XX}\ (V^2/Hz)$','interpreter','latex');
xlabel('Frequency (Hz)','interpreter','latex');
legend('Accel 1','Accel 2','Accel 3','Accel 4','Accel 5','Accel 6','Accel 7','Accel 8','Accel 9','Accel 10','interpreter','latex');
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Children(4).TickLabelInterpreter = 'latex';
f1.Children(5).TickLabelInterpreter = 'latex';

% 02 - Spectral POD of WF
M = size(WF,1,2);
N = size(WF,3);
WF = reshape(WF,prod(M),N);
DataLocation = ~isnan(WF);
WF = reshape(WF(DataLocation),[],N);
[Phi,~] = eig(1/N*(WF*WF'),'vector');
Phi = flip(Phi,2);
a = (WF'*Phi);

% f2 = figure(2);
% subplot(3,1,1);
% loglog(freq(end/2+1:end),abs(a(end/2+1:end,1:5))');
% grid on;
% title('Wavefront - POD Mode Coefficients','interpreter','latex');
% ylabel('$S_{XX}\ (\mu m^2/Hz)$','interpreter','latex');
% subplot(3,1,2);
% loglog(freq(end/2+1:end),Y(1:6,end/2+1:end));
% grid on;
% title('Microphones','interpreter','latex');
% ylabel('$S_{XX}\ (V^2/Hz)$','interpreter','latex');
% subplot(3,1,3);
% loglog(freq(end/2+1:end),Y(7:16,end/2+1:end));
% grid on;
% title('Accelerometers','interpreter','latex');
% ylabel('$S_{XX}\ (V^2/Hz)$','interpreter','latex');
% xlabel('Frequency (Hz)','interpreter','latex');
% f2.Children(1).TickLabelInterpreter = 'latex';
% f2.Children(2).TickLabelInterpreter = 'latex';
% f2.Children(3).TickLabelInterpreter = 'latex';

% 03 - LSE
L = (a'*Y')/(Y*Y');
a_lse = (L*Y)';
a_ao = a-a_lse;

% f3 = figure(3);
% subplot(3,1,1);
% loglog(freq(end/2+1:end),abs(a(end/2+1:end,1))',freq(end/2+1:end),abs(a_lse(end/2+1:end,1))',freq(end/2+1:end),abs(a_ao(end/2+1:end,1))');
% grid on;
% legend('original','lse','ao');


% % Step 1 - FFT
% 
% WF = fft(wf,nPtns,3);
% Y = fft(y(sensors,:),nPtns,2);
% freq = (-0.5:1/nPtns:0.5-1/nPtns)*sampleRate;
% 
% % Step 2 - Compute Spectral POD
% M = size(WF,1,2);
% N = size(WF,3);
% WF = reshape(WF,prod(M),N);
% DataLocation = ~isnan(WF);
% WF = reshape(WF(DataLocation),[],N);
% [Phi,~] = eig(1/N*(WF*WF'),'vector');
% Phi = flip(Phi,2);
% a = (WF'*Phi);
% 
% % Step 3 - LSE
% L = (a'*conj(Y'))/(Y*conj(Y'));
% 
% a_lse = (L*Y)';
% a_ao = a-a_lse;
% 
% 
% 
% 
% figure(1);
% loglog(freq,(abs(fftshift(a(:,2)))).^2/nPtns/sampleRate,freq,(abs(fftshift(a_lse(:,2)))).^2/nPtns/sampleRate,freq,(abs(fftshift(a_ao(:,2)))).^2/nPtns/sampleRate);
% grid on;
% legend('original','lse','ao');
% 

