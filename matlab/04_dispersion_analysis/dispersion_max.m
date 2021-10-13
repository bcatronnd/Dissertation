close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[wf,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);

blockSize = 2^12;
sampleRate = [277.7778 49020];

wf = squeeze(wf(round(size(wf,1)/2),:,:));

[sxx1,freq1] = computeSXX(wf,'average',1,'blocksize',blockSize,'dim',2,'positiveonly',1,'samplerate',sampleRate(2),'window',@hann);

% 2-D Dispersion
wf = wf(:,1:floor(size(wf,2)/blockSize)*blockSize);
wf = reshape(wf,size(wf,1),blockSize,size(wf,2)/blockSize);
window = hann(size(wf,1)).*hann(size(wf,2))';
cw = 1/sqrt(sum(window.^2,'all')/blockSize);
sxx2 = mean(cw*fftshift(fftshift((abs(fft2(wf.*window))).^2,1),2)/prod(size(wf,[1 2]))/prod(sampleRate),3);
for aa=1:2
    freq2{aa} = (-0.5:1/size(wf,aa):0.5-1/size(wf,aa))*sampleRate(aa);
end


f1 = figure(1);
loglog(freq1,sxx1,freq2{2},sampleRate(1)*max(sxx2,[],1));
grid on;
xlim([1e2 2e4]);
xlabel('Temporal Frequency, $f$ ($Hz$)','interpreter','latex');
ylabel('$S_{xx}$ ($\mu m^2/Hz$)','interpreter','latex');
legend('Power Spectra - Time','max(2-D Dispersion)$\cdot f_{s,x}$','interpreter','latex');
f1.Units = 'inches';
f1.Position = [1 1 5.5 3.25];
f1.Children(2).TickLabelInterpreter = 'latex';

saveas(f1,'dispersion_max.eps','epsc');