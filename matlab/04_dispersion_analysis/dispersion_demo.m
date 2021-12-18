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
subplot(3,1,1);
loglog(freq1,sxx1);
grid on;
xlim([1e2 2e4]);
ylabel('$S_{xx}$ ($\mu m^2/Hz$)','interpreter','latex');
subplot(3,1,2);
surf(freq2{2},freq2{1},log10(sxx2),'linestyle','none','facecolor','interp');
view(2);
xlim([1e2 2e4]);
ylim(100*[-1 1]);
ylabel('Spatial Frequency, $\xi_x$ ($m^{-1}$)','interpreter','latex');
subplot(3,1,3);
surf(freq2{2},(freq2{2}'./freq2{1})',log10(sxx2),'linestyle','none','facecolor','interp');
view(2);
xlim([1e2 2e4]);
ylim(400*[-1 1]);
xlabel('Temporal Frequency, $f$ ($Hz$)','interpreter','latex');
ylabel('Assumed Velocity ($m/s$)','interpreter','latex');
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Children(3).TickLabelInterpreter = 'latex';
f1.Children(1).XScale = 'log';
f1.Children(2).XScale = 'log';
f1.Children(1).Layer = 'top';
f1.Children(2).Layer = 'top';
f1.Units = 'inches';
f1.Position = [1 1 5.5 7.5];

% saveas(f1,'dispersion_demo.eps','epsc');