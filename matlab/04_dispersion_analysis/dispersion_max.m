close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[wf,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);

blockSize = 2^10;
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

%%
df = diff(freq2{1}(1:2));
f1 = figure(1);
loglog(freq1,sxx1,freq2{2},sampleRate(1)*max(sxx2,[],1),freq2{2},trapz(freq2{1},sxx2,1));
grid on;
xlim([1e2 2e4]);
% ylim(10.^[-10 -6]);
xlabel('Temporal Frequency, $f$ ($Hz$)','interpreter','latex');
ylabel('$S_{xx}$ ($\mu m^2/Hz$)','interpreter','latex');
legend('$S_{xx}^1$','max($S_{xx}^2$)$f_{s,x}$','$\int S_{xx}^2df_{s,x}$','interpreter','latex','location','southwest');
f1.Units = 'inches';
f1.Position = [1 1 5.5 3.25];
f1.Children(2).TickLabelInterpreter = 'latex';

saveas(f1,'dispersion_max.eps','epsc');
% saveas(f1,'dispersion_max.png','png');
clc;
disp(['Error (max): ' num2str(sum((sxx1-sampleRate(1)*max(sxx2(:,end/2+1:end),[],1)).^2))]);
disp(['Error (sum): ' num2str(sum((sxx1-sampleRate(1)*sum(sxx2(:,end/2+1:end),1)).^2))]);
disp(['Error (rms): ' num2str(sum((sxx1-sampleRate(1)*rms(sxx2(:,end/2+1:end),1)).^2))]);
disp(['Error (mean): ' num2str(sum((sxx1-sampleRate(1)*mean(sxx2(:,end/2+1:end),1)).^2))]);