close all; clc; clearvars;

load('../05_synthetic_wavefront/synthetic_wavefront.mat');
PlotColors = linspecer(3);

wf0 = mean(nanrms(reshape(wf.wf,size(wf.wf,1)*size(wf.wf,2),size(wf.wf,3))),2);
wf20 = mean(nanrms(reshape(wf.AO,size(wf.AO,1)*size(wf.AO,2),size(wf.AO,3))),2);
highpassFilters = logspace(2,4);
for aa=1:length(highpassFilters)
    WF1 = WFfilter(wf.wf,'time-highpass',highpassFilters(aa)/wf.sampleRate(3));
    WF2 = WFfilter(wf.AO,'time-highpass',highpassFilters(aa)/wf.sampleRate(3));
    wf1(aa) = mean(nanrms(reshape(WF1,size(WF1,1)*size(WF1,2),size(WF1,3))),2);
    wf2(aa) = mean(nanrms(reshape(WF2,size(WF2,1)*size(WF2,2),size(WF2,3))),2);
    clear WF1 WF2;
end
% 
BlockSize = 2^10;
window = permute(hann(BlockSize),[3 2 1]);
cw = 1/sqrt(sum(window.^2,'all')/BlockSize);
freq = (-0.5:1/BlockSize:0.5-1/BlockSize)*wf.sampleRate(3);
wf.wf = wf.wf(:,:,1:floor(size(wf.wf,3)/BlockSize)*BlockSize);
wf.AO = wf.AO(:,:,1:floor(size(wf.AO,3)/BlockSize)*BlockSize);
wf.wf = reshape(wf.wf,size(wf.wf,1),size(wf.wf,2),BlockSize,size(wf.wf,3)/BlockSize);
wf.AO = reshape(wf.AO,size(wf.AO,1),size(wf.AO,2),BlockSize,size(wf.AO,3)/BlockSize);
wf.wf = permute(mean(abs(cw*fftshift(fft(wf.wf.*window,[],3))).^2/BlockSize/wf.sampleRate(3),[1 2 4],'omitnan'),[1 3 2]);
wf.AO = permute(mean(abs(cw*fftshift(fft(wf.AO.*window,[],3))).^2/BlockSize/wf.sampleRate(3),[1 2 4],'omitnan'),[1 3 2]);





%%
close all;
f1 = figure(1);
subplot(2,1,1);
semilogx(highpassFilters,wf1/wf20,'linewidth',1.25,'color',PlotColors(1,:));
hold on;
semilogx(highpassFilters,wf2./wf20,'linewidth',1.25,'color',PlotColors(2,:));
grid on;
xlim([100 10000]);
% ylim([0 1]);
xlabel('High-Pass Filter Cutoff Frequency (Hz)','interpreter','latex');
ylabel('$OPD_{RMS,F}/OPD_{RMS,AO}$','interpreter','latex');
legend('Total Wavefront','AO-Only Wavefront','interpreter','latex');
f1.Children(2).TickLabelInterpreter = 'latex';

subplot(2,1,2);
loglog(freq,wf.wf,'linewidth',1.25,'color',PlotColors(1,:));
hold on;
loglog(freq,wf.AO,'linewidth',1.25,'color',PlotColors(2,:));
grid on;
xlim([100 10000]);
xlabel('Frequency (Hz)','interpreter','latex');
ylabel('$S_{xx}$ ($\mu m^2/Hz$)','interpreter','latex');
f1.Children(1).TickLabelInterpreter = 'latex';

f1.Units = 'inches';
f1.Position = [2 2 6 4.5];
saveas(f1,'filter_temporal.eps','epsc');