close all; clc; clearvars;

load('synthetic_wavefront.mat');

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


f1 = figure(1);
semilogx(highpassFilters,wf1/wf20,highpassFilters,wf2./wf20);
grid on;
xlim([100 10000]);
% ylim([0 1]);
xlabel('High-Pass Filter Cut-Off Frequency (Hz)','interpreter','latex');
ylabel('$OPD_{F,rms}/OPD_{AO,rms}$','interpreter','latex');
legend('Total Wavefront','AO-Only Wavefront','interpreter','latex');
f1.Units = 'inches';
f1.Position = [2 2 6 3.25];
f1.Children(end).TickLabelInterpreter = 'latex';

saveas(f1,'filter_temporal.eps','epsc');