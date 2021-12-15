close all; clc; clearvars;

% Inputs
Directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
TestPoint = {'20210901001' '20210901002' '20210901003'};
BlockLength = 2^12;
OverlapFactor = 4;
SensorSelection{1} = 1:2;
SensorSelection{2} = 3:6;
SensorSelection{3} = 8:13;
SensorSelection{4} = [7 14:16];
SampleRate = 49000;
wfSize = 64000;

StepSize = BlockLength*2^-OverlapFactor;
BlockIndex = 1+(0:StepSize:wfSize-BlockLength);
BlockIndex(2,:) = BlockIndex(1,:)+BlockLength-1;
BlockNumber = size(BlockIndex,2);
FrequencyN = permute((-0.5:1/BlockLength:0.5-1/BlockLength),[1 3 2]);
Frequency = SampleRate*FrequencyN;
% Windowing Functions
Window.y.w = reshape(hann(BlockLength),1,BlockLength);
Window.y.cw = 1/sqrt(sum(Window.y.w.^2,'all')/numel(Window.y.w));
PlotColors = linspecer(length(TestPoint));

f1 = figure(1);
for aa=1:length(TestPoint)
    load([Directory TestPoint{aa} '_CDAQ.mat'],'scanData');
    load([Directory TestPoint{aa} '_RL.mat']);
    blockageRatio = 34.462/36^2;
    RunLog.u = RunLog.u*(1+blockageRatio);
    % Sensor PSD
    for bb=1:length(SensorSelection)
        scanData2 = scanData(1:BlockIndex(2,end),SensorSelection{bb})';
        Y = complex(zeros(size(scanData2,1),BlockLength,BlockNumber));
        for cc=1:BlockNumber
            Y(:,:,cc) = fftshift(fft(Window.y.w.*scanData2(:,BlockIndex(1,cc):BlockIndex(2,cc)),BlockLength,2),2);
        end
        Y = mean(abs(Window.y.cw*Y).^2,3)/BlockLength/SampleRate;
        clear cc scanData2;
        
        subplot(length(SensorSelection),1,bb);
        loglog(squeeze(Frequency)/RunLog.u,Y,'color',PlotColors(aa,:));
        hold on;
        grid on;
    end
end
subplot(4,1,1);
title('Ambient Microphones','interpreter','latex');
xlim(10.^[0 2]);
legend('M=0.3','M=0.4','M=0.5','interpreter','latex');
subplot(4,1,2);
title('Test-Section Microphones','interpreter','latex');
xlim(10.^[0 2]);
subplot(4,1,3);
title('Test-Section Accelerometers','interpreter','latex');
xlim(10.^[0 2]);
subplot(4,1,4);
title('Other Accelerometers','interpreter','latex');
xlim(10.^[0 2]);
xlabel('St/l ($m^{-1}$)','interpreter','latex');
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Children(3).TickLabelInterpreter = 'latex';
f1.Children(5).TickLabelInterpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 5.5 6];

saveas(f1,'sensor_spectra.eps','epsc');

