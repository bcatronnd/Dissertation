close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/ProcessedWavefronts/';
testPoint = {'20210417002' '20210417003' '20210417004' '20210417005' '20210417006' '20210417007' '20210417008'};

load([directory '20210417009' '_RL.mat']);
MicHeight = str2double(extract(RunLog.notes{3},digitsPattern))/1e3;
BeamHeight = str2double(extract(RunLog.notes{4},digitsPattern))/1e3;

fileID = fopen('spherical_measurement.txt','w');
fprintf(fileID,'\\begin{tabular}{c c c c c c c}\n');
fprintf(fileID,'\\hline \n');
fprintf(fileID,'$f_{speaker}$ & $V_{speaker}$ & $p_{rms}$ & OPD$_{RMS}$ & $|A_0|_{mic}$ & $|A_0|_{wf}$ & Diff \\\\ \n');
fprintf(fileID,'(Hz) & (mV) & (Pa) & ($\\mu m$) & ($kg/s^2$) & ($kg/s^2$) & (\\%%) \\\\ \n');
fprintf(fileID,'\\hline \n');

for aa=1:length(testPoint)
    [WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint{aa} '_WF.mat'],'Scale',1e6);%,'ZernikeRemoval',1:3);
    load([directory testPoint{aa} '_CDAQ.mat']);
    load('spherical_sample.mat');
    runStats = str2double(extract(RunLog.notes{2},digitsPattern));
    Lambda = RunLog.c/runStats(1);
    
    % Mic
    [MicSxx,MicFreq] = computeSXX(1e3*data(:,2)','blocksize',2^12,'positiveonly',1,'samplerate',1/timeStamps(2),'window',@hann);
    [~,MicInd] = min(abs(MicFreq-runStats(1)));
    MicMax = MicSxx(MicInd);
    prms = sqrt(MicMax*diff(MicFreq(1:2)));
    A0mic = prms*sqrt(2)*MicHeight;
    
    % WF
    R = mean(BeamHeight);
    Ap = diff(BeamHeight);
    opdrms = mean(nanrms(reshape(WFfilter(WF,'time-bandpass',(runStats(1)+10*[-1 1])/RunLog.samplerate(3)),prod(size(WF,[1 2])),size(WF,3)),1));
    A0wf = opdrms*sqrt(R/Ap)/fitresult(Lambda/Ap);
    
    
    fprintf(fileID,[num2str(runStats(1)) ' & ' num2str(runStats(2)) ' & ' num2str(prms,'%0.2f') ' & ' num2str(opdrms,'%0.3e') ' & ' num2str(A0mic,'%0.2f') ' & ' num2str(A0wf,'%0.2f') ' & ' num2str(diff([A0mic A0wf])/mean([A0mic A0wf])*100,'%0.2f') '\\\\ \n']);
end
fprintf(fileID,'\\hline \n');
fprintf(fileID,'\\end{tabular}\n');

% f1 = figure(1);
% loglog(MicFreq,MicSxx);
% grid on;


