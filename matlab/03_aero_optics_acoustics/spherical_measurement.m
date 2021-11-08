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
    load('spherical_sample_win.mat');
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
%     A0wf = opdrms/calcOPDA0(R,Lambda,Ap)/1e6;
    
%     disp(num2str(Lambda/Ap,'%0.2f'));
%     disp(num2str(R/Ap,'%0.2f'));
%     disp(' ');
    fprintf(fileID,[num2str(runStats(1)) ' & ' num2str(runStats(2)) ' & ' num2str(prms,'%0.2f') ' & ' num2str(opdrms,'%0.3e') ' & ' num2str(A0mic,'%0.2f') ' & ' num2str(A0wf,'%0.2f') ' & ' num2str(diff([A0mic A0wf])/mean([A0mic A0wf])*100,'%0.2f') '\\\\ \n']);
end
fprintf(fileID,'\\hline \n');
fprintf(fileID,'\\end{tabular}\n');

% f1 = figure(1);
% loglog(MicFreq,MicSxx);
% grid on;


function [opd_a0] = calcOPDA0(R,LAMBDA,AP)
nLenslets = 32;
zMax = 5;
dz = 0.05;
c0 = 340;
kgd = 2.27e-4;
phaseSteps = 25;

[x,y,z] = meshgrid(0.975*AP*(-0.5:1/(nLenslets-1):0.5),0.975*AP*(-0.5:1/(nLenslets-1):0.5),-zMax:dz:zMax);
r = sqrt(x(:,:,1).^2+y(:,:,1).^2);
t = atan2(y(:,:,1),x(:,:,1));
mask = ones(nLenslets);
mask(r>0.5*AP) = NaN;
phase = (0:phaseSteps-1)/phaseSteps*2*pi;

k = 2*pi/LAMBDA;
R2 = sqrt((x+R).^2+y.^2+z.^2);
OPD = kgd/c0^2*sum(1./R2.*exp(-1i*k*R2),3)*dz;
stats = zeros(1,phaseSteps);
for cc=1:phaseSteps
    wf = mask.*real(OPD.*exp(1i*phase(cc)));
    wf = ZernikeRemoval(1:3,wf,r.*mask,t.*mask,'noll');
    stats(cc) = nanrms(wf(:));
end
opd_a0 = mean(stats);

end