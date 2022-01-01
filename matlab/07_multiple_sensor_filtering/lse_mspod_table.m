close all; clc; clearvars;

% Inputs
Directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
TestPoint = {'20210901001' '20210901002' '20210901003'};
BlockLength = 2^10;
OverlapFactor = 0;
SensorSelection{1} = 7:16;
SensorSelection{2} = 1:2;
SensorSelection{3} = 3:6;
SensorSelection{4} = [3:6 8:13];
SensorSelection{5} = 1:16;

DataSummary = zeros(4+length(SensorSelection),length(TestPoint));
DataSummaryLabels = {'Original','Tip/Tilt Removal','Accelerometers','Ambient Microphones','Duct Microphones','Test-Section Sensors','All Sensors','Velocity Filter','+All Sensors'};





for aa=1:length(TestPoint)
    % Load Data
    [wf,CineInfo,RunLog,WFInfo] = loadWF([Directory TestPoint{aa} '_WF.mat'],'Scale',1e6);
    load([Directory TestPoint{aa} '_CDAQ.mat'],'scanData');
    
    % OPDrms - Original
    wfSize = size(wf);
    DataSummary(1,aa) = mean(nanrms(reshape(wf,prod(wfSize(1:2)),wfSize(3))));
    
    % OPDrms - Z123 Removal
    wf = ZernikeRemoval(1:3,wf,WFInfo.rho,WFInfo.theta,'noll');
    DataSummary(2,aa) = mean(nanrms(reshape(wf,prod(wfSize(1:2)),wfSize(3))));
    
    % Basic Setup
    StepSize = BlockLength*2^-OverlapFactor;
    BlockIndex = 1+(0:StepSize:size(wf,3)-BlockLength);
    BlockIndex(2,:) = BlockIndex(1,:)+BlockLength-1;
    BlockSize = [size(wf,[1 2]) BlockLength];
    BlockNumber = size(BlockIndex,2);
    DataLocation = reshape(WFInfo.Mask_WF==1,prod(BlockSize(1:2)),1);
    LensletNumber = prod(BlockSize(1:2));
    Frequency.yn = (-0.5:1/BlockSize(1):0.5-1/BlockSize(1))';
    Frequency.xn = (-0.5:1/BlockSize(2):0.5-1/BlockSize(2));
    Frequency.tn = permute((-0.5:1/BlockSize(3):0.5-1/BlockSize(3)),[1 3 2]);
    Frequency.y = Frequency.yn*RunLog.samplerate(1);
    Frequency.x = Frequency.xn*RunLog.samplerate(2);
    Frequency.t = permute(Frequency.tn*RunLog.samplerate(3),[1 3 2]);
    
    % Windowing Functions
    Window.wf.s = createSpatialWindow(WFInfo.Mask_WF);
    Window.wf.t = reshape(hann(BlockLength),1,1,BlockLength);
    Window.wf.w = Window.wf.s.*Window.wf.t;
    Window.wf.cw = 1/sqrt(sum(Window.wf.w.^2,'all')/numel(Window.wf.w));
    Window.y.w = reshape(hann(BlockLength),1,BlockLength);
    Window.y.cw = 1/sqrt(sum(Window.y.w.^2,'all')/numel(Window.y.w));
    
    % Wavefront PSD In Time
    wf = wf(:,:,1:BlockIndex(2,end));
    wf(isnan(wf)) = 0;
    WF = complex(zeros([BlockSize,BlockNumber]));
    for bb=1:BlockNumber
        WF(:,:,:,bb) = Window.wf.cw*fftshift(fftn(Window.wf.w.*wf(:,:,BlockIndex(1,bb):BlockIndex(2,bb))));
    end
    WF = reshape(WF,prod(BlockSize(1:2)),BlockSize(3),BlockNumber);
    WF = permute(WF,[1 3 2]);
    clear bb;
    
    for bb=1:length(SensorSelection)
        % Sensor PSD
        scanData2 = scanData(1:BlockIndex(2,end),SensorSelection{bb})';
        Y = complex(zeros(size(scanData2,1),BlockLength,BlockNumber));
        for cc=1:BlockNumber
            Y(:,:,cc) = fftshift(fft(Window.y.w.*scanData2(:,BlockIndex(1,cc):BlockIndex(2,cc)),BlockLength,2),2);
        end
        Y = permute(Y,[1 3 2]);
        clear cc scanData2;
        
        % LSE-SPOD
        WFfiltered = complex(zeros(size(WF)));
        for cc=1:BlockLength
            Q = WF(:,:,cc);
            y = Y(:,:,cc);
            [PSI,~] = eig(Q'*Q,'vector');
            PHI = Q*PSI;
            L = (PSI*y')/(y*y');
            PSIlse = L*y;
            PSIao = PSI-PSIlse;
            WFfiltered(:,:,cc) = (PSIao*PHI')';
        end
        clear cc PSI PHI PSIao PSIlse Q y;
        
        % Inverse FFT
        WFfiltered = ipermute(WFfiltered,[1 3 2]);
        WFfiltered = reshape(WFfiltered,[BlockSize BlockNumber]);
        for cc=1:BlockNumber
            WFfiltered(:,:,:,cc) = ifftn(ifftshift(WFfiltered(:,:,:,cc)));
        end
        DataSummary(2+bb,aa) = mean(nanrms(reshape(real(WFfiltered).*WFInfo.Mask_WF,prod(BlockSize(1:2)),BlockSize(3)*BlockNumber)));
        clear Y WFfiltered;
    end
    
    % Velocity Filter
    wf = WFfilter(wf,'velocity-lowpass',[RunLog.u*0.85*RunLog.samplerate(2)/RunLog.samplerate(3),RunLog.u*0.02*RunLog.samplerate(1)/RunLog.samplerate(3),2,0.025]);
    DataSummary(end-1,aa) = mean(nanrms(reshape(wf.*WFInfo.Mask_WF,prod(wfSize(1:2)),[])));
    
    % Velocity+Sensor
    WF = complex(zeros([BlockSize,BlockNumber]));
    for bb=1:BlockNumber
        WF(:,:,:,bb) = Window.wf.cw*fftshift(fftn(Window.wf.w.*wf(:,:,BlockIndex(1,bb):BlockIndex(2,bb))));
    end
    WF = reshape(WF,prod(BlockSize(1:2)),BlockSize(3),BlockNumber);
    WF = permute(WF,[1 3 2]);
    clear bb;
    % Sensor PSD
    scanData2 = scanData(1:BlockIndex(2,end),SensorSelection{end})';
    Y = complex(zeros(size(scanData2,1),BlockLength,BlockNumber));
    for cc=1:BlockNumber
        Y(:,:,cc) = fftshift(fft(Window.y.w.*scanData2(:,BlockIndex(1,cc):BlockIndex(2,cc)),BlockLength,2),2);
    end
    Y = permute(Y,[1 3 2]);
    clear cc scanData2;
    % LSE-SPOD
    WFfiltered = complex(zeros(size(WF)));
    for cc=1:BlockLength
        Q = WF(:,:,cc);
        y = Y(:,:,cc);
        [PSI,~] = eig(Q'*Q,'vector');
        PHI = Q*PSI;
        L = (PSI*y')/(y*y');
        PSIlse = L*y;
        PSIao = PSI-PSIlse;
        WFfiltered(:,:,cc) = (PSIao*PHI')';
    end
    clear cc PSI PHI PSIao PSIlse Q y;
    % Inverse FFT
    WFfiltered = ipermute(WFfiltered,[1 3 2]);
    WFfiltered = reshape(WFfiltered,[BlockSize BlockNumber]);
    for cc=1:BlockNumber
        WFfiltered(:,:,:,cc) = ifftn(ifftshift(WFfiltered(:,:,:,cc)));
    end
    DataSummary(end,aa) = mean(nanrms(reshape(real(WFfiltered).*WFInfo.Mask_WF,prod(BlockSize(1:2)),BlockSize(3)*BlockNumber)));
    clear Y WFfiltered;
end

save('lse_mspod_table.mat','DataSummary','DataSummaryLabels','SensorSelection','TestPoint');

% fileID = fopen('lse_mspod_table.txt','w');
% fprintf(fileID,'\\begin{tabular}{c c c c}\n');
% fprintf(fileID,' & M=0.3 & M=0.4 & M=0.5 \\\\ \n');
% fprintf(fileID,'\\hline \\hline \n');
% fprintf(fileID,[DataSummaryLabels{1} ' & ' num2str(DataSummary(1,1),'%0.4f') ' & ' num2str(DataSummary(1,2),'%0.4f') ' & ' num2str(DataSummary(1,3),'%0.4f') ' \\\\ \n']);
% fprintf(fileID,[DataSummaryLabels{2} ' & ' num2str(DataSummary(2,1),'%0.4f') '(' num2str((1-DataSummary(2,1)/DataSummary(1,1))*100,'%0.1f') '\\%%) & ' num2str(DataSummary(2,2),'%0.4f') '(' num2str((1-DataSummary(2,2)/DataSummary(1,2))*100,'%0.1f') '\\%%) & ' num2str(DataSummary(2,3),'%0.4f') '(' num2str((1-DataSummary(2,3)/DataSummary(1,3))*100,'%0.1f') '\\%%) \\\\ \n']);
% fprintf(fileID,'\\hline \n');
% for aa=3:size(DataSummary,1)-2
%     fprintf(fileID,[DataSummaryLabels{aa} ' & ' num2str(DataSummary(aa,1),'%0.4f') '(' num2str((1-DataSummary(aa,1)/DataSummary(1,1))*100,'%0.1f') '\\%%) & ' num2str(DataSummary(aa,2),'%0.4f') '(' num2str((1-DataSummary(aa,2)/DataSummary(1,2))*100,'%0.1f') '\\%%) & ' num2str(DataSummary(aa,3),'%0.4f') '(' num2str((1-DataSummary(aa,3)/DataSummary(1,3))*100,'%0.1f') '\\%%) \\\\ \n']);
% end
% fprintf(fileID,'\\hline \n');
% for aa=size(DataSummary,1)-1:size(DataSummary,1)
%     fprintf(fileID,[DataSummaryLabels{aa} ' & ' num2str(DataSummary(aa,1),'%0.4f') '(' num2str((1-DataSummary(aa,1)/DataSummary(1,1))*100,'%0.1f') '\\%%) & ' num2str(DataSummary(aa,2),'%0.4f') '(' num2str((1-DataSummary(aa,2)/DataSummary(1,2))*100,'%0.1f') '\\%%) & ' num2str(DataSummary(aa,3),'%0.4f') '(' num2str((1-DataSummary(aa,3)/DataSummary(1,3))*100,'%0.1f') '\\%%) \\\\ \n']);
% end
% fprintf(fileID,'\\end{tabular}\n');


