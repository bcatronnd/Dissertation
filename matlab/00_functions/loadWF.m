function [WF,CineInfo,RunLog,WFInfo] = loadWF(filename,varargin)

% Remove '_WF.mat' From Filename If It Exists
if strcmp(filename(end-6:end),'_WF.mat')
    filename = filename(1:end-7);
end
% Default Options
option.permute = [];
option.scale = 1;
option.zernike = [];
option.reaperature = [];
option.lensletnoise = 0;
% Parse Varargin Inputs
for aa=1:2:nargin-1
    switch lower(varargin{aa})
        case 'permute'
            option.permute = varargin{aa+1};
        case 'scale'
            option.scale = varargin{aa+1};
        case 'zernikeremoval'
            option.zernike = varargin{aa+1};
        case 'reaperature'
            option.reaperature = varargin{aa+1};
        case 'lensletnoise'
            option.lensletnoise = varargin{aa+1};
    end
end
clear aa;
% Load Wavefront File
load([filename '_WF.mat'],'CineInfo','Mask_WF','WF','WFmean');
if ~exist('CineInfo','var')
    CineInfo = [];
end
if exist('Mask_WF','var')
    WFInfo.Mask_WF = Mask_WF;
else
    WFInfo.Mask_WF = rms(WF,3)>(rms(WF(:))*1e-3);
end
if exist('WFmean','var')
    WFInfo.WFmean = WFmean;
else
    WFInfo.WFmean = mean(WF,3);
end
clear Mask_WF WFmean;
% Check For Run Log File And Load
if exist([filename '_RL.mat'],'file')
    load([filename '_RL.mat'],'RunLog');
else
    RunLog = [];
end
% NaN Mask and Apply
WFInfo.Mask_WF(WFInfo.Mask_WF==0) = NaN;
WF = WF.*WFInfo.Mask_WF;
WFInfo.WFmean = WFInfo.WFmean.*WFInfo.Mask_WF;
% Remove NaN Rows and Columns
WF(~any(WF,2),:,:) = [];
WF(:,~any(WF,1),:) = [];
WFInfo.Mask_WF(~any(WFInfo.Mask_WF,2),:) = [];
WFInfo.Mask_WF(:,~any(WFInfo.Mask_WF,1)) = [];
WFInfo.WFmean(~any(WFInfo.WFmean,2),:) = [];
WFInfo.WFmean(:,~any(WFInfo.WFmean,1)) = [];
% Permute And Flip
if isempty(option.permute) && isfield(RunLog,'permute')
    option.permute = RunLog.permute;
end
if ~isempty(option.permute)
    WF = permute(WF,abs(option.permute));
    WFInfo.Mask_WF = permute(WFInfo.Mask_WF,abs(option.permute));
    WFInfo.WFmean = permute(WFInfo.WFmean,abs(option.permute));
    if any(option.permute<0)
        option.flip = (option.permute<0).*(1:length(option.permute));
        option.flip(option.flip==0) = [];
        for aa=1:length(option.flip)
            WF = flip(WF,option.flip(aa));
            WFInfo.Mask_WF = flip(WFInfo.Mask_WF,option.flip(aa));
            WFInfo.WFmean = flip(WFInfo.WFmean,option.flip(aa));
        end
    end
end
% Calculate Rho and Theta
[WFInfo.rho,WFInfo.theta] = ZernikeRT(WFInfo.Mask_WF,0.975);
% Reaperature
if ~isempty(option.reaperature)
    WF(WFInfo.rho>option.reaperature) = NaN;
    WFInfo.Mask_WF(WFInfo.rho>option.reaperature) = NaN;
    WFInfo.WFmean(WFInfo.rho>option.reaperature) = NaN;
    WFInfo.theta(WFInfo.rho>option.reaperature) = NaN;
    WFInfo.rho(WFInfo.rho>option.reaperature) = NaN;
    % Remove NaN Rows and Columns
    WF(~any(WF,2),:,:) = [];
    WF(:,~any(WF,1),:) = [];
    WFInfo.Mask_WF(~any(WFInfo.Mask_WF,2),:) = [];
    WFInfo.Mask_WF(:,~any(WFInfo.Mask_WF,1)) = [];
    WFInfo.WFmean(~any(WFInfo.WFmean,2),:) = [];
    WFInfo.WFmean(:,~any(WFInfo.WFmean,1)) = [];
    WFInfo.theta(~any(WFInfo.theta,2),:) = [];
    WFInfo.theta(:,~any(WFInfo.theta,1)) = [];
    WFInfo.rho(~any(WFInfo.rho,2),:) = [];
    WFInfo.rho(:,~any(WFInfo.rho,1)) = [];
end
% Zernike Removal
if ~isempty(option.zernike)
    WF = ZernikeRemoval(option.zernike,WF,WFInfo.rho,WFInfo.theta,'noll');
end
% Scale Output
WF = WF*option.scale;
% Lenslet Noise
if option.lensletnoise
    opdrms = rms(WF,3);
    opdrms_avg = mean(opdrms,'all','omitnan');
    opdrms_std = std(opdrms,[],'all','omitnan');
    WFInfo.MaskNoise = double(opdrms<(opdrms_avg+option.lensletnoise*opdrms_std));
    WFInfo.MaskNoise(WFInfo.MaskNoise==0) = NaN;
    WF = WF.*WFInfo.MaskNoise;
end
end