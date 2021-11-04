close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);

BlockSize = 2^12;
WF = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));
[WF,freq] = simpleDispersion(WF,'BlockSize',BlockSize,'SampleRate',RunLog.samplerate,'Output','lin');


BlockSize = size(WF);
Frequency.y = reshape(-1/2:1/BlockSize(1):1/2-1/BlockSize(1),[],1);
Frequency.x = reshape(-1/2:1/BlockSize(2):1/2-1/BlockSize(2),1,[]);
Frequency.t = reshape(-1/2:1/BlockSize(3):1/2-1/BlockSize(3),1,1,[]);
Frequency.rho = sqrt(Frequency.x.^2+Frequency.y.^2);

blockageRatio = 34.462/36^2;
RunLog.u = RunLog.u*(1+blockageRatio);

u = 150:1:155;
v = 0:1:5;

width = 0.025;
order = 2;
% High-Pass Rho
gain = sqrt(1./(1+(Frequency.rho/0.1).^(-2*2)));
WF = WF.*gain.^2;
for aa=1:length(u)
    for bb=1:length(v)
        U = u(aa)*RunLog.samplerate(2)/RunLog.samplerate(3);
        V = v(bb)*RunLog.samplerate(1)/RunLog.samplerate(3);
        dist = abs(U*Frequency.x+V*Frequency.y-Frequency.t)/sqrt(U^2+V^2+1);
        gain = sqrt(1./(1+(dist/width).^(+2*order)));
%         wf(aa,bb) = mean(nanrms(reshape(real(ifftn(ifftshift(WF.*gain))),size(WF,1)*size(WF,2),[])));
        wf(aa,bb) = sum(WF.*gain.^2,'all');
    end
end

[wfm,id] = max(wf,[],'all','linear');
[idr,idc] = ind2sub(size(wf),id);
%%
close all;
f1 = figure(1);
mesh(v,u,wf/wfm,'FaceAlpha','0.5','FaceColor','interp','edgecolor','none');
hold on;
plot3(v(idc),u(idr),1,'ko','markerfacecolor','k','markersize',10);
view(45,45);
grid on;
xlabel('v (m/s)','interpreter','latex');
ylabel('u (m/s)','interpreter','latex');
zlabel('Normalized Power','interpreter','latex');
f1.Children(1).TickLabelInterpreter = 'latex';

diary filter_velocity_real.txt;
diary on;
disp(['u: ' num2str(u(idr)) 'm/s - ' num2str(u(idr)/RunLog.u,'%0.2f') 'u' char(8734)]);
disp(['v: ' num2str(v(idc)) 'm/s - ' num2str(v(idc)/RunLog.u,'%0.2f') 'u' char(8734)]);
disp(['Flow Angle: ' num2str(atand(v(idc)/u(idr)),'%0.1f') char(176)]);
diary off;

saveas(f1,'filter_velocity_real.eps','epsc');
saveas(f1,'filter_velocity_real.png','png');