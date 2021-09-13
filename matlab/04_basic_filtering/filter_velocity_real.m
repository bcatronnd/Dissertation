close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/ProcessedWavefronts/';
filename = '20180104022';
[WF,CineInfo,RunLog,WFInfo] = loadWF([directory filename '_WF.mat'],'Permute',[2 1 3],'Scale',1e6);

BlockSize = 2^10;
WF = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));
[WF,freq] = simpleDispersion(WF,'BlockSize',BlockSize,'SampleRate',RunLog.samplerate,'Output','lin');


BlockSize = size(WF);
Frequency.y = reshape(-1/2:1/BlockSize(1):1/2-1/BlockSize(1),[],1);
Frequency.x = reshape(-1/2:1/BlockSize(2):1/2-1/BlockSize(2),1,[]);
Frequency.t = reshape(-1/2:1/BlockSize(3):1/2-1/BlockSize(3),1,1,[]);
Frequency.rho = sqrt(Frequency.x.^2+Frequency.y.^2);



u = 190:1:220;
v = -25:1:-10;

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
zlabel('Normalized Energy','interpreter','latex');
f1.Children(1).TickLabelInterpreter = 'latex';

disp(['u: ' num2str(u(idr))]);
disp(['v: ' num2str(v(idc))]);

saveas(f1,'filter_velocity_real.eps','epsc');