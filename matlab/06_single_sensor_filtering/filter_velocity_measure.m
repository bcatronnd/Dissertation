close all; clc; clearvars;

load('../05_synthetic_wavefront/synthetic_wavefront.mat');

c = 340;
M = 0.6;
uBL_u = 0.7:0.025:0.9;

BlockSize = 2^10;
[WF,freq] = simpleDispersion(wf.wf,'BlockSize',BlockSize,'SampleRate',wf.sampleRate,'Output','lin');

BlockSize = size(WF);
Frequency.y = reshape(-1/2:1/BlockSize(1):1/2-1/BlockSize(1),[],1);
Frequency.x = reshape(-1/2:1/BlockSize(2):1/2-1/BlockSize(2),1,[]);
Frequency.t = reshape(-1/2:1/BlockSize(3):1/2-1/BlockSize(3),1,1,[]);
Frequency.rho = sqrt(Frequency.x.^2+Frequency.y.^2);

u = c*M*uBL_u;
v = 0;

width = 0.0125;
order = 1;
% High-Pass Rho
gain = sqrt(1./(1+(Frequency.rho/0.1).^(-2*2)));
WF = WF.*gain.^2;
for aa=1:length(u)
    for bb=1:length(v)
        U = u(aa)*wf.sampleRate(2)/wf.sampleRate(3);
        V = v(bb)*wf.sampleRate(1)/wf.sampleRate(3);
        dist = abs(U*Frequency.x+V*Frequency.y-Frequency.t)/sqrt(U^2+V^2+1);
        gain = sqrt(1./(1+(dist/width).^(+2*order)));
%         wf(aa,bb) = mean(nanrms(reshape(real(ifftn(ifftshift(WF.*gain))),size(WF,1)*size(WF,2),[])));
        wf_f(aa,bb) = sum(WF.*gain.^2,'all');
    end
end

[wfm,id] = max(wf_f,[],'all','linear');
[idr,idc] = ind2sub(size(wf_f),id);

%%%%% Plot
f1 = figure(1);
plot(uBL_u,wf_f/wfm,'k-o');
grid on;
xlabel('$u_{BL}/U$','interpreter','latex');
ylabel('Normalized Power','interpreter','latex');
xlim([0.7 0.9]);
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 5.5 3.5];

saveas(f1,'filter_velocity_measure.eps','epsc');

