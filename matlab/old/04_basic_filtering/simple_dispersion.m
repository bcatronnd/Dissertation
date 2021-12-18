close all; clc; clearvars;


u = [250 -150]; % Velocity Vector (m/s)
g0 = [1 0.5];
fcut = [6000 3000];
fmax = 20000;

order = 4;

fx = 200;       % Samples/m
ft = 15000;     % Samples/sec

nx = 2^7;       % Number of Points in x
nxp = 2^6;
nt = 2^12;      % Number of Points in t
ntp = 2^10;

fstep = ft/nt/2;
%
[t,x] = meshgrid((0:nt-1)/ft,(0:nx-1)/fx);
opd = (2*rand(size(x))-1)*1e2;
freq = 0:fstep:fmax;
for aa=1:length(u)
%     opd = opd+sum(sqrt(g0(aa)/(1+(freq/fcut(aa)).^(2*order))).*real(exp(1i*2*pi*freq/u(aa).*x).*exp(1i*2*pi*freq.*t+2*pi*rand(size(freq)))),3);
    for bb=1:length(freq)
        opd = opd+sqrt(g0(aa)/(1+(freq(bb)/fcut(aa))^(2*order)))*real(exp(1i*2*pi*freq(bb)/u(aa)*x).*exp(1i*2*pi*freq(bb)*t+2*pi*rand(1)));
    end
end
opd = imresize(opd,[nxp nt]);
opd = reshape(opd,nxp,ntp,[]);


OPD = fftshift(mean(abs((fft(fft(opd,[],1),[],2)).^2),3));
f1 = figure;
surf((-ntp/2:ntp/2-1)/ntp,(-nxp/2:nxp/2-1)/nxp,log10(OPD),'linestyle','none');%,'facecolor','interp');
view(2);
xlim([0 0.5]);
ylim([-0.4 0.4]);
% colorbar;
caxis(log10([mean(OPD,'all') max(OPD,[],'all')]));
xlabel('Normalized Temporal Frequency','interpreter','latex');
ylabel('Normalized Spacial Frequency','interpreter','latex');
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Children(1).Layer = 'top';
f1.Units = 'inches';
f1.Position = [1 1 5.5 3.5];

f2 = figure;
freqt = (-ntp/2:ntp/2-1)/ntp;
freqx = (-nxp/2:nxp/2-1)'/nxp;
surf(freqt,freqt./freqx,log10(OPD),'linestyle','none','facecolor','interp');
view(2);
xlim([1e-2 0.5]);
ylim(2.5*[-1 1]);
caxis(log10([mean(OPD,'all') max(OPD,[],'all')]));
xlabel('Normalized Temporal Frequency','interpreter','latex');
ylabel({'Normalized Velocity' '(Temporal Freq./Spatial Freq.)'},'interpreter','latex');
f2.Children(1).TickLabelInterpreter = 'latex';
f2.Children(1).Layer = 'top';
f2.Units = 'inches';
f2.Position = [1 1 5.5 3.5];
f2.Children(1).XScale = 'log';

if isempty(mfilename)
    saveas(f1,'simple_dispersion.eps','epsc');
    saveas(f2,'simple_dispersion_vel.eps','epsc');
else
    saveas(f1,[mfilename '.eps'],'epsc');
    saveas(f2,[mfilename '_vel.eps'],'epsc');
end