close all; clc; clearvars; %#ok<*UNRCH>

sampleRate = [200*[1 1] 30000];
nSamples = 2.^[6 6 13];
c = 340;
M = 0.6;
uBL_u = 0.8;
surfaceStrength = -14.5;
nMakePlots = 0;     % 0: off, 1: plot, 2: plot and save, 3: Combo Only

% Frequency Space
freq.x = (-0.5:1/nSamples(2):0.5-1/nSamples(2))*sampleRate(2);
freq.y = reshape((-0.5:1/nSamples(1):0.5-1/nSamples(1))*sampleRate(1),nSamples(1),1,1);
freq.t = reshape((-0.5:1/nSamples(3):0.5-1/nSamples(3))*sampleRate(3),1,1,nSamples(3));
freq.rho = sqrt(freq.x.^2+freq.y.^2);
freq.theta = atan2(freq.y,freq.x);

%%%%% Aero-Optics Signal
AO.ellipsoid = [8 90 175];
AO.strength = -14.5;
AO.slope = -0.13;
% Calculations
AO.speed = c*M*uBL_u;
b = 1/2/AO.strength*(AO.strength^2-1/AO.slope^2);
AO.WF = zeros(nSamples);
XR = freq.x*cos(atan(-1/AO.speed))+freq.t*sin(atan(-1/AO.speed));
YR = freq.y;
TR = -freq.t*sin(atan(-1/AO.speed))+freq.x*cos(atan(-1/AO.speed));
R = sqrt((XR/AO.ellipsoid(1)).^2+(YR/AO.ellipsoid(2)).^2+(TR/AO.ellipsoid(3)).^2);
AO.WF = 10.^(b-sqrt(R.^2/AO.slope^2+b^2));
clear b XR YR TR R;

%%%%% White-Noise Stationary Signal
SN.rho0 = 5;
SN.strength = -14.5;
SN.slope = -0.175;
% Calculations
b = 1/2/SN.strength*(SN.strength^2-1/SN.slope^2);
SN.WF = 10.^(b-sqrt(repmat((freq.rho./(SN.rho0.*sqrt(10-6*cos(2*freq.theta+pi)))).^2,1,1,nSamples(3))/SN.slope^2+b^2));
clear b;

%%%%% Blade Pass Frequency Contamination
BPF.freq = 500;
BPF.harmonic = 0.5:0.5:5;
BPF.rho0 = 20;
BPF.tThickness = 100;
BPF.strength = -14;
BPF.slope = -0.13;
BPF.cutoff = 500;
BPF.aspectRatio = 1;
% Calculations
b = 1/2/BPF.strength*(BPF.strength^2-1/BPF.slope^2);
BPF.WF = zeros(nSamples);
for aa=1:length(BPF.harmonic)
    R = sqrt((sqrt(freq.x.^2+(BPF.aspectRatio*freq.y).^2)./(BPF.rho0/sqrt(1+((BPF.freq*BPF.harmonic(aa)-BPF.freq)/BPF.cutoff).^2)*sqrt(10-6*cos(2*freq.theta+pi)))).^2+((freq.t-BPF.freq*BPF.harmonic(aa))/BPF.tThickness).^2);
    BPF.WF = BPF.WF+10.^(b-sqrt(R.^2/BPF.slope^2+b^2));
    R = sqrt((sqrt(freq.x.^2+(BPF.aspectRatio*freq.y).^2)./(BPF.rho0/sqrt(1+((BPF.freq*BPF.harmonic(aa)-BPF.freq)/BPF.cutoff).^2)*sqrt(10-6*cos(2*freq.theta+pi)))).^2+((freq.t+BPF.freq*BPF.harmonic(aa))/BPF.tThickness).^2);
    BPF.WF = BPF.WF+10.^(b-sqrt(R.^2/BPF.slope^2+b^2));
end
clear b R;

%%%%% Zero Frequency Contamination
ZERO.rho0 = 25;
ZERO.tThickness = 50;
ZERO.strength = -14.5;
ZERO.slope = -0.5;
ZERO.aspectRatio = 0.55;
% Calculations
b = 1/2/ZERO.strength*(ZERO.strength^2-1/ZERO.slope^2);
R = sqrt((sqrt(freq.x.^2+(ZERO.aspectRatio*freq.y).^2)./(ZERO.rho0*sqrt(10-6*cos(2*freq.theta+pi)))).^2+(freq.t/ZERO.tThickness).^2);
ZERO.WF = 10.^(b-sqrt(R.^2/ZERO.slope^2+b^2));
clear b R;

%%%%% Acoustic Cone Signal
CONE.strength = [-13 -16];
CONE.slope = -0.3;
CONE.thickness = 8;
CONE.lowPassRho = 200;
CONE.lowPassX = 115;
% Calculations
freqMod.x0 = sin(0.5*atan(1/c/(M+1))+0.5*atan(1/c/(M-1)))*freq.t;
freqMod.y0 = 0;
freqMod.ax = sin(0.5*atan(1/c/(M+1))-0.5*atan(1/c/(M-1)))*freq.t;
freqMod.ay = sin(atan(1/c))*freq.t;
freqMod.theta = atan2(freq.y-freqMod.y0,freq.x-freqMod.x0);
freqMod.rho = (sqrt((freq.x-freqMod.x0).^2+(freq.y-freqMod.y0).^2)./sqrt((freqMod.ax.*cos(freqMod.theta)).^2+(freqMod.ay.*sin(freqMod.theta)).^2)-1).*sqrt((freqMod.ax.*cos(freqMod.theta)).^2+(freqMod.ay.*sin(freqMod.theta)).^2)/CONE.thickness;
freqMod.rho(:,:,end/2+1) = freq.rho(:,:)/CONE.thickness;
b1 = ((CONE.strength(2)-CONE.strength(1))/(sampleRate(3)/2)*abs(freq.t)+CONE.strength(1));
b = 1/2./b1.*(b1.^2-1/CONE.slope^2);
CONE.WF = 10.^(b-sqrt(freqMod.rho.^2/CONE.slope^2+b.^2));
CONE.WF = CONE.WF.*sqrt(1./(1+(freqMod.rho/CONE.lowPassRho).^2));
CONE.WF = CONE.WF.*sqrt(1./(1+(freq.x/CONE.lowPassX).^2));
clear freqMod b b1;

%%%%% Background Noise
BACK.strength = -18;
BACK.deviation = 0.75;
% Calculations
BACK.WF = 10.^(randn(nSamples)*BACK.deviation+BACK.strength);
BACK.WF(2:nSamples(1),2:nSamples(2),nSamples(3)/2+2:nSamples(3)) = flip(flip(flip(BACK.WF(2:nSamples(1),2:nSamples(2),2:nSamples(3)/2),1),2),3);

%%%%% Sound and Vibration
SV.WF = BPF.WF+ZERO.WF+CONE.WF;

%%%%% Plot
views = [-125 25; -55 25; 180 0; 270 0];
f1 = figure(1);
for aa=1:4
    subplot(2,2,aa)
    patch(isosurface(freq.x,freq.y,freq.t,AO.WF,10^surfaceStrength),'edgecolor','none','facecolor','red','facelighting','gouraud');
    patch(isosurface(freq.x,freq.y,freq.t,SN.WF,10^surfaceStrength),'edgecolor','none','facecolor','blue','facelighting','gouraud');
    patch(isosurface(freq.x,freq.y,freq.t,BPF.WF,10^surfaceStrength),'edgecolor','none','facecolor','green','facelighting','gouraud');
    patch(isosurface(freq.x,freq.y,freq.t,ZERO.WF,10^surfaceStrength),'edgecolor','none','facecolor','yellow','facelighting','gouraud');
    patch(isosurface(freq.x,freq.y,freq.t,CONE.WF,10^surfaceStrength),'edgecolor','none','facecolor','magenta','facelighting','gouraud');
    patch(isosurface(freq.x,freq.y,freq.t,BACK.WF,10^surfaceStrength),'edgecolor','none','facecolor','cyan','facelighting','gouraud');
    grid on;
    hold on;
    daspect([1 1 sampleRate(3)/sampleRate(1)/3]);
    xlim(sampleRate(1)/2*[-1 1]);
    ylim(sampleRate(2)/2*[-1 1]);
    zlim(sampleRate(3)/2*[0 1]);
    camlight;
    xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
    ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
    zlabel('$f\ (Hz)$','Interpreter','Latex');
%     title('Synthetic Signal','Interpreter','Latex');
    f1.Children(1).TickLabelInterpreter = 'latex';
    view(views(aa,:));
    camlight;
end
f1.Units = 'inches';
f1.Position = [1 1 6 8];
% Save Plot
saveas(f1,'synthetic_wavefront.eps','epsc');

%%%%% Animation
nFrames = 150;
theta = (0:nFrames-1)/(nFrames)*4*pi;
az = rad2deg(theta);
el = 25*sin(theta/2);

f2 = figure(2);
scolor = parula(2);
patch(isosurface(freq.x,freq.y,freq.t(end/2+1:end),AO.WF(:,:,end/2+1:end),10^surfaceStrength),'edgecolor','none','facecolor','red','facelighting','gouraud');
patch(isocaps(freq.x,freq.y,freq.t(end/2+1:end),AO.WF(:,:,end/2+1:end),10^surfaceStrength,'all'),'edgecolor','none','facecolor','red','facelighting','gouraud');
patch(isosurface(freq.x,freq.y,freq.t(end/2+1:end),SN.WF(:,:,end/2+1:end),10^surfaceStrength),'edgecolor','none','facecolor','blue','facelighting','gouraud');
patch(isocaps(freq.x,freq.y,freq.t(end/2+1:end),SN.WF(:,:,end/2+1:end),10^surfaceStrength,'all'),'edgecolor','none','facecolor','blue','facelighting','gouraud');
patch(isosurface(freq.x,freq.y,freq.t(end/2+1:end),BPF.WF(:,:,end/2+1:end),10^surfaceStrength),'edgecolor','none','facecolor','green','facelighting','gouraud');
patch(isocaps(freq.x,freq.y,freq.t(end/2+1:end),BPF.WF(:,:,end/2+1:end),10^surfaceStrength,'all'),'edgecolor','none','facecolor','green','facelighting','gouraud');
patch(isosurface(freq.x,freq.y,freq.t(end/2+1:end),ZERO.WF(:,:,end/2+1:end),10^surfaceStrength),'edgecolor','none','facecolor','yellow','facelighting','gouraud');
patch(isocaps(freq.x,freq.y,freq.t(end/2+1:end),ZERO.WF(:,:,end/2+1:end),10^surfaceStrength,'all'),'edgecolor','none','facecolor','yellow','facelighting','gouraud');
patch(isosurface(freq.x,freq.y,freq.t(end/2+1:end),CONE.WF(:,:,end/2+1:end),10^surfaceStrength),'edgecolor','none','facecolor','magenta','facelighting','gouraud');
patch(isocaps(freq.x,freq.y,freq.t(end/2+1:end),CONE.WF(:,:,end/2+1:end),10^surfaceStrength,'all'),'edgecolor','none','facecolor','magenta','facelighting','gouraud');
patch(isosurface(freq.x,freq.y,freq.t(end/2+1:end),BACK.WF(:,:,end/2+1:end),10^surfaceStrength),'edgecolor','none','facecolor','cyan','facelighting','gouraud');
patch(isocaps(freq.x,freq.y,freq.t(end/2+1:end),BACK.WF(:,:,end/2+1:end),10^surfaceStrength,'all'),'edgecolor','none','facecolor','cyan','facelighting','gouraud');
grid on;
daspect([1 1 50]);
xlim(sampleRate(1)/2*[-1 1]);
ylim(sampleRate(2)/2*[-1 1]);
zlim(sampleRate(3)/2*[0 1]);
xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
zlabel('$f\ (Hz)$','Interpreter','Latex');
f2.Children(1).TickLabelInterpreter = 'latex';
f2.Units = 'inches';
f2.Position = [1 1 5.5 6.25];
cl = camlight;

filename = 'synthetic_wavefront.gif';
frameRate = 15;
for aa=1:nFrames-1
    view(az(aa),el(aa));
    camlight(cl);
    drawnow;
    frame = getframe(f2);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if aa==1
        imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',1/frameRate);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/frameRate);
    end
end

%%%%% Make Wavefronts
[wf.x,wf.y] = meshgrid(0.975*(-1:2/(nSamples(1)-1):1),0.975*(-1:2/(nSamples(2)-1):1));
wf.rho = sqrt(wf.x.^2+wf.y.^2);
wf.theta = atan2(wf.y,wf.x);
wf.mask = ones(size(wf.x));
wf.mask(wf.rho>1) = NaN;
wf.x = wf.x.*wf.mask;
wf.y = wf.y.*wf.mask;
wf.rho = wf.rho.*wf.mask;
wf.theta = wf.theta.*wf.mask;
wf.sampleRate = sampleRate;
% Aero-Optics Signal
phase = pi*(2*rand(nSamples)-1);
phase(nSamples(1)/2+1,nSamples(2)/2+1,nSamples(3)/2+1) = 0;
phase(2:nSamples(1),2:nSamples(2),nSamples(3)/2+2:nSamples(3)) = -flip(flip(flip(phase(2:nSamples(1),2:nSamples(2),2:nSamples(3)/2),1),2),3);
wf.AO = real(ifftn(ifftshift(sqrt((AO.WF)*prod(sampleRate)*numel(AO.WF)).*exp(1i*phase))));
wf.AO = wf.AO.*wf.mask;
clear phase;
% White-Noise Stationary Signal
phase = pi*(2*rand(nSamples)-1);
phase(nSamples(1)/2+1,nSamples(2)/2+1,nSamples(3)/2+1) = 0;
phase(2:nSamples(1),2:nSamples(2),nSamples(3)/2+2:nSamples(3)) = -flip(flip(flip(phase(2:nSamples(1),2:nSamples(2),2:nSamples(3)/2),1),2),3);
wft.SN = real(ifftn(ifftshift(sqrt((SN.WF)*prod(sampleRate)*numel(SN.WF)).*exp(1i*phase))));
wft.SN = wft.SN.*wf.mask;
clear phase;
% Background Noise
phase = pi*(2*rand(nSamples)-1);
phase(nSamples(1)/2+1,nSamples(2)/2+1,nSamples(3)/2+1) = 0;
phase(2:nSamples(1),2:nSamples(2),nSamples(3)/2+2:nSamples(3)) = -flip(flip(flip(phase(2:nSamples(1),2:nSamples(2),2:nSamples(3)/2),1),2),3);
wft.BACK = real(ifftn(ifftshift(sqrt((BACK.WF)*prod(sampleRate)*numel(BACK.WF)).*exp(1i*phase))));
wft.BACK = wft.BACK.*wf.mask;
clear phase;
% Sound and Vibration
phase = pi*(2*rand(nSamples)-1);
phase(nSamples(1)/2+1,nSamples(2)/2+1,nSamples(3)/2+1) = 0;
phase(2:nSamples(1),2:nSamples(2),nSamples(3)/2+2:nSamples(3)) = -flip(flip(flip(phase(2:nSamples(1),2:nSamples(2),2:nSamples(3)/2),1),2),3);
wft.SV = real(ifftn(ifftshift(sqrt((SV.WF)*prod(sampleRate)*numel(SV.WF)).*exp(1i*phase))));
wft.SV = wft.SV.*wf.mask;
clear phase;
% Total
wf.wf = wf.AO+wft.SN+wft.BACK+wft.SV;
% Save
% save('synthetic_wavefront.mat','wf');
% disp('File Saved');