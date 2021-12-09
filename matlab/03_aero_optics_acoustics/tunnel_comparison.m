close all; clc; clearvars;


% Measured
dir = '/media/briancatron/ResearchData/TUNNEL_NOISE/FULL_WAVEFRONTS/EMPTY_TUNNEL/';
test_point = 'Jan0418022';
pFrame = 1000;

filter_freq = 536;
filter_freq_half_width = 25;
sample_rate = 34000;

load([dir test_point '_WF.mat']);
load([dir test_point '_WFdata.mat'],'Mask_WF');
Mask_WF(Mask_WF==0) = NaN;
WF1 = 1e6*WF.*Mask_WF;
clear WF;
WF1 = WF1-mean(WF1,3);
% Permute WF
WF1 = permute(WF1,[2 1 3]);
Mask_WF = permute(Mask_WF,[2 1]);
% POD Filter
WF1 = podfilter(WF1,[],'sigma',3);
% Remove piston/tip/tilt
[R,T] = ZernikeRT(WF1,0.975);
WF1 = ZernikeRemoval(1:3,WF1,R,T,'noll');
WF1 = WF1.*Mask_WF;
% Apply Frequency Filters
WF1 = filterWF(WF1,(filter_freq+filter_freq_half_width*[-1 1])/sample_rate,'bandpass').*Mask_WF;
% wf2 = filterWF(WF,(filter_freq+filter_freq_half_width*[-1 1])/sample_rate,'stop').*Mask_WF;

%% Simulated Beam
% Fluid Specs
c = 340;
M = 0.6;
f = 536;
kgd = 2.27e-4;
% Duct Specs
lx = 36*0.0254;
ly = 36*0.0254;
dx = 0.1*0.0254;
% Aperture Specs
ap = 10*0.0254;
nap = 55;
gamma = 90;
yoff = 12.5*0.0254;
% Duct Modes
duct_modes.m = [0 1 2 2 3];
duct_modes.n = [2 3 0 2 1];
duct_modes.c = [0.0191 0.0145 0.0194 0.0363 0.0060];
duct_modes.p = [0 0.5 1 1 1.5]*pi;
duct_modes.d = [-1 -1 -1 -1 -1];

%%%%%
phase = (0:0.01:1)*2*pi;
u = c*M;
beta = 1-M^2;
k0 = 2*pi*f/c;
kx = duct_modes.m*pi/(lx);
ky = duct_modes.n*pi/(ly);
kz = (-duct_modes.d*M*k0+sqrt(k0^2-beta*(kx.^2+ky.^2)))/beta;
%%%%%
% Correct RMS=1 Modes
[x,y] = meshgrid(0:0.01:1);
for aa=1:length(duct_modes.m)
    phirms(aa) = rms(reshape(cos(duct_modes.m(aa)*pi*x).*cos(duct_modes.n(aa)*pi*y),1,[]));
end
duct_modes.c = duct_modes.c.*phirms;
%%%%% Aperture OPD
% Beam Coordinates
xb = (0:dx:lx)/sind(gamma);
yb = linspace(-ap/2,ap/2,nap)*0.975;
zb = linspace(-ap/2,ap/2,nap)*0.975;
[X1,Y1,Z1] = meshgrid(xb,yb,zb);
R1 = sqrt(Y1.^2+Z1.^2);
Z1 = Z1/sind(gamma);
% Aperture Coordinates
[Za,Ya] = meshgrid(zb,yb);
R = sqrt(Za.^2+Ya.^2)/(ap/2);
T = atan2(Ya,Za);
% Duct Coordinates
X0 = X1*sind(gamma);
Y0 = Y1+yoff;
Z0 = Z1-X1*cosd(gamma);
% NaN Outside Aperture
X0(R1>ap/2) = NaN;
Y0(R1>ap/2) = NaN;
Z0(R1>ap/2) = NaN;
X1(R1>ap/2) = NaN;
Y1(R1>ap/2) = NaN;
Z1(R1>ap/2) = NaN;
% Beam OPD
for aa=1:length(duct_modes.m)
    ncp(:,:,:,aa) = kgd/c^2*cos(kx(aa)*X0).*cos(ky(aa)*Y0).*exp(-duct_modes.d(aa)*1i*kz(aa)*Z0).*exp(1i*duct_modes.p(aa));
end
opd = real(squeeze(sum(repmat(squeeze(trapz(ncp,2)*dx/sind(gamma)),1,1,1,length(phase)).*repmat(reshape(exp(1i*phase),1,1,1,[]),size(ncp,1),size(ncp,3),size(ncp,4),1),3)));
opd = ZernikeRemoval(1:3,opd,R,T,'noll');
%%%%% Plot
% for aa=1:length(phase)-1
%     f1 = figure(1);
%     surf(opd(:,:,aa),'linestyle','none');
%     view(2);
%     colormap(redblue);
%     axis image tight off;
%     caxis(max(abs(opd),[],'all')*[-1 1]);
%     f1.Color = [1 1 1 1];
%     drawnow;
%     %%%%% Make Gif
%     frame = getframe(f1);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if aa==1
%         imwrite(imind,cm,'simulated_beam.gif','gif','LoopCount',Inf,'DelayTime',0.1);
%     else
%         imwrite(imind,cm,'simulated_beam.gif','gif','WriteMode','append','DelayTime',0.1);
%     end
%
%
% end

mask_opd = [NaN*ones(55,5) Mask_WF];


%%
pFrame = 1300;
pPhase = 36;

f1 = figure(1);
colormap(redblue);
subplot(1,2,1);
surf(WF1(:,:,pFrame),'linestyle','none');
view(2);
axis image tight off;
% colorbar;
caxis(max(abs(WF1),[],'all')*[-1 1]);
title('Measured Wavefront','interpreter','latex');

subplot(1,2,2);
surf(mask_opd.*fliplr(opd(:,:,pPhase)),'linestyle','none');
view(2);
axis image tight off;
caxis(max(abs(opd),[],'all')*[-1 1]);
title('Simulated Wavefront','interpreter','latex');

f1.Units = 'inches';
f1.Position = [1 1 6 3];

saveas(f1,'tunnel_comparison.eps','epsc');







