close all; clc; clearvars;

nLenslets = 32;
lambda = 0.25:0.025:5;
dist = 5:5:25;
zMax = 5;
dz = 0.05;
c0 = 340;
kgd = 2.27e-4;
phaseSteps = 25;
AP = 0.25;

dist = dist*AP;
lambda = lambda*AP;

[x,y,z] = meshgrid(0.975*AP*(-0.5:1/(nLenslets-1):0.5),0.975*AP*(-0.5:1/(nLenslets-1):0.5),-zMax:dz:zMax);
r = sqrt(x(:,:,1).^2+y(:,:,1).^2);
t = atan2(y(:,:,1),x(:,:,1));
mask = ones(nLenslets);
mask(r>0.5*AP) = NaN;
phase = (0:phaseSteps-1)/phaseSteps*2*pi;


OPDrms = zeros(length(lambda),length(dist));
OPDrmsMax = OPDrms;
OPDrmsMin = OPDrms;

for aa=1:length(lambda)
    k = 2*pi/lambda(aa);
    for bb=1:length(dist)
        R = sqrt((x+dist(bb)).^2+y.^2+z.^2);
        OPD = kgd/c0^2*sum(1./R.*exp(-1i*k*R),3)*dz;
        stats = zeros(1,phaseSteps);
        for cc=1:phaseSteps
            wf = mask.*real(OPD.*exp(1i*phase(cc)));
            wf = ZernikeRemoval(1:3,wf,r.*mask,t.*mask,'noll');
            stats(cc) = nanrms(wf(:));
        end
        OPDrms(aa,bb) = mean(stats);
        OPDrmsMax(aa,bb) = max(stats);
        OPDrmsMin(aa,bb) = min(stats);
    end
end




%%
close all;

a = 1e-3;
b = -0.5;

f1 = figure(1);
plot(lambda/AP,OPDrms*1e6.*sqrt(dist/AP));
grid on;
xlabel('$\Lambda/Ap$','interpreter','latex');
ylabel('$\frac{OPD_{RMS}\sqrt{R/Ap}}{A_O}\ (\frac{\mu m}{kg/s^2})$','interpreter','latex','fontsize',14);
for aa=1:length(dist)
    sLegend{aa} = ['R/Ap = ' num2str(dist(aa)/AP)];
end
legend(sLegend,'interpreter','latex');
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 5.5 3.5];

saveas(f1,'spherical_sample.eps','epsc');

X = lambda/AP;
Y = mean(OPDrms*1e6.*sqrt(dist/AP),2);
% f(x) = (p1*x^3 + p2*x^2 + p3*x + p4) / 
%        (x^2 + q1*x + q2)
% p1 =  -2.205e-05  (-2.763e-05, -1.647e-05)
% p2 =   0.0001551  (0.0001131, 0.000197)
% p3 =   0.0002423  (0.0001491, 0.0003355)
% p4 =   0.0003664  (0.0003376, 0.0003952)
% q1 =      -1.161  (-1.206, -1.116)
% q2 =      0.8607  (0.8365, 0.885)
% Goodness of fit:
%   SSE: 1.475e-08
%   R-square: 0.9993
%   Adjusted R-square: 0.9992
%   RMSE: 8.93e-06