close all; clc; clearvars;

% Input
M = 0.5;
delta = 0.022;
rho = 1.04;
gamma = 90;
Re = 21e-3;
Cf = 2.25e-3;
Ap = 5*0.0254;
Kgd = 2.27e-4;

% Fit Values
a = 0.7795;
b = -0.9616;
c = -0.1382;

% Calculation
B = 0.19/sind(gamma);
G = 1-0.19*M^2+0.03*M^4;
Ap_corr = 1-a*exp(b*Ap/delta)-(1-a)*exp(c*Ap/delta);
OPDrms = B*Kgd*rho*M^2*delta*sqrt(Cf)*G*Ap_corr;
OPDrms2 = sqrt(2*OPDrms^2);

disp(['Single BL OPDrms: ' num2str(OPDrms*1e6,'%0.4f') ' ' char(956) 'm'])
disp(['Double BL OPDrms: ' num2str(OPDrms2*1e6,'%0.4f') ' ' char(956) 'm'])
