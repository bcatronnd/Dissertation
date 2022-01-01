close all; clc; clearvars;

% Inputs
M = 0.5;
Kgd = 2.27e-4;
gamma = 90;
delta = 0.022;
rho = 1.04;
cf = 2.25e-3;

% Code
G = 1-0.19*M^2+0.03*M^4;
B = 0.19/sind(gamma);

opdrms = B*Kgd*rho*M^2*delta*sqrt(cf)*G*1e6


