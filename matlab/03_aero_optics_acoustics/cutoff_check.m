close all; clc; clearvars;

% Input
f = 536;
c = 350;
M = 0.6;
lx = 36*0.0254;
ly = 31*0.0254;
m = 1;
n = 3;


%
w = 2*pi*f;
k0 = w/c;
kx = m*pi/lx;
ky = n*pi/ly;
km = sqrt(kx^2+ky^2);
kzp = (-M*k0+sqrt(k0^2-(1-M^2)*km^2))/(1-M^2);
kzm = (+M*k0+sqrt(k0^2-(1-M^2)*km^2))/(1-M^2);
fcuton = c/2/pi*sqrt((1-M^2)*km^2);

z = 0:0.01:1;
p = exp(1i*kzm*z);

f1 = figure(1);
colororder(linspecer(5));
plot(z,abs(p),'linewidth',1.25);
grid on;
xlabel('z (m)','interpreter','latex');
ylabel('$|\hat{p}(z)|/|\hat{p}(0)|$','interpreter','latex'); 
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 5.5 3.5];

% saveas(f1,'cutoff_check.eps','epsc');
saveas(f1,'cutoff_check.png');
