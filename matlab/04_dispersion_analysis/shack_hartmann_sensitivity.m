close all; clc; clearvars;

xPoints = 100;
pPoints = 50;
% lambda = logspace(-1,1,250);
step = -7;
lambda = 1./(2^step:2^step:10);

phi = (0:pPoints-1)'/(pPoints-1)*2*pi;
x = (0:xPoints-1)/(xPoints-1);

for aa=1:length(lambda)
    yp = -2*pi/lambda(aa)*sin(x*2*pi/lambda(aa)+phi);
    ypmax(aa) = max(mean(yp,2));
    yprms(aa) = rms(mean(yp,2));
end
clear yp;

f1 = figure(1);
semilogx(lambda,yprms);
grid on;
xlabel('$\Lambda/Ap$','interpreter','latex');
ylabel('$\tan(\theta)_{RMS} \approx d_{RMS}/f$','interpreter','latex');
title('Shack-Hartmann Sensitivity','interpreter','latex');
% xlim([1e-1 1e1]);
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 5.5 4];

saveas(f1,'shack_hartmann_sensitivity.eps','epsc');