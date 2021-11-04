close all; clc; clearvars;

nIterations = 10;

x = 1:nIterations;
yTunnel = 4*ones(1,nIterations);
% yModel = fliplr(1:nIterations).^1.05+(2*rand(1,nIterations)-1);
yModel = [11.5627132839231,9.92239853147905,9.54355696772071,8.25301071267418,5.89684837734810,6.14295289107599,5.26683815740818,3.19824883866002,2.83909189393667,1.17605211061700];
yReq = 6*ones(1,nIterations);


f1 = figure(1);
plot(x,sqrt(yModel.^2+yTunnel.^2),'-ok',x,yModel,':ok',x,yTunnel,'b-',x,yReq,'r-');
grid on;
yticks([]);
xlabel('Design Iteration','interpreter','latex');
ylabel('OPD$_{RMS}$','interpreter','latex');
legend('Measured Performance','True Performance','Environmental Contamination','Design Requirements','interpreter','latex');
f1.Children(2).TickLabelInterpreter = 'latex';

saveas(f1,'design_iteration.eps','epsc');