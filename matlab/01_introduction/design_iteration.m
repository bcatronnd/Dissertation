close all; clc; clearvars;

nIterations = 10;

x = 1:nIterations;
yTunnel = 4*ones(1,nIterations);
% yModel = fliplr(1:nIterations).^1.05+(2*rand(1,nIterations)-1);
yModel = [11.5627132839231,9.92239853147905,9.54355696772071,8.25301071267418,5.89684837734810,6.14295289107599,5.26683815740818,3.19824883866002,2.83909189393667,1.17605211061700];
yReq = 6*ones(1,nIterations);
PlotColors = linspecer(4);

f1 = figure(1);
plot(x,yTunnel,'linewidth',1.25,'color',PlotColors(4,:));
hold on;
plot(x,yReq,'linewidth',1.25,'color',PlotColors(2,:));
plot(x,sqrt(yModel.^2+yTunnel.^2),'linewidth',1.25,'marker','o','color',PlotColors(1,:),'markerfacecolor',PlotColors(1,:));
plot(x,yModel,'linewidth',1.25,'marker','d','color',PlotColors(3,:),'markerfacecolor',PlotColors(3,:));

xlim([1 nIterations]);
grid on;
xticks(x);
yticks([]);
xlabel('Design Iteration','interpreter','latex');
ylabel('OPD$_{RMS}$','interpreter','latex');
legend('Environmental Contamination','Design Requirements','Measured Performance','True Performance','interpreter','latex');
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 5 3];

saveas(f1,'design_iteration.eps','epsc');