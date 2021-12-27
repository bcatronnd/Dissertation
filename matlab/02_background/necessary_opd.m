close all; clc; clearvars;

lambda = (1:0.1:10)*1e-6;
SR = [0.9 0.7 0.5];
[lam,sr] = meshgrid(lambda,SR);
OPD = lam/2/pi.*sqrt(-log(sr));
PlotColors = linspecer(3);

%%%%% Plot
f1 = figure(1);
for aa=1:length(SR)
    plot(lambda*1e6,OPD(aa,:)*1e6,'linewidth',1.25,'color',PlotColors(aa,:));
    hold on;
end
grid on;
xlim([1 10]);
xlabel('$\lambda$ ($\mu m$)','interpreter','latex');
ylabel('$\textrm{OPD}_{RMS}$ ($\mu m$)','interpreter','latex');
leg_sr = [repmat('SR = ',3,1) num2str(SR')];
legend(leg_sr,'location','northwest')
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Children(1).Interpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 5 3];

saveas(f1,'necessary_opd.eps','epsc')