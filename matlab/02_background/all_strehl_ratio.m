close all; clc; clearvars;

lambda = (1:0.1:10)*1e-6;
OPD = 3.821e-7;
SR = exp(-1*(2*pi*OPD./lambda).^2);
PlotColors = linspecer(2);

%%%%% Plot
f1 = figure(1);
plot(lambda*1e6,SR,'linewidth',1.25,'color',PlotColors(2,:)); grid on;
xlabel('$\lambda$ ($\mu m$)','interpreter','latex');
ylabel('Strehl Ratio, SR ($I/I_0$)','interpreter','latex');
xlim([0 10]);
f1.Children(1).Layer = 'top';
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 3 3];

saveas(f1,'all_strehl_ratio.eps','epsc')