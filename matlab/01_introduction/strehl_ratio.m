close all; clc; clearvars;

lambda = (1:0.1:10)*1e-6;
OPD = 3.821e-7;
SR = exp(-1*(2*pi*OPD./lambda).^2);

%%%%% Plot
f1 = figure;
plot(lambda*1e6,SR,'linewidth',1.25,'color',linspecer(1));
grid on;
xlabel('$\lambda$ ($\mu m$)','interpreter','latex');
ylabel('Strehl Ratio, SR ($I/I_0$)','interpreter','latex');
xlim([0 10]);
f1.Children(1).Layer = 'top';
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 5 3];
saveas(f1,'strehl_ratio.eps','epsc')