close all; clc; clearvars;

lambda = (350:840)/1e3;
R = 287.058;
kgd = 0.776*R*(1+0.00753./lambda.^2)*1e-6;

%%%%% Plot
f1 = figure(1);
plot(lambda,kgd,'k','linewidth',1.25); 
xlim([360 830]/1e3); grid on; axis manual;
ylabel({'Gladstone-Dale Constant';'$K_{GD}$ (m$^3$/kg)'},'interpreter','latex');
spectrumLabel(gca);
xlabel('Wavelength of Light, $\lambda$ (nm)','interpreter','latex');
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Children(2).XTick = (400:50:800)/1e3;
f1.Children(2).XTickLabel = {};
f1.Units = 'inches';
f1.Position = [1 1 5 3];
f1.Children(2).Position = [0.1300 0.1100+0.05 0.7750 0.8150-0.05];
f1.Children(1).Position = [0.1300 0.0624+0.05 0.7750 0.0476];
saveas(f1,'gladstone_dale_wavelength.eps','epsc')