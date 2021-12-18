close all; clc; clearvars;

lambda = (1:0.1:10)*1e-6;
k = 2*pi./lambda;
z = 1000;
l = 1;
I = (k*l/8/z).^2;

%%%%% Plot
f1 = figure;
plot(lambda*1e6,I/I(1),'linewidth',1.25,'color',linspecer(1));
grid on;
xlabel('$\lambda$ ($\mu m$)','interpreter','latex');
ylabel('$I_0(\lambda)/I_0(1 \mu m)$','interpreter','latex');
% title('Diffraction Limited Far-Field Performance','interpreter','latex');
xlim([0 10]);
f1.Children(1).Layer = 'top';
f1.Children(1).TickLabelInterpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 5 3];
saveas(f1,'farfield_intensity.eps','epsc')









