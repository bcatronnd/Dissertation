clearvars; close all; clc;

% Lp = linspace(80,150,100);
% l_n = 1;
% D = 1;
% LAMBDA = [0.1 1 10];
% THETA = 0;
% 
% f = waitbar(0,'Please wait...');
% for aa=1:length(Lp)
%     for bb = 1:length(LAMBDA)
%         [~,results(aa,bb)] = PotentialWaveOPD(D,LAMBDA(bb),THETA,l_n,Lp(aa));
%     end
%     waitbar(aa/length(Lp),f)
% end
% close(f)
% clear aa bb
% 
% save([mfilename '.mat'],'Lp','LAMBDA','results')
load([mfilename '.mat'])

fig = figure;
axes('ColorOrder',brewermap(10,'Set1'),'NextPlot','replacechildren')
semilogy(Lp,results*1e6,'linewidth',1.5)
grid on
xlabel('SPL (dB)')
ylabel('$\overline {OPD_{rms}}/m$ ($\mu m/m$)')
for aa=1:length(LAMBDA)
    s_legend{aa} = ['$\Lambda/Ap=' num2str(LAMBDA(aa)) '$'];
end
legend(s_legend)
ylim([1e-6 1e0])

%%%%% Print Formatting
% figProps = struct(fig);
fig.Units = 'inches';
fig.Position = [1 1 5 2.5];
fig.CurrentAxes.XLabel.Interpreter = 'latex';
fig.CurrentAxes.XLabel.FontSize = 12;
fig.CurrentAxes.YLabel.Interpreter = 'latex';
fig.CurrentAxes.YLabel.FontSize = 12;
fig.CurrentAxes.TickLabelInterpreter = 'latex';
fig.CurrentAxes.XAxis.FontSize = 10;
fig.CurrentAxes.YAxis.FontSize = 10;
fig.CurrentAxes.Legend.Interpreter = 'latex';
fig.CurrentAxes.Legend.FontSize = 10;
fig.CurrentAxes.Legend.Location = 'northwest';
fig.CurrentAxes.YMinorGrid = 'off';
saveas(fig,'planar_sample_calc_3.eps','epsc')