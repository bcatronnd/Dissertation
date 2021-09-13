clearvars; close all; clc;

% freq = 267*2;
% c = 340;
% M = 0.6;
% u = c*M;
% 
% Lp = 20*log10(1/20e-6);
% l_n = 1;
% D = 6*0.0254;
% LAMBDA = [(u+c) abs(u-c) u]/freq;
% THETA = 0:0.25:45;
% 
% f = waitbar(0,'Please wait...');
% for aa=1:length(THETA)
%     for bb = 1:length(LAMBDA)
%         [~,results(aa,bb)] = PotentialWaveOPD(D,LAMBDA(bb),THETA(aa),l_n,Lp);
%     end
%     waitbar(aa/length(THETA),f)
% end
% close(f)
% clear aa bb
% 
% save([mfilename '_u.mat'],'THETA','results')
load([mfilename '.mat'])

fig = figure;
axes('ColorOrder',brewermap(10,'Set1'),'NextPlot','replacechildren')
plot(THETA,results*1e6,'linewidth',1.5)
grid on
xlabel('$\theta$ (deg)')
ylabel('$\overline {OPD_{rms}}$ ($\mu m$)')
legend('u+c','u-c')
xlim([0 45])

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
fig.CurrentAxes.Legend.Location = 'northeast';
saveas(fig,'planar_sample_calc_2.eps','epsc')

