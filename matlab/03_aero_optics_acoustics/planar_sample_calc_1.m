clearvars; close all; clc;

% Lp = 20*log10(1/20e-6);
% l_n = 1;
% D = 1;
% LAMBDA = logspace(-1,1,250);
THETA = [0 10 20];
% 
% f = waitbar(0,'Please wait...');
% for aa=1:length(LAMBDA)
%     for bb = 1:length(THETA)
%         [~,results(aa,bb)] = PotentialWaveOPD(D,LAMBDA(aa),THETA(bb),l_n,Lp);
%     end
%     waitbar(aa/length(LAMBDA),f)
% end
% clear aa bb
% close(f)
% close all
% 
% save([mfilename '.mat'],'LAMBDA','D','results')
load([mfilename '.mat'])

fig = figure;
axes('ColorOrder',brewermap(10,'Set1'),'NextPlot','replacechildren')
semilogx(LAMBDA/D,results*1e6,'linewidth',1.5)
grid on
xlabel('$\Lambda/Ap$')
ylabel('$OPD_{RMS}$ ($\mu m$)')
for aa=1:length(THETA)
    s_legend{aa} = ['$\theta=' num2str(THETA(aa)) '^\circ$'];
end
legend(s_legend);

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
fig.CurrentAxes.Title.Interpreter = 'latex';
fig.CurrentAxes.Title.FontSize = 12;
saveas(fig,'planar_sample_calc_1.eps','epsc')

% fig.Position = [0 0 800 600];
% saveas(fig,[mfilename '.png'],'png')
% pause(10)
% close all
