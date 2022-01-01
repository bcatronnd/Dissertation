close all; clc; clearvars;

load('../07_multiple_sensor_filtering/lse_mspod_table.mat');
M = [0.3 0.4 0.5];

SensorCount = zeros(1,length(SensorSelection));
for aa=1:length(SensorSelection)
    SensorCount(aa) = length(SensorSelection{aa});
end


close(findobj('type','figure','number',1));
f1 = figure(1);
colors = linspecer(5);

sLegend = cell(1,3);
for aa=1:length(M)
    plot(SensorCount,(1-DataSummary(3:7,aa)./DataSummary(2,aa))*100,'o','markerfacecolor',colors(aa,:),'markeredgecolor',colors(aa,:));
    hold on;
    sLegend{aa} = ['M=' num2str(M(aa),'%0.1f')];
end
% for aa=1:length(M)
%     p = polyfit(SensorCount,(1-DataSummary(3:7,aa)./DataSummary(2,aa))*100,1);
%     plot([min(SensorCount) max(SensorCount)],polyval(p,[min(SensorCount) max(SensorCount)]),'linewidth',1.25,'color',colors(aa,:));
% end
grid on;
xlabel('Number of Additional Sensors','interpreter','latex');
ylabel('Reduction in $OPD_{RMS}$ (\%)','interpreter','latex');
legend(sLegend,'interpreter','latex','location','northwest');
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 5.5 3.25];

saveas(f1,'lse_summary.eps','epsc');


