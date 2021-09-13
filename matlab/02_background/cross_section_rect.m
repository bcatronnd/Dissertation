close all; clc; clearvars;

m = 0:2;
n = 0:3;
l_x = 1;
l_y = 1;
[x,y] = meshgrid((0:100)/100*l_x,(0:100)/100*l_y);

%%%%% Plot
f1 = figure(1);
sgtitle({'Characteristic Functions for a Rectangular Duct'},'interpreter','latex')
colormap(redblue)
for aa=1:length(m)
    for bb=1:length(n)
        subplot(length(m),length(n),(aa-1)*length(n)+bb);
        surf(x,y,cos(m(aa)*pi*x/l_x).*cos(n(bb)*pi*y/l_y),'linestyle','none'); view(2); hold on;
        k = sqrt((m(aa)*pi/l_x)^2+(n(bb)*pi/l_y)^2);
        axis equal;
        plot3([0 1 1 0 0],[0 0 1 1 0],ones(1,5),'-k','linewidth',0.1);
        if m(aa)==0 && n(bb)==0
            title({['$m=$' num2str(m(aa)) '  $n=$' num2str(n(bb))],'$k_m=0$'},'interpreter','latex');
        else
            contour(x,y,cos(m(aa)*pi*x/l_x).*cos(n(bb)*pi*y/l_y),[0 0],'-.k','linewidth',0.25);
            title({['$m=$' num2str(m(aa)) '  $n=$' num2str(n(bb))],['$k_m=$' num2str(k,'%.3f')]},'interpreter','latex');
        end
        xticks([0 1]); xticklabels({'$0$','$l_x$'});
        yticks([0 1]); yticklabels({'$0$','$l_y$'});
        f1.Children(1).TickLabelInterpreter = 'latex';
    end
end
f1.Units = 'inches';
f1.Position = [1 1 5 4];
saveas(f1,'cross_section_rect.eps','epsc')