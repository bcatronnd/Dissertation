close all; clc; clearvars;

m = 0:2;
n = 0:3;
[r,t] = meshgrid((0:100)/100,(0:100)/100*2*pi);
x = r.*cos(t);
y = r.*sin(t);

%%%%% Plot
f1 = figure(1);
sgtitle({'Characteristic Functions for a Circular Duct'},'interpreter','latex');
colormap(redblue);
for aa=1:length(m)
    for bb=1:length(n)
        subplot(length(m),length(n),(aa-1)*length(n)+bb);
        k = max(besseljdzero(m(aa),n(bb)+1,0==m(aa)));
%         if m(aa)==0 && n(bb)==0
%             k = 0;
%         elseif m(aa)==0
%             k = BessDerivZerosBisect2(0,n(bb));
%         else
%             k = BessDerivZerosBisect2(m(aa),n(bb)+1);
%         end
        surf(x,y,besselj(m(aa),k*r).*real(exp(1i*(m(aa)*t))),'linestyle','none'); view(2); axis off equal; hold on;
        plot3(cos((0:100)/100*2*pi),sin((0:100)/100*2*pi),1*ones(1,101),'-k','linewidth',0.1);
        if m(aa)==0 && n(bb)==0
            title({['$m=\ $' num2str(m(aa)) '$\quad n=\ $' num2str(n(bb))],'$k_{m}=\ 0$'},'interpreter','latex');
        else
            contour(x,y,besselj(m(aa),k*r).*real(exp(1i*(m(aa)*t))),[0 0],'-.k','linewidth',0.25);
            caxis([-1 1]);
            title({['$m=\ $' num2str(m(aa)) '$\quad n=\ $' num2str(n(bb))],['$k_{m}=\ $' num2str(k,'%.3f')]},'interpreter','latex');
        end
    end
end
f1.Units = 'inches';
f1.Position = [1 1 5 4];

saveas(f1,'cross_section_circ.eps','epsc')