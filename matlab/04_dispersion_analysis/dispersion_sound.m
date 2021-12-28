close all; clc; clearvars;

mach = [0 0.3 0.6];
c = 340;

m = 0:3;
n = 0:3;
fRange = 0:10:500;

% CODE
PlotColors = linspecer(length(mach));

f1 = figure(1);
subplot(1,2,1);
hold on;
for aa=1:length(mach)
    plot(NaN,NaN,'linewidth',1.25,'color',PlotColors(aa,:));
    slegend{aa} = ['M=' num2str(mach(aa),'%0.2f')];
end
for aa=1:length(mach)
    for bb=1:length(m)
        kx = m(bb)*pi;
        ky = 0;
        f = [fRange c/2/pi*sqrt((1-mach(aa)^2).*(kx^2+ky^2))]+1e-6;
        f = unique(f);
        k0 = 2*pi*f/c;
        kzp = (-mach(aa)*k0+sqrt(k0.^2-(1-mach(aa)^2)*(kx^2+ky^2)))./(1-mach(aa)^2);
        kzm = (+mach(aa)*k0+sqrt(k0.^2-(1-mach(aa)^2)*(kx^2+ky^2)))./(1-mach(aa)^2);
        kzp(imag(kzp)>0) = NaN;
        kzm(imag(kzm)>0) = NaN;
        
        plot3(kzp/2/pi,repmat(kx/2/pi,1,length(f)),f,'linewidth',1.25,'color',PlotColors(aa,:));
%         plot3(kzp/2/pi,repmat(-kx/2/pi,1,length(f)),f,'linewidth',1.25,'color',PlotColors(aa,:),'linestyle',':');
        plot3(-kzm/2/pi,repmat(kx/2/pi,1,length(f)),f,'linewidth',1.25,'color',PlotColors(aa,:));
%         plot3(-kzm/2/pi,repmat(-kx/2/pi,1,length(f)),f,'linewidth',1.25,'color',PlotColors(aa,:),'linestyle',':');
    end
end
grid on;
view(130,30);
% xlim(5*[-1 1]);
daspect([1 1 60]);
zlim([0 fRange(end)]);
xlabel('$\xi_x\ (m^{-1})$','interpreter','latex');
ylabel('$\xi_y\ (m^{-1})$','interpreter','latex');
zlabel('$f\ (Hz)$','interpreter','latex');
title(['$m=0,\ n\in [0,' num2str(m(end)) ']$'],'interpreter','latex');
f1.Children(1).TickLabelInterpreter = 'latex';
% legend(slegend,'interpreter','latex','location','southoutside','orientation','horizontal');

subplot(1,2,2);
hold on;
for aa=1:length(mach)
    plot(NaN,NaN,'linewidth',1.25,'color',PlotColors(aa,:));
    slegend{aa} = ['M=' num2str(mach(aa),'%0.1f')];
end
for aa=1:length(mach)
    for bb=1:length(n)
        kx = 0;
        ky = n(bb)*pi;
        f = [fRange c/2/pi*sqrt((1-mach(aa)^2).*(kx^2+ky^2))]+1e-6;
        f = unique(f);
        k0 = 2*pi*f/c;
        kzp = (-mach(aa)*k0+sqrt(k0.^2-(1-mach(aa)^2)*(kx^2+ky^2)))./(1-mach(aa)^2);
        kzm = (+mach(aa)*k0+sqrt(k0.^2-(1-mach(aa)^2)*(kx^2+ky^2)))./(1-mach(aa)^2);
        kzp(imag(kzp)>0) = NaN;
        kzm(imag(kzm)>0) = NaN;
        
        plot3(kzp/2/pi,repmat(kx/2/pi,1,length(f)),f,'linewidth',1.25,'color',PlotColors(aa,:));
        plot3(kzp/2/pi,repmat(-kx/2/pi,1,length(f)),f,'linewidth',1.25,'color',PlotColors(aa,:));
        plot3(-kzm/2/pi,repmat(kx/2/pi,1,length(f)),f,'linewidth',1.25,'color',PlotColors(aa,:));
        plot3(-kzm/2/pi,repmat(-kx/2/pi,1,length(f)),f,'linewidth',1.25,'color',PlotColors(aa,:));
    end
end
grid on;
view(180,0);
% xlim(5*[-1 1]);
zlim([0 fRange(end)]);
xlabel('$\xi_x\ (m^{-1})$','interpreter','latex');
ylabel('$\xi_y\ (m^{-1})$','interpreter','latex');
zlabel('$f\ (Hz)$','interpreter','latex');
title(['$m\in [0,' num2str(n(end)) '],\ n=0$'],'interpreter','latex');
f1.Children(1).TickLabelInterpreter = 'latex';
legend(slegend,'interpreter','latex','location','southeast');%'southoutside','orientation','horizontal');
f1.Children(2).Position(4) = f1.Children(2).Position(4)-0.025;
f1.Children(2).Position(2) = f1.Children(2).Position(2)+0.025;
f1.Children(3).Position(4) = f1.Children(3).Position(4)-0.025;
f1.Children(3).Position(2) = f1.Children(3).Position(2)+0.025;

f1.Units = 'inches';
f1.Position = [1 1 5.5 3];

% saveas(f1,'dispersion_sound.eps','epsc');










