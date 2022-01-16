close all; clc; clearvars;

mach = 0.6;
c = 340;
fRange = 0:5:500;
lx = 1;
ly = 1;
gamma = 75:15:120;

[m,n] = meshgrid(0:5);
m = m(:);
n = n(:);
kx = m*pi/lx;
ky = n*pi/ly;
km = sqrt(kx.^2+ky.^2);
validModes = (c/2/pi*sqrt((1-mach^2)*km.^2))<max(fRange);
m(validModes==0) = [];
n(validModes==0) = [];
kx = m*pi/lx;
ky = n*pi/ly;
km = sqrt(kx.^2+ky.^2);

fRange = [fRange reshape(c/2/pi*sqrt((1-mach^2).*(kx.^2+ky.^2)),1,[])]+1e-12;
fRange = unique(fRange);
k0 = 2*pi*fRange/c;
kzp = (-mach*k0+sqrt(k0.^2-(1-mach^2)*(kx.^2+ky.^2)))./(1-mach^2);
kzm = (+mach*k0+sqrt(k0.^2-(1-mach^2)*(kx.^2+ky.^2)))./(1-mach^2);
kzp(imag(kzp)>0) = NaN;
kzm(imag(kzm)>0) = NaN;

PlotColors = linspecer(max(m)+1);
PlotLines = {'-','--','-.',':'};


%%
close(findobj('type','figure','number',1));
f1 = figure(1);

for aa=1:length(gamma)
    khp = kzp*sind(gamma(aa))-kx*cosd(gamma(aa));
    khm = -kzm*sind(gamma(aa))+kx*cosd(gamma(aa));
    
    subplot(2,2,aa);
    plot3(NaN,NaN,NaN,'color',PlotColors(1,:),'linewidth',1.25);
    hold on;
    plot3(NaN,NaN,NaN,'color',PlotColors(2,:),'linewidth',1.25);
    plot3(NaN,NaN,NaN,'color',PlotColors(3,:),'linewidth',1.25);
    plot3(NaN,NaN,NaN,'color',PlotColors(4,:),'linewidth',1.25);
    plot3(NaN,NaN,NaN,'color','k','linestyle',PlotLines{1},'linewidth',1.25);
    plot3(NaN,NaN,NaN,'color','k','linestyle',PlotLines{2},'linewidth',1.25);
    plot3(NaN,NaN,NaN,'color','k','linestyle',PlotLines{3},'linewidth',1.25);
    plot3(NaN,NaN,NaN,'color','k','linestyle',PlotLines{4},'linewidth',1.25);
    for bb=1:length(m)
        plot3(khp(bb,:)/2/pi,repmat(ky(bb)/2/pi,1,length(fRange)),fRange,khp(bb,:)/2/pi,repmat(-ky(bb)/2/pi,1,length(fRange)),fRange,khm(bb,:)/2/pi,repmat(ky(bb)/2/pi,1,length(fRange)),fRange,khm(bb,:)/2/pi,repmat(-ky(bb)/2/pi,1,length(fRange)),fRange,'color',PlotColors(m(bb)+1,:),'linestyle',PlotLines{n(bb)+1},'linewidth',1.25);
    end
    xlim([-4 2]);
    ylim([-2 2]);
    zlim([0 500]);
    
    grid on;
    view(-15,15);
    daspect([1 1 60]);
    xlabel('$\xi_x\ (m^{-1})$','interpreter','latex');
    ylabel('$\xi_y\ (m^{-1})$','interpreter','latex');
    zlabel('$f\ (Hz)$','interpreter','latex');
    title(['$\gamma=$' num2str(gamma(aa))],'interpreter','latex');
    f1.Children(1).TickLabelInterpreter = 'latex';
end

legend('$m=0$','$m=1$','$m=2$','$m=3$','$n=0$','$n=1$','$n=2$','$n=3$','interpreter','latex','location','south','orientation','horizontal','numcolumns',4);
f1.Units = 'inches';
f1.Position = [1 1 5.5 6];
f1.Children(5).Position(2) = f1.Children(5).Position(2)+0.0375;
f1.Children(4).Position(2) = f1.Children(4).Position(2)+0.0375;
f1.Children(3).Position(2) = f1.Children(3).Position(2)+0.05;
f1.Children(2).Position(2) = f1.Children(2).Position(2)+0.05;
f1.Children(1).Position(2) = f1.Children(1).Position(2)-0.15;
f1.Children(1).Position(1) = (1-f1.Children(1).Position(3))/2;

saveas(f1,'dispersion_sound_angle.eps','epsc');


