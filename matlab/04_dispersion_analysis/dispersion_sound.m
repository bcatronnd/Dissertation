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
        plot3(kzp/2/pi,repmat(-kx/2/pi,1,length(f)),f,'linewidth',1.25,'color',PlotColors(aa,:));
        plot3(-kzm/2/pi,repmat(kx/2/pi,1,length(f)),f,'linewidth',1.25,'color',PlotColors(aa,:));
        plot3(-kzm/2/pi,repmat(-kx/2/pi,1,length(f)),f,'linewidth',1.25,'color',PlotColors(aa,:));
    end
end
grid on;
view(145,15);
% xlim(5*[-1 1]);
zlim([0 fRange(end)]);
xlabel('$\xi_x\ (m^{-1})$','interpreter','latex');
ylabel('$\xi_y\ (m^{-1})$','interpreter','latex');
zlabel('$f\ (Hz)$','interpreter','latex');
title(['$m\in [0,' num2str(m(end)) '],\ n=0$'],'interpreter','latex');
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
title(['$m=0,\ n\in [0,' num2str(n(end)) ']$'],'interpreter','latex');
f1.Children(1).TickLabelInterpreter = 'latex';
legend(slegend,'interpreter','latex','location','southeast');%'southoutside','orientation','horizontal');

f1.Units = 'inches';
f1.Position = [1 1 5.5 3];

saveas(f1,'dispersion_sound.eps','epsc');












% kx = m*pi;
% ky = 0;
% f = [0 reshape(c/2/pi*sqrt((1-mach.^2).*(kx(:).^2+ky^2)),1,[])];
% f = permute(unique(sort(f)),[1 3 2]);
% f(2:end) = f(2:end)+1e-6;
% k0 = 2*pi*f/c;
% kzp = (-mach.*k0+sqrt(k0.^2-(1-mach.^2).*(kx'.^2+ky^2)))./(1-mach.^2);
% kzm = (+mach.*k0+sqrt(k0.^2-(1-mach.^2).*(kx'.^2+ky^2)))./(1-mach.^2);
% kzp(imag(kzp)>0) = NaN;
% kzm(imag(kzm)>0) = NaN;
% 
% hold on;
% for aa=1:length(mach)
% plot3(repmat(kx,length(f),1)',repmat(squeeze(f),1,length(kx))',squeeze(kzp(:,aa,:)),'color',PlotColors(aa,:))
% end



% M = 0.5;
% modeRangeM = [0 20];
% modeRangeN = 0:40;
% modeRangeNplot = 0:5:40;
% c = 340;
% 
% [X,Y,Z] = meshgrid(-30:10,-20:20,0:50:5000);
% 
% % Code
% for aa=1:length(modeRangeM)
%     [m,n] = meshgrid(modeRangeM(aa),modeRangeN);
%     kx = m*pi;
%     ky = n*pi;
%     f = [0 permute(c/2/pi*sqrt((1-M^2)*(kx(:).^2+ky(:).^2)),[2 1])+1e-6];
%     f = permute(unique(sort(f)),[1 3 2]);
%     w = 2*pi*f;
%     k0 = w/c;
%     kzp = (-M*k0+sqrt(k0.^2-(1-M^2)*(kx.^2+ky.^2)))/(1-M^2);
%     kzm = (+M*k0+sqrt(k0.^2-(1-M^2)*(kx.^2+ky.^2)))/(1-M^2);
%     kzp(imag(kzp)>0) = NaN;
%     kzm(imag(kzm)>0) = NaN;
%     
%     x = [kzp(:); -kzm(:); -kzm(:); kzp(:)]/2/pi;
%     y = [reshape(repmat(ky,1,1,length(f)),[],1); reshape(repmat(ky,1,1,length(f)),[],1); reshape(repmat(-ky,1,1,length(f)),[],1); reshape(repmat(-ky,1,1,length(f)),[],1)]/2/pi;
%     z = repmat(reshape(repmat(f,length(modeRangeN),length(modeRangeM(aa)),1),[],1),4,1);
%     xnan = isnan(x);
%     x(xnan) = [];
%     y(xnan) = [];
%     z(xnan) = [];
%     
%     shp = alphaShape(x,y,z);
%     C{aa} = inShape(shp,X,Y,Z);
% end
% C = double(xor(C{1},C{2}));
% 
% m = 0;
% n = modeRangeNplot';
% kx = m*pi;
% ky = n*pi;
% f = [0:25:5000 permute(c/2/pi*sqrt((1-M^2)*(kx(:).^2+ky(:).^2)),[2 1])+1e-6];
% f = unique(sort(f));
% w = 2*pi*f;
% k0 = w/c;
% kzp = (-M*k0+sqrt(k0.^2-(1-M^2)*(kx.^2+ky.^2)))/(1-M^2);
% kzm = (+M*k0+sqrt(k0.^2-(1-M^2)*(kx.^2+ky.^2)))/(1-M^2);
% kzp(imag(kzp)>0) = NaN;
% kzm(imag(kzm)>0) = NaN;
% 
% 
% 
% %%
% 
% 
% close all;
% f1 = figure(1);
% patch(isocaps(X,Y,Z,C,1-1e-6,'all'),'facecolor','blue','edgecolor','none','facelighting','gouraud');
% patch(isosurface(X,Y,Z,C,1-1e-6),'edgecolor','none','facecolor','blue','facelighting','gouraud');%,'specularstrength',0.375);
% hold on;
% for aa=1:2
%     m = modeRangeM(aa);
%     n = modeRangeNplot';
%     kx = m*pi;
%     ky = n*pi;
%     f = [0:25:5000 permute(c/2/pi*sqrt((1-M^2)*(kx(:).^2+ky(:).^2)),[2 1])+1e-6];
%     f = unique(sort(f));
%     w = 2*pi*f;
%     k0 = w/c;
%     kzp = (-M*k0+sqrt(k0.^2-(1-M^2)*(kx.^2+ky.^2)))/(1-M^2);
%     kzm = (+M*k0+sqrt(k0.^2-(1-M^2)*(kx.^2+ky.^2)))/(1-M^2);
%     kzp(imag(kzp)>0) = NaN;
%     kzm(imag(kzm)>0) = NaN;
%     plot3(kzp'/2/pi,repmat(ky/2/pi,1,length(f))',repmat(f,length(n),1)','m',kzp'/2/pi,repmat(-ky/2/pi,1,length(f))',repmat(f,length(n),1)','m',-kzm'/2/pi,repmat(ky/2/pi,1,length(f))',repmat(f,length(n),1)','m',-kzm'/2/pi,repmat(-ky/2/pi,1,length(f))',repmat(f,length(n),1)','m');
% end
% view(3);
% grid on;
% daspect([1 1 75]);
% xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
% ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
% zlabel('$f\ (Hz)$','Interpreter','Latex');
% xlim([-30 10]);
% ylim([-20 20]);
% zlim([0 5000]);
% material dull;
% camlight;


