close all; clc; clearvars;

% Inputs
tsMach = 0.6;
freq = 519.2;
fanMode.m = 8;
fanMode.n = 0;
fanMaxDp = 1;
phaseSteps = 50;
showPlots = 0;
nPointsRect = 51;
filename = ['tunnel_acoustic_against_' num2str(tsMach) '_' num2str(fanMode.m) '_' num2str(fanMode.n) '_tsMid.gif'];
frameRate = 15;
plotLocations = [81];

% Code
load(['tunnel_discretization_' num2str(tsMach,'%0.2f') '.mat']);
phase = (0:phaseSteps-1)/phaseSteps*2*pi;
tsEntr = find(WindTunnel.testsection,1,'first');
tsExit = find(WindTunnel.testsection,1,'last');
tsMid = round(tsEntr/2+tsExit/2);
w = 2*pi*freq;
tsCoeff = complex(zeros(tsExit-tsEntr+1,63,phaseSteps));
beam = complex(zeros(nPointsRect,phaseSteps));

for aa=1:length(phase)
    % Fan Inlet
    load(['transition_modes/eta_' num2str(WindTunnel.eta(1),'%0.4f') '.mat']);
    km = max(besseljdzero(fanMode.m,fanMode.n+1,0==fanMode.m))/sqrt(WindTunnel.area(1)/pi);
    k0 = w/WindTunnel.c(1);
    if imag((WindTunnel.mach(1)*k0+sqrt(k0^2-(1-WindTunnel.mach(1)^2)*km.^2))/(1-WindTunnel.mach(1)^2))~=0
        disp('warning: cutoff')
    end
    fanAcoustic = (besselj(fanMode.m,km*transition_mode.meshd.NodesR/sqrt(1/pi)*sqrt(WindTunnel.area(1)/pi)).*cos(fanMode.m*transition_mode.meshd.NodesT+phase(aa))/besselj(fanMode.m,km*sqrt(WindTunnel.area(1)/pi))*fanMaxDp/2)';
    circle_mode = transition_mode;
    P = fanAcoustic;
    quarryPoints = true(size(P));
    % Spatial Stepping
    cc = 1;
    for bb=1:tsMid
%         if showPlots && any(plotLocations==bb)
%             close all; %#ok<*UNRCH>
%             f1 = figure(1);
%             if WindTunnel.eta(bb)==1 && WindTunnel.eta(bb-1)==1
%                 pdeplot(rectangle_mode.mesh.Nodes,rectangle_mode.mesh.Elements,'XYData',real(P));
%             else
%                 pdeplot(transition_mode.meshd.Nodes,transition_mode.meshd.Elements,'XYData',real(P));
%             end
%             if WindTunnel.cruciform(bb)==1
%                 hold on;
%                 plot([0.5 0.5 NaN 0 1],[0 1 NaN 0.5 0.5],'k-');
%             end
%             colormap(redblue);
%             axis equal off;
%             caxis(fanMaxDp/2*[-1 1]);
%             colorbar('off');
%             [N,D] = rat(phase(aa)/pi);
%             title(['$\theta$=' num2str(N) '$\pi$/' num2str(D)],'interpreter','latex');
%             f1.Units = 'inches';
%             f1.Position = [1 1 3 3];
%             saveas(f1,['tunnel_slices/tunnel_acoustic_against_' num2str(tsMach) '_' num2str(fanMode.m) '_' num2str(fanMode.n) '_' num2str(phase(aa)) '_' num2str(bb,'%0.3u') '.eps'],'epsc')
%             saveas(f1,['tunnel_slices/tunnel_acoustic_against_' num2str(tsMach) '_' num2str(fanMode.m) '_' num2str(fanMode.n) '_' num2str(phase(aa)) '_' num2str(bb,'%0.3u') '.png'],'png')
%         end
%         if bb==tsMid
%             beam(:,aa) = sum(reshape(P,nPointsRect,nPointsRect),2);
%         end
        switch floor(WindTunnel.eta(bb))
            case 0 % Circle to Square Transition
                load(['transition_modes/eta_' num2str(WindTunnel.eta(bb),'%0.4f') '.mat']);
                quarryPoints = and(quarryPoints,~isnan(transition_mode.eigvd(:,1)));
                coeff = (transition_mode.eigvd(quarryPoints,:)\P(quarryPoints))';
                k0 = w/WindTunnel.c(bb);
                km = sqrt(transition_mode.k2'/WindTunnel.area(bb));
                kzm = (WindTunnel.mach(bb)*k0+sqrt(k0^2-(1-WindTunnel.mach(bb)^2)*km.^2))/(1-WindTunnel.mach(bb)^2);
                P = WindTunnel.pta(bb)*(exp(1i*kzm*WindTunnel.z(bb))*(coeff.*transition_mode.eigvd)')';
            case 1 % Square Duct
                switch floor(WindTunnel.eta(bb-1))
                    case 0 % Transition to Square Grid
                        load(['transition_modes/rect_' num2str(nPointsRect) '.mat']);
                        quarryPoints = and(quarryPoints,~isnan(transition_mode.eigvd(:,1)));
                        coeff = (transition_mode.eigvd(quarryPoints,:)\P(quarryPoints))';
                        k0 = w/WindTunnel.c(bb);
                        km = sqrt(rectangle_mode.k2'/WindTunnel.area(bb));
                        kzm = (WindTunnel.mach(bb)*k0+sqrt(k0^2-(1-WindTunnel.mach(bb)^2)*km.^2))/(1-WindTunnel.mach(bb)^2);
                        P = WindTunnel.pta(bb)*(exp(1i*kzm*WindTunnel.z(bb))*(coeff.*rectangle_mode.eigv)')';
                    case 1 % Square Duct
                        switch WindTunnel.cruciform(bb)
                            case 0 % No Cruciform
                                coeff = (rectangle_mode.eigv\P)';
                                k0 = w/WindTunnel.c(bb);
                                km = sqrt(rectangle_mode.k2'/WindTunnel.area(bb));
                                kzm = (WindTunnel.mach(bb)*k0+sqrt(k0^2-(1-WindTunnel.mach(bb)^2)*km.^2))/(1-WindTunnel.mach(bb)^2);
                                P = WindTunnel.pta(bb)*(exp(1i*kzm*WindTunnel.z(bb))*(coeff.*rectangle_mode.eigv)')';
                            case 1 % Cruciform
                                k0 = w/WindTunnel.c(bb);
                                km = sqrt(rectangle_mode.k2'/(WindTunnel.area(bb)/4));
                                kzm = (WindTunnel.mach(bb)*k0+sqrt(k0^2-(1-WindTunnel.mach(bb)^2)*km.^2))/(1-WindTunnel.mach(bb)^2);
                                switch WindTunnel.cruciform(bb-1)
                                    case 0 % Transition to Cruciform
                                        P1 = P;
                                        P1(rectangle_mode.mesh.Nodes(1,:)'<0.5) = NaN;
                                        P1(rectangle_mode.mesh.Nodes(2,:)'<0.5) = NaN;
                                        coeff = (reshape(rectangle_mode.eigv1(repmat(~isnan(P1),1,transition_mode.n)),[],transition_mode.n)\P1(~isnan(P1)))';
                                        P1 = WindTunnel.pta(bb)*(exp(1i*kzm*WindTunnel.z(bb))*(coeff.*rectangle_mode.eigv1)')';
                                        P2 = P;
                                        P2(rectangle_mode.mesh.Nodes(1,:)'>0.5) = NaN;
                                        P2(rectangle_mode.mesh.Nodes(2,:)'<0.5) = NaN;
                                        coeff = (reshape(rectangle_mode.eigv2(repmat(~isnan(P2),1,transition_mode.n)),[],transition_mode.n)\P2(~isnan(P2)))';
                                        P2 = WindTunnel.pta(bb)*(exp(1i*kzm*WindTunnel.z(bb))*(coeff.*rectangle_mode.eigv2)')';
                                        P3 = P;
                                        P3(rectangle_mode.mesh.Nodes(1,:)'>0.5) = NaN;
                                        P3(rectangle_mode.mesh.Nodes(2,:)'>0.5) = NaN;
                                        coeff = (reshape(rectangle_mode.eigv3(repmat(~isnan(P3),1,transition_mode.n)),[],transition_mode.n)\P3(~isnan(P3)))';
                                        P3 = WindTunnel.pta(bb)*(exp(1i*kzm*WindTunnel.z(bb))*(coeff.*rectangle_mode.eigv3)')';
                                        P4 = P;
                                        P4(rectangle_mode.mesh.Nodes(1,:)'<0.5) = NaN;
                                        P4(rectangle_mode.mesh.Nodes(2,:)'>0.5) = NaN;
                                        coeff = (reshape(rectangle_mode.eigv4(repmat(~isnan(P4),1,transition_mode.n)),[],transition_mode.n)\P4(~isnan(P4)))';
                                        P4 = WindTunnel.pta(bb)*(exp(1i*kzm*WindTunnel.z(bb))*(coeff.*rectangle_mode.eigv4)')';
                                    case 1 % Cruciform
                                        P1(rectangle_mode.mesh.Nodes(1,:)'<0.5) = NaN;
                                        P1(rectangle_mode.mesh.Nodes(2,:)'<0.5) = NaN;
                                        coeff = (reshape(rectangle_mode.eigv1(repmat(~isnan(P1),1,transition_mode.n)),[],transition_mode.n)\P1(~isnan(P1)))';
                                        P1 = WindTunnel.pta(bb)*(exp(1i*kzm*WindTunnel.z(bb))*(coeff.*rectangle_mode.eigv1)')';
                                        P2(rectangle_mode.mesh.Nodes(1,:)'>0.5) = NaN;
                                        P2(rectangle_mode.mesh.Nodes(2,:)'<0.5) = NaN;
                                        coeff = (reshape(rectangle_mode.eigv2(repmat(~isnan(P2),1,transition_mode.n)),[],transition_mode.n)\P2(~isnan(P2)))';
                                        P2 = WindTunnel.pta(bb)*(exp(1i*kzm*WindTunnel.z(bb))*(coeff.*rectangle_mode.eigv2)')';
                                        P3(rectangle_mode.mesh.Nodes(1,:)'>0.5) = NaN;
                                        P3(rectangle_mode.mesh.Nodes(2,:)'>0.5) = NaN;
                                        coeff = (reshape(rectangle_mode.eigv3(repmat(~isnan(P3),1,transition_mode.n)),[],transition_mode.n)\P3(~isnan(P3)))';
                                        P3 = WindTunnel.pta(bb)*(exp(1i*kzm*WindTunnel.z(bb))*(coeff.*rectangle_mode.eigv3)')';
                                        P4(rectangle_mode.mesh.Nodes(1,:)'<0.5) = NaN;
                                        P4(rectangle_mode.mesh.Nodes(2,:)'>0.5) = NaN;
                                        coeff = (reshape(rectangle_mode.eigv4(repmat(~isnan(P4),1,transition_mode.n)),[],transition_mode.n)\P4(~isnan(P4)))';
                                        P4 = WindTunnel.pta(bb)*(exp(1i*kzm*WindTunnel.z(bb))*(coeff.*rectangle_mode.eigv4)')';
                                        P = mean([P1 P2 P3 P4],2,'omitnan');
                                end
                        end
                end
        end
        if WindTunnel.testsection(bb)==1
            tsCoeff(cc,:,aa) = coeff;
            cc = cc+1;
        end
    end
    if showPlots
        f1 = figure(1); %#ok<*UNRCH>
%         subplot(2,1,1);
%         pdeplot(circle_mode.meshd.Nodes,circle_mode.meshd.Elements,'XYData',fanAcoustic);
%         colormap(redblue);
%         axis equal off tight;
%         caxis(fanMaxDp/2*[-1 1]);
%         colorbar('off');
%         subplot(2,1,2);
        if WindTunnel.eta(bb)==1
            pdeplot(rectangle_mode.mesh.Nodes,rectangle_mode.mesh.Elements,'XYData',real(P));
        else
            pdeplot(transition_mode.meshd.Nodes,transition_mode.meshd.Elements,'XYData',real(P));
        end
        colormap(redblue);
        axis equal off tight;
        caxis(fanMaxDp/2*[-1 1]);
        colorbar('off');
        f1.Units = 'inches';
        f1.Position = [1 1 2.5 2.5];
        f1.Color = [1 1 1];
        drawnow;
        frame = getframe(f1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if aa==1
            imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',1/frameRate);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/frameRate);
        end
    end
end
tsCoeff = cat(3,tsCoeff,tsCoeff(:,:,1));
phase = cat(2,phase,2*pi);


%%
close all;
tol = 1e-2;
pModes = any(abs(tsCoeff)>tol,[1 3]);
mModes = rectangle_mode.m(pModes);
nModes = rectangle_mode.n(pModes);
legString = cell(1,3);
for aa=1:sum(pModes)
    legString{aa} = ['m=' num2str(mModes(aa)) '  n=' num2str(nModes(aa))];
end
PlotColors = linspecer(length(legString));
pModes = find(pModes);

f2 = figure(2);
for aa=1:length(legString)
    plot(abs(phase)/pi,squeeze(real(tsCoeff(6,pModes(aa),:)))','linewidth',1.25,'color',PlotColors(aa,:));
    hold on;
end
grid on;
xlabel('Initial Phase ($\pi$ rad)','interpreter','latex');
ylabel('$C_m$','interpreter','latex');
legend(legString,'interpreter','latex');
f2.Children(2).TickLabelInterpreter = 'latex';
f2.Units = 'inches';
f2.Position = [1 1 5 3];
% f2.Color = 'none';
saveas(f2,['tunnel_acoustic_against_' num2str(tsMach) '_' num2str(fanMode.m) '_' num2str(fanMode.n) '.eps'],'epsc')
% saveas(f2,['tunnel_acoustic_against_' num2str(tsMach) '_' num2str(fanMode.m) '_' num2str(fanMode.n) '.png'],'png')
% exportgraphics(f2,['tunnel_acoustic_against_' num2str(tsMach) '_' num2str(fanMode.m) '_' num2str(fanMode.n) '.png'],'BackgroundColor','none');

% f2 = figure(2);
% subplot(2,1,1);
% plot(phase/pi,squeeze(abs(tsCoeff(6,pModes,:)))');
% grid on;
% ylabel('$|C_m|$','interpreter','latex');
% legend(legString,'interpreter','latex');
% subplot(2,1,2);
% plot(phase/pi,squeeze(angle(tsCoeff(6,pModes,:)))'/pi);
% grid on;
% xlabel('Initial Phase ($\pi$ rad)','interpreter','latex');
% ylabel('$\angle C_m$ ($\pi$ rad)','interpreter','latex');
% % sgtitle({['M=' num2str(tsMach)];['$\hat{p}^0$(m=' num2str(fanMode.m) ',n=' num2str(fanMode.n) ')']},'interpreter','latex');
% f2.Children(1).TickLabelInterpreter = 'latex';
% f2.Children(3).TickLabelInterpreter = 'latex';
% f2.Units = 'inches';
% f2.Position = [1 1 6 4];
% 
% % saveas(f2,['tunnel_acoustic_against_' num2str(tsMach) '_' num2str(fanMode.m) '_' num2str(fanMode.n) '.eps'],'epsc')
% % saveas(f2,['tunnel_acoustic_against_' num2str(tsMach) '_' num2str(fanMode.m) '_' num2str(fanMode.n) '.png'],'png')
% 
% %%
% f3 = figure(3);
% surf(real(beam),'linestyle','none');
% view(2);
% colormap(redblue);
% caxis(max(abs(beam),[],'all')*[-1 1]);
% axis off;
