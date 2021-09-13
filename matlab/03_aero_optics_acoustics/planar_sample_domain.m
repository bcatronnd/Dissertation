clearvars; close all; clc;

LAMBDA = 0.75;
SIN_max = 0.1;
APERTURE = 0.4;
THETA = 7.5;

% Field
x = 0.5*linspace(-1,1,200);
y = linspace(0,1,100);
[X,Y] = meshgrid(x,y);
Z = SIN_max*sin(X*2*pi/LAMBDA);

% Beam
theta = linspace(0,2*pi,100);
[Tb,Yb] = meshgrid(theta,y);
Xb = APERTURE*cos(Tb)/2/cosd(THETA)-sind(2*THETA)*Yb;
Zb = APERTURE*sin(Tb)/2;

R = [cosd(-2*THETA) -sind(-2*THETA) 0; sind(-2*THETA) cosd(-2*THETA) 0; 0 0 1];


fig = figure;
hold on
contourf(X,Y,Z,25,'linestyle','none')
surf(Xb,Yb,Zb,'linestyle','none','facecolor','g','facealpha',0.5)
colormap(redblue)
axis equal off
% view(0,25)
view(2)

% Lambda
arrow([-LAMBDA/2 0.05 5],[LAMBDA/2 0.05 5],'ends','both','linewidth',1.5);
t1 = text(0.25,0.09,5,'$\Lambda$','fontsize',16,'horizontalalignment','center','Interpreter','latex','fontsize',12);
% Ap
arrow([-APERTURE/2 0.35 5]*R,[APERTURE/2 0.35 5]*R,'ends','both','linewidth',1.5);
t2 = text(-APERTURE/2,0.3625,5,'$Ap$','fontsize',16,'horizontalalignment','center','Interpreter','latex','fontsize',12);
% ln
arrow([0 0 5],[0 1 5],'ends','both','linewidth',1.5);
t3 = text(0.05,0.6,5,'$l_n$','fontsize',16,'horizontalalignment','center','Interpreter','latex','fontsize',12);
% theta
arrow([0 0.9 5],[0 0.9 5]*R,'ends','both','linewidth',1.5);
t4 = text(-0.125,0.925,5,'$\theta$','fontsize',16,'horizontalalignment','center','Interpreter','latex','fontsize',12);
plot3(-sind(2*THETA)*y,y,0*y+5,'--k','linewidth',1.5)


%%%%% Print Formatting
% figProps = struct(fig);
fig.Units = 'inches';
fig.Position = [1 1 5 2.5];
saveas(fig,'planar_sample_domain.eps','epsc')
