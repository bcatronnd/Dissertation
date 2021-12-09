close all; clc; clearvars;

k2Max = 675;
nMax = 10;
nPoints = 51;

load('transition_modes/eta_1.0000.mat');

% Modes to Generate
[rectangle_mode.m,rectangle_mode.n] = meshgrid(0:nMax);
rectangle_mode.m = rectangle_mode.m(:);
rectangle_mode.n = rectangle_mode.n(:);
rectangle_mode.m = rectangle_mode.m(2:end);
rectangle_mode.n = rectangle_mode.n(2:end);
rectangle_mode.k2 = pi^2*(rectangle_mode.m.^2+rectangle_mode.n.^2);
[rectangle_mode.k2,Ind] = sort(rectangle_mode.k2);
rectangle_mode.m = rectangle_mode.m(Ind);
rectangle_mode.n = rectangle_mode.n(Ind);
rectangle_mode.m(rectangle_mode.k2>k2Max) = [];
rectangle_mode.n(rectangle_mode.k2>k2Max) = [];
rectangle_mode.k2(rectangle_mode.k2>k2Max) = [];

% Rectangle Mesh
[x,y] = meshgrid(0:1/(nPoints-1):1);
DT = delaunayTriangulation(x(:),y(:));
rectangle_mode.mesh.Nodes = DT.Points';
rectangle_mode.mesh.Elements = DT.ConnectivityList';

% Replace Transition Mode Values
transition_mode.k2 = rectangle_mode.k2;
xmin = min(transition_mode.mesh.Nodes(1,:));
ymin = min(transition_mode.mesh.Nodes(2,:));
transition_mode.eigv = (cos(rectangle_mode.m*pi*(transition_mode.mesh.Nodes(1,:)-xmin)).*cos(rectangle_mode.n*pi*(transition_mode.mesh.Nodes(2,:)-ymin)))';
transition_mode.eigvd = (cos(rectangle_mode.m*pi*(transition_mode.meshd.Nodes(1,:)-xmin)).*cos(rectangle_mode.n*pi*(transition_mode.meshd.Nodes(2,:)-ymin)))';

% Generate Rectangle Mode Values
rectangle_mode.eigv = (cos(rectangle_mode.m*pi*rectangle_mode.mesh.Nodes(1,:)).*cos(rectangle_mode.n*pi*rectangle_mode.mesh.Nodes(2,:)))';

% Generate Cruciform Mode Values
rectangle_mode.k2c = 4*rectangle_mode.k2;
rectangle_mode.eigv1 = (cos(rectangle_mode.m*pi*(rectangle_mode.mesh.Nodes(1,:)-0.5)*2).*cos(rectangle_mode.n*pi*(rectangle_mode.mesh.Nodes(2,:)-0.5)*2))';
rectangle_mode.eigv1(repmat(rectangle_mode.mesh.Nodes(1,:)'<0.5,1,transition_mode.n)) = NaN;
rectangle_mode.eigv1(repmat(rectangle_mode.mesh.Nodes(2,:)'<0.5,1,transition_mode.n)) = NaN;
rectangle_mode.eigv2 = (cos(rectangle_mode.m*pi*(rectangle_mode.mesh.Nodes(1,:)-0.0)*2).*cos(rectangle_mode.n*pi*(rectangle_mode.mesh.Nodes(2,:)-0.5)*2))';
rectangle_mode.eigv2(repmat(rectangle_mode.mesh.Nodes(1,:)'>0.5,1,transition_mode.n)) = NaN;
rectangle_mode.eigv2(repmat(rectangle_mode.mesh.Nodes(2,:)'<0.5,1,transition_mode.n)) = NaN;
rectangle_mode.eigv3 = (cos(rectangle_mode.m*pi*(rectangle_mode.mesh.Nodes(1,:)-0.0)*2).*cos(rectangle_mode.n*pi*(rectangle_mode.mesh.Nodes(2,:)-0.0)*2))';
rectangle_mode.eigv3(repmat(rectangle_mode.mesh.Nodes(1,:)'>0.5,1,transition_mode.n)) = NaN;
rectangle_mode.eigv3(repmat(rectangle_mode.mesh.Nodes(2,:)'>0.5,1,transition_mode.n)) = NaN;
rectangle_mode.eigv4 = (cos(rectangle_mode.m*pi*(rectangle_mode.mesh.Nodes(1,:)-0.5)*2).*cos(rectangle_mode.n*pi*(rectangle_mode.mesh.Nodes(2,:)-0.0)*2))';
rectangle_mode.eigv4(repmat(rectangle_mode.mesh.Nodes(1,:)'<0.5,1,transition_mode.n)) = NaN;
rectangle_mode.eigv4(repmat(rectangle_mode.mesh.Nodes(2,:)'>0.5,1,transition_mode.n)) = NaN;

% Save Data
% save(['transition_modes/rect_' num2str(nPoints) '.mat'],'transition_mode','rectangle_mode');




k = 6;

figure(1);
subplot(3,2,1);
pdeplot(transition_mode.meshd.Nodes,transition_mode.meshd.Elements,'XYData',transition_mode.eigvd(:,k));
axis equal off;
colorbar('off');
rectangle('Position',[-0.5 -0.5 1 1]);
subplot(3,2,2);
pdeplot(rectangle_mode.mesh.Nodes,rectangle_mode.mesh.Elements,'XYData',rectangle_mode.eigv(:,k));
axis equal off;
colorbar('off');
rectangle('Position',[0 0 1 1]);
subplot(3,2,3);
pdeplot(rectangle_mode.mesh.Nodes,rectangle_mode.mesh.Elements,'XYData',rectangle_mode.eigv1(:,k));
axis equal off;
colorbar('off');
rectangle('Position',[0 0 1 1]);
subplot(3,2,4);
pdeplot(rectangle_mode.mesh.Nodes,rectangle_mode.mesh.Elements,'XYData',rectangle_mode.eigv2(:,k));
axis equal off;
colorbar('off');
rectangle('Position',[0 0 1 1]);
subplot(3,2,5);
pdeplot(rectangle_mode.mesh.Nodes,rectangle_mode.mesh.Elements,'XYData',rectangle_mode.eigv3(:,k));
axis equal off;
colorbar('off');
rectangle('Position',[0 0 1 1]);
subplot(3,2,6);
pdeplot(rectangle_mode.mesh.Nodes,rectangle_mode.mesh.Elements,'XYData',rectangle_mode.eigv4(:,k));
axis equal off;
colorbar('off');
rectangle('Position',[0 0 1 1]);
colormap(redblue);
sgtitle(num2str(k));
