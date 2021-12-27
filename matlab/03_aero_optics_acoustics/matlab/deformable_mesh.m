close all; clc; clearvars;

f1 = figure(1);

load('eta_0.0000.mat');
subplot(1,3,1);
pdemesh(transition_mode.meshd.Nodes,transition_mode.meshd.Elements);
grid on;
axis equal tight off;
rectangle('Position',[-0.5 -0.5 1 1]);
rectangle('Position',sqrt(1/pi)*[-1 -1 2 2],'Curvature',[1 1]);

load('eta_0.5000.mat');
subplot(1,3,2);
pdemesh(transition_mode.meshd.Nodes,transition_mode.meshd.Elements);
grid on;
axis equal tight off;
rectangle('Position',[-0.5 -0.5 1 1]);
rectangle('Position',sqrt(1/pi)*[-1 -1 2 2],'Curvature',[1 1]);

load('eta_1.0000.mat');
subplot(1,3,3);
pdemesh(transition_mode.meshd.Nodes,transition_mode.meshd.Elements);
grid on;
axis equal tight off;
rectangle('Position',[-0.5 -0.5 1 1]);
rectangle('Position',sqrt(1/pi)*[-1 -1 2 2],'Curvature',[1 1]);

f1.Units = 'inches';
f1.Position = [1 1 6 2];

saveas(f1,'deformable_mesh.eps','epsc');