close all; clc; clearvars;

load('tunnel_discretization_0.50.mat','wteta');

% General Mesh Inputs
hmax = 1/32;
evr = [1e-3,675];
order = 6;
theta = (0:2^order-1)/2^(order-1)*pi;
% Defroming Mesh Inputs000
hmaxD = hmax;
etaD = 0.5;

%%%%% Generate The Deformable Mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = cos(theta).*(etaD/2./max(abs([cos(theta); sin(theta)]))+sqrt(1/pi)*(1-etaD));
y = sin(theta).*(etaD/2./max(abs([cos(theta); sin(theta)]))+sqrt(1/pi)*(1-etaD));
model = createpde();
tr = triangulation(polyshape(x,y));
geometryFromMesh(model,tr.Points',tr.ConnectivityList');
specifyCoefficients(model,'m',-1,'d',0,'c',-1,'a',0','f',0);
applyBoundaryCondition(model,'neumann','Edge',1:model.Geometry.NumEdges,'g',0,'q',0);
generateMesh(model,'Hmax',hmaxD,'GeometricOrder','linear');
MeshD = model.Mesh;
% f1 = figure(1);
% pdemesh(MeshD); grid on; axis equal;
% xticks(-0.5:0.25:0.5); yticks(-0.5:0.25:0.5);
% rectangle('Position',0.5*[-1 -1 2 2]);
% rectangle('Position',sqrt(1/pi)*[-1 -1 2 2],'Curvature',[1 1]);
% title('Deformable Mesh','interpreter','latex');
% f1.WindowStyle = 'docked';
% f1.Color = [1 1 1];
% f1.Children.TickLabelInterpreter = 'latex';

for aa=1:length(wteta)
    if isfile(['transition_modes/eta_' num2str(wteta(aa),'%0.4f') '.mat'])
        disp(['ETA = ' num2str(wteta(aa),'%0.4f') ' Already Generated']);
    else
        clear model;
        x = cos(theta).*(wteta(aa)/2./max(abs([cos(theta); sin(theta)]))+sqrt(1/pi)*(1-wteta(aa)));
        y = sin(theta).*(wteta(aa)/2./max(abs([cos(theta); sin(theta)]))+sqrt(1/pi)*(1-wteta(aa)));
        
        model = createpde();
        tr = triangulation(polyshape(x,y));
        geometryFromMesh(model,tr.Points',tr.ConnectivityList');
        specifyCoefficients(model,'m',-1,'d',0,'c',-1,'a',0','f',0);
        applyBoundaryCondition(model,'neumann','Edge',1:model.Geometry.NumEdges,'g',0,'q',0);
        generateMesh(model,'Hmax',hmax);
        MeshN = transformMesh(MeshD,wteta(aa),etaD);
        
        %         f2 = figure(2);
                
        %         pdemesh(MeshN.Nodes,MeshN.Elements); grid on; axis equal;
        %         xticks(-0.5:0.25:0.5); yticks(-0.5:0.25:0.5); ylim(sqrt(1/pi)*[-1 1]);
        %         rectangle('Position',0.5*[-1 -1 2 2]);
        %         rectangle('Position',sqrt(1/pi)*[-1 -1 2 2],'Curvature',[1 1]);
        %         title({'';'Deformed Mesh';['$\eta$ = ' num2str(eta(aa),'%.3f')]},'interpreter','latex');
        %         if aa==1
        %             axis_fix = axis;
        %         else
        %             axis(axis_fix);
        %         end
        %         f2.WindowStyle = 'docked';
        %         f2.Color = [1 1 1];
        %         f2.Children.TickLabelInterpreter = 'latex';
        %         drawnow;
        
        results = solvepdeeig(model,evr);
        Eigenvectors = zeros(size(MeshD.Nodes,2),size(results.Eigenvectors,2));
        for bb=1:size(results.Eigenvectors,2)
            uintrp = interpolateSolution(results,MeshN.Nodes,bb);
            Eigenvectors(:,bb) = uintrp/rms(uintrp,'omitnan');
        end
        
        transition_mode.eta = wteta(aa);
        transition_mode.n = length(results.Eigenvalues);
        transition_mode.k2 = abs(results.Eigenvalues);
        transition_mode.mesh = results.Mesh;
        transition_mode.eigv = results.Eigenvectors;
        transition_mode.meshd = MeshN;
        transition_mode.eigvd = Eigenvectors;
        
        save(['transition_modes/eta_' num2str(wteta(aa),'%0.4f') '.mat'],'transition_mode');
    end
    
    
    
    
    
end






function [MeshN] = transformMesh(MeshN1,etaN,etaN1)
MeshN.NodesT = atan2(MeshN1.Nodes(2,:),MeshN1.Nodes(1,:));
MeshN.NodesR = sqrt((MeshN1.Nodes(1,:)).^2+(MeshN1.Nodes(2,:)).^2).*(etaN/2./max(abs([cos(MeshN.NodesT); sin(MeshN.NodesT)]))+sqrt(1/pi)*(1-etaN))./(etaN1/2./max(abs([cos(MeshN.NodesT); sin(MeshN.NodesT)]))+sqrt(1/pi)*(1-etaN1));
MeshN.Nodes = [MeshN.NodesR.*cos(MeshN.NodesT); MeshN.NodesR.*sin(MeshN.NodesT)];
MeshN.Elements = MeshN1.Elements;
end