clearvars; close all; clc;

[X,Y] = meshgrid(linspace(-1,1,55)*0.975);
R = sqrt(X.^2+Y.^2);
T = atan2(Y,X);
X(R>1) = NaN;
Y(R>1) = NaN;
T(R>1) = NaN;
R(R>1) = NaN;

%%%%% ZernikeName.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Noll Index')
for jj=1:10
    disp(strjoin([' j=' num2str(jj) ':' ZernikeName(jj,'noll')],''))
end
disp('ANSI Index')
for jj=1:10
    disp(strjoin([' j=' num2str(jj) ':' ZernikeName(jj,'ansi')],''))
end
clear jj

%%%%% Turning Off Warnings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('on','all');

ZernikeName(3);
warning('off','MATLAB:ZernikeScheme');
ZernikeName(3);

ZernikeName(4,'nool');
warning('off','MATLAB:ZernikeName');
ZernikeName(4,'nool');
clear ans

%%%%% ZernikeOrder.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = 5;
[n,m] = ZernikeOrder(index,'noll');
disp(['Noll Index ' num2str(index) ': n=' num2str(n) ',m=' num2str(m)])
[n,m] = ZernikeOrder(index,'ansi');
disp(['ANSI Index ' num2str(index) ': n=' num2str(n) ',m=' num2str(m)])
clear index n m

%%%%% ZernikePolynomial.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
for aa=1:20
    surf(ZernikePolynomial(aa,R,T,'noll'),'linestyle','none')
    view(2)
    axis equal tight off
    title({['$Z_{j}=' num2str(aa) '$'];strjoin(ZernikeName(aa,'noll'))},'interpreter','latex','fontsize',20)
    caxis(5*[-1 1])
    colorbar('ticklabelinterpreter','latex','fontsize',10)
    drawnow
    pause(0.5)
end

%%%%% ZernikeReconstruct.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nWF = 25;
nMODES = 10;
WF = ZernikeReconstruct(1:nMODES,2*rand(nMODES,nWF)-1,R,T,'noll');

%%%%% ZernikeCoefficient.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PISTON = ZernikeCoefficient(1,WF,R,T,'noll');
hTILT = ZernikeCoefficient(2,WF,R,T,'noll');
vTILT = ZernikeCoefficient(3,WF,R,T,'noll');

%%%%% ZernikeRemoval.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WF2 = ZernikeRemoval(1:3,WF,R,T,'noll');

%%%%% ZernikeKeep.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WF3 = ZernikeKeep(1:3,WF,R,T,'noll');

%%%%% Sample Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
for aa=1:nWF
    subplot(1,3,1)
    surf(WF(:,:,aa),'linestyle','none')
    view(2)
    axis equal tight off
    caxis(max(WF,[],'all')*[-1 1])
    title('Original','interpreter','latex','fontsize',15)
    
    subplot(1,3,2)
    surf(WF3(:,:,aa),'linestyle','none')
    view(2)
    axis equal tight off
    caxis(max(WF,[],'all')*[-1 1])
    title('Piston/Tip/Tilt','interpreter','latex','fontsize',15)
    colorbar('south','position',[0.1 0.15 0.8 0.05],'ticklabelinterpreter','latex','fontsize',10)
    
    subplot(1,3,3)
    surf(WF2(:,:,aa),'linestyle','none')
    view(2)
    axis equal tight off
    caxis(max(WF,[],'all')*[-1 1])
    title('Removed','interpreter','latex','fontsize',15)
    
    drawnow
    pause(0.5)
end
clear aa

%%%%% Annular Aperature %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = 1/3;
X(R<epsilon) = NaN;
Y(R<epsilon) = NaN;
T(R<epsilon) = NaN;
R(R<epsilon) = NaN;

figure(3)
for aa=1:20
    surf(ZernikeReconstruct(1:25,2*rand(25,1)-1,R,T,'annular',epsilon),'linestyle','none')
    view(2)
    axis equal tight off
    caxis(5*[-1 1])
    colorbar('ticklabelinterpreter','latex','fontsize',10)
    title('Annular Aperature','interpreter','latex','fontsize',15)
    drawnow
    pause(0.5)
end
clear aa
