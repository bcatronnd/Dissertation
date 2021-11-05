close all; clc; clearvars;

nLenslets = 32;
lambda = 0.25:0.025:10;
dist = 5:5:25;
zMax = 5;
dz = 0.05;
c0 = 340;
kgd = 2.27e-4;
phaseSteps = 25;
AP = 0.25;

dist = dist*AP;
lambda = lambda*AP;

[x,y,z] = meshgrid(0.975*AP*(-0.5:1/(nLenslets-1):0.5),0.975*AP*(-0.5:1/(nLenslets-1):0.5),-zMax:dz:zMax);
r = sqrt(x(:,:,1).^2+y(:,:,1).^2);
t = atan2(y(:,:,1),x(:,:,1));
mask = ones(nLenslets);
mask(r>0.5*AP) = NaN;
phase = (0:phaseSteps-1)/phaseSteps*2*pi;


OPDrms = zeros(length(lambda),length(dist));
OPDrmsMax = OPDrms;
OPDrmsMin = OPDrms;

for aa=1:length(lambda)
    k = 2*pi/lambda(aa);
    for bb=1:length(dist)
        R = sqrt((x+dist(bb)).^2+y.^2+z.^2);
        OPD = kgd/c0^2*sum(1./R.*exp(-1i*k*R),3)*dz;
        stats = zeros(1,phaseSteps);
        for cc=1:phaseSteps
            wf = mask.*real(OPD.*exp(1i*phase(cc)));
            wf = ZernikeRemoval(1:3,wf,r.*mask,t.*mask,'noll');
            stats(cc) = nanrms(wf(:));
        end
        OPDrms(aa,bb) = mean(stats);
        OPDrmsMax(aa,bb) = max(stats);
        OPDrmsMin(aa,bb) = min(stats);
    end
end




%%
close all;

%fit
X = repmat(lambda'/AP,1,5);
Y = OPDrms*1e6.*sqrt(dist/AP);
[fitresult, gof] = createFit(X, Y);

f1 = figure(1);
plot(lambda/AP,OPDrms*1e6.*sqrt(dist/AP),'o');
hold on;
plot(lambda/AP,fitresult(lambda/AP),'k-');
hold off;
grid on;
xlabel('$\Lambda/Ap$','interpreter','latex');
ylabel('$\frac{OPD_{RMS}\sqrt{R/Ap}}{|A_O|}\ (\frac{\mu m}{kg/s^2})$','interpreter','latex','fontsize',14);
for aa=1:length(dist)
    sLegend{aa} = ['R/Ap = ' num2str(dist(aa)/AP)]; %#ok<SAGROW>
end
sLegend{aa+1} = 'Fitted Curve';
legend(sLegend,'interpreter','latex');
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 5.5 3.5];

saveas(f1,'spherical_sample.eps','epsc');
% save('spherical_sample.mat','fitresult');

% fileID = fopen('spherical_sample.txt','w');
% fprintf(fileID,'\\begin{tabular}{c r}\n');
% fprintf(fileID,'\\hline \n');
% fprintf(fileID,'Coefficent & Value \\\\ \n');
% fprintf(fileID,'\\hline \n');
% coeffs = coeffnames(fitresult);
% for aa=1:length(coeffs)
%     fprintf(fileID,['$' coeffs{aa}(1) '_' coeffs{aa}(2) '$ & ' num2str(fitresult.(coeffs{aa}),'%0.3e') ' \\\\ \n']);
% end
% fprintf(fileID,'\\hline \n');
% fprintf(fileID,'\\end{tabular}\n');
% 
% diary spherical_sample_disp.txt;
% diary on;
% disp(fitresult);
% disp(' ');
% disp(gof);
% diary off;


function [fitresult, gof] = createFit(X, Y)
[xData, yData] = prepareCurveData( X, Y );

% Set up fittype and options.
ft = fittype( 'rat22' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.49809429119639 0.900852488532005 0.574661219130188 0.845178185054037 0.738640291995402];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
end



