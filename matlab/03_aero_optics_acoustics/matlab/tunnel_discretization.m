close all; clc; clearvars;

tsMach = 0.6;
tsTemp = 302;
tsIndex = 7;
wtZ = convertUnits('inch2meter',cumsum([0 95.2 84.5 160.8 81.2 142.6 144 108 91.8 48.5 88.1 122.7 87.6 290.4]));
wtS = convertUnits('inch2meter',[sqrt(pi/4)*100 85.2 85.2 74.4 74.4 55.3 36 36 91.8 91.8 91.8 111.6 111.6 sqrt(pi/4)*96]).^2;
wtdZ = [1 0 1 0 1 1 0 6 0 0 1 0 1];
wtETA = [0 1 1 1 1 1 1 1 1 1 1 1 1 0];
wtCROSS = [0 0 0 0 0 1 1 0 0 0 0 0 0 0];
dzM = 0.25;
tol = 1e-6;
maxIt = 10;

% Interpolate Tunnel Area
WindTunnel.z = [];
WindTunnel.area = [];
WindTunnel.eta = [];
WindTunnel.cruciform = [];
WindTunnel.testsection = [];
for aa=1:length(wtdZ)
    N = ceil((wtZ(aa+1)-wtZ(aa))/dzM);
    WindTunnel.z = [WindTunnel.z (0:N-1)/(N)*(wtZ(aa+1)-wtZ(aa))+wtZ(aa)];
    if wtdZ(aa)==0
        S = wtS(aa)*ones(1,N);
    elseif wtdZ(aa)==1
        S = ((0:N-1)/(N)*(sqrt(wtS(aa+1))-sqrt(wtS(aa)))+sqrt(wtS(aa))).^2;
    elseif wtdZ(aa)==6
        x = (0:N-1)/(N);
        S = ((10*x.^3-15*x.^4+6*x.^5)*(sqrt(wtS(aa+1))-sqrt(wtS(aa)))+sqrt(wtS(aa))).^2;
    else
        S = zeros(1,N);
    end
    WindTunnel.area = [WindTunnel.area S];
    WindTunnel.eta = [WindTunnel.eta (0:N-1)/(N)*(wtETA(aa+1)-wtETA(aa))+wtETA(aa)]; 
    WindTunnel.cruciform = [WindTunnel.cruciform wtCROSS(aa)*ones(1,N)];
    if aa==tsIndex
        WindTunnel.testsection = [WindTunnel.testsection ones(1,N)];
    else
        WindTunnel.testsection = [WindTunnel.testsection zeros(1,N)];
    end
end
clear aa N x;
WindTunnel.z = [WindTunnel.z wtZ(end)];
WindTunnel.area = [WindTunnel.area wtS(end)];
WindTunnel.eta = [WindTunnel.eta wtETA(end)];
WindTunnel.cruciform = [WindTunnel.cruciform wtCROSS(end)];
WindTunnel.testsection = [WindTunnel.testsection 0];

% Mach Number
Sstar = wtS(tsIndex)*216/125*tsMach*(1+tsMach^2/5)^-3;
WindTunnel.mach = zeros(1,length(WindTunnel.z));
for aa=1:length(WindTunnel.z)
    M = tsMach/2;
    for bb=1:maxIt
        error = Sstar/WindTunnel.area(aa)-216/125*M*(1+M^2/5)^-3;
        if abs(error)<tol
            break
        end
        M = M-error/(1080*(M^2-1)/(M^2+5)^4);
    end
    if bb==maxIt
        disp('Not Converged');
    end
    WindTunnel.mach(aa) = M;
end
clear M error aa bb;

% Temperature, Speed of Sound, and Velocity
T0 = tsTemp/(1+tsMach^2/5)^-1;
WindTunnel.temp = T0*(1+WindTunnel.mach.^2/5).^-1;
WindTunnel.c = sqrt(1.4*287*WindTunnel.temp);
WindTunnel.u = WindTunnel.c.*WindTunnel.mach;

% Sound with Flow
Ma = fliplr(WindTunnel.mach(2:end));
Mb = fliplr(WindTunnel.mach(1:end-1));
Xaa = 1+0.2*Ma.*Ma;
Xab = 1+0.2*Ma.*Mb;
Xbb = 1+0.2*Mb.*Mb;
Yab = 1-0.2*Ma.*Mb;
WindTunnel.ptw = fliplr((1+Ma)./(1+Mb).*(2*Mb)./(Mb+Ma).*Xaa./Xab.*(Xaa./Xbb).^2.5);
WindTunnel.prw = fliplr((1+Ma)./(1-Ma).*(Mb-Ma)./(Mb+Ma).*Yab./Xab);
clear Ma Mb Xaa Xab Xbb Yab;

% Sound against Flow
Ma = WindTunnel.mach(1:end-1);
Mb = WindTunnel.mach(2:end);
Xaa = 1+0.2*Ma.*Ma;
Xab = 1+0.2*Ma.*Mb;
Xbb = 1+0.2*Mb.*Mb;
Yab = 1-0.2*Ma.*Mb;
WindTunnel.pta = (1-Ma)./(1-Mb).*(2*Mb)./(Ma+Mb).*Xaa./Xab.*(Xaa./Xbb).^2.5;
WindTunnel.pra = (1-Ma)./(1+Ma).*(Mb-Ma)./(Mb+Ma).*Yab./Xab;



f1 = figure(1);
subplot(4,1,1);
plot(wtZ/wtZ(end),wtS,'ok',WindTunnel.z/wtZ(end),WindTunnel.area,'k-');%,wtz/wtZ(end),wteta,'b-');
grid on;
ylabel('Area ($m^2$)','interpreter','latex');
subplot(4,1,2);
plot(WindTunnel.z/wtZ(end),WindTunnel.mach,'k-');
grid on;
ylabel('Mach Number','interpreter','latex');
subplot(4,1,3);
plot(WindTunnel.z(2:end)/wtZ(end),WindTunnel.ptw,'k-',WindTunnel.z(2:end)/wtZ(end),WindTunnel.prw,'k:');
grid on;
ylabel({'Sound Transmission';'With Flow'},'interpreter','latex');
subplot(4,1,4);
plot(WindTunnel.z(1:end-1)/wtZ(end),WindTunnel.pta,'k-',WindTunnel.z(1:end-1)/wtZ(end),WindTunnel.pra,'k:');
grid on;
ylabel({'Sound Transmission';'Against Flow'},'interpreter','latex');
sgtitle(['$M_{TS}=' num2str(tsMach) '$'],'interpreter','latex');
f1.Children(2).TickLabelInterpreter = 'latex';
f1.Children(3).TickLabelInterpreter = 'latex';
f1.Children(4).TickLabelInterpreter = 'latex';
f1.Children(5).TickLabelInterpreter = 'latex';
f1.Units = 'inches';
f1.Position = [1 1 5.5 7.5];

saveas(f1,['tunnel_discretization_' num2str(tsMach,'%0.2f') '.eps'],'epsc');
save(['tunnel_discretization_' num2str(tsMach,'%0.2f') '.mat'],'WindTunnel');