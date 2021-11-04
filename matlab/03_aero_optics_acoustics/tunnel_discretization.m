close all; clc; clearvars;

tsMach = 0.3;
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
wtz = [];
wts = [];
wteta = [];
wtcross = [];
wtts = [];
for aa=1:length(wtdZ)
    N = ceil((wtZ(aa+1)-wtZ(aa))/dzM);
    wtz = [wtz (0:N-1)/(N)*(wtZ(aa+1)-wtZ(aa))+wtZ(aa)]; %#ok<AGROW>
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
    wts = [wts S]; %#ok<AGROW>
    wteta = [wteta (0:N-1)/(N)*(wtETA(aa+1)-wtETA(aa))+wtETA(aa)]; %#ok<AGROW>
    wtcross = [wtcross wtCROSS(aa)*ones(1,N)]; %#ok<AGROW>
    if aa==tsIndex
        wtts = [wtts ones(1,N)]; %#ok<AGROW>
    else
        wtts = [wtts zeros(1,N)]; %#ok<AGROW>
    end
end
clear aa N x;
wtz = [wtz wtZ(end)];
wts = [wts wtS(end)];
wteta = [wteta wtETA(end)];
wtcross = [wtcross wtCROSS(end)];
wtts = [wtts 0];

% Mach Number
Sstar = wtS(tsIndex)*216/125*tsMach*(1+tsMach^2/5)^-3;
wtm = zeros(1,length(wtz));
for aa=1:length(wtz)
    M = tsMach/2;
    for bb=1:maxIt
        error = Sstar/wts(aa)-216/125*M*(1+M^2/5)^-3;
        if abs(error)<tol
            break
        end
        M = M-error/(1080*(M^2-1)/(M^2+5)^4);
    end
    if bb==maxIt
        disp('Not Converged');
    end
    wtm(aa) = M;
end
clear M error aa bb;

% Sound with Flow
Ma = fliplr(wtm(2:end));
Mb = fliplr(wtm(1:end-1));
Xaa = 1+0.2*Ma.*Ma;
Xab = 1+0.2*Ma.*Mb;
Xbb = 1+0.2*Mb.*Mb;
Yab = 1-0.2*Ma.*Mb;
ptw = fliplr((1+Ma)./(1+Mb).*(2*Mb)./(Mb+Ma).*Xaa./Xab.*(Xaa./Xbb).^2.5);
prw = fliplr((1+Ma)./(1-Ma).*(Mb-Ma)./(Mb+Ma).*Yab./Xab);
clear Ma Mb Xaa Xab Xbb Yab;

% Sound against Flow
Ma = wtm(1:end-1);
Mb = wtm(2:end);
Xaa = 1+0.2*Ma.*Ma;
Xab = 1+0.2*Ma.*Mb;
Xbb = 1+0.2*Mb.*Mb;
Yab = 1-0.2*Ma.*Mb;
pta = fliplr((1-Ma)./(1-Mb).*(2*Mb)./(Ma+Mb).*Xaa./Xab.*(Xaa./Xbb).^2.5);
pra = fliplr((1-Ma)./(1+Ma).*(Mb-Ma)./(Mb+Ma).*Yab./Xab);



f1 = figure(1);
subplot(4,1,1);
plot(wtZ/wtZ(end),wtS,'ok',wtz/wtZ(end),wts,'k-');%,wtz/wtZ(end),wteta,'b-');
grid on;
ylabel('Area ($m^2$)','interpreter','latex');
subplot(4,1,2);
plot(wtz/wtZ(end),wtm,'k-');
grid on;
ylabel('Mach Number','interpreter','latex');
subplot(4,1,3);
plot(wtz(2:end)/wtZ(end),ptw,'k-',wtz(2:end)/wtZ(end),prw,'k:');
grid on;
ylabel({'Sound Transmission';'With Flow'},'interpreter','latex');
subplot(4,1,4);
plot(wtz(1:end-1)/wtZ(end),pta,'k-',wtz(1:end-1)/wtZ(end),pra,'k:');
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
save(['tunnel_discretization_' num2str(tsMach,'%0.2f') '.mat'],'wtcross','wteta','wtm','wts','wtz','wtts','pra','prw','pta','ptw');