close all; clc; clearvars;

% Inputs
tsMach = 0.5;
freq = 519.2;
fanMode.m = 16;
fanMode.n = 0;
fanMaxDp = 1;
phaseSteps = 25;
showPlots = 1;
c = 340;


% Code
load(['tunnel_discretization_' num2str(tsMach,'%0.2f') '.mat']);
phase = (0:phaseSteps-1)/phaseSteps*2*pi;
tsEntr = find(wtts,1,'first');
tsExit = find(wtts,1,'last')+1;
w = 2*pi*freq;
k0 = w/c;

for aa=1%:length(phase)
    % Fan Inlet 
    load(['transition_modes/eta_' num2str(wteta(1),'%0.4f') '.mat']);
    km = max(besseljdzero(fanMode.m,fanMode.n+1,0==fanMode.m))/sqrt(wts(1)/pi);
    fanAcoustic = (besselj(fanMode.m,km*transition_mode.meshd.NodesR/sqrt(1/pi)*sqrt(wts(1)/pi)).*cos(fanMode.m*transition_mode.meshd.NodesT+phase(aa))/besselj(fanMode.m,km*sqrt(wts(1)/pi))*fanMaxDp/2)';
    if showPlots
        f1 = figure(1); %#ok<*UNRCH>
        subplot(1,2,1);
        pdeplot(transition_mode.meshd.Nodes,transition_mode.meshd.Elements,'XYData',fanAcoustic);
        colormap(redblue);
        axis equal off;
        caxis(fanMaxDp/2*[-1 1]);
        colorbar('off');
        drawnow;
    end
    %
    P = fanAcoustic;
%     bb=1;
%     kzm = (wtm(bb)*k0+sqrt(k0^2-(1-wtm(bb)^2)*km))/(1-wtm(bb)^2);
%     P = real(P*pta(bb)*exp(-1i*wtz(bb+1)*kzm));
    
    for bb=1%:tsExit
        kzm = (wtm(bb)*k0+sqrt(k0^2-(1-wtm(bb)^2)*km))/(1-wtm(bb)^2); % Add - in front of first wtm(bb) for with
        P = P*pta(bb);
        if bb==1
            
            
        end
        
        
        
        % P = real(P*pta(bb)*exp(-1i*wtz(bb+1)*kzm));
        if wteta(bb)~=1
            load(['transition_modes/eta_' num2str(wteta(bb+1),'%0.4f') '.mat']);
            coeff = mldivide(transition_mode.eigvd(~isnan(transition_mode.eigvd(:,1)),:),P(~isnan(transition_mode.eigvd(:,1))))';
            km = sqrt(transition_mode.k2'/wts(bb+1));
        end
        
        
        
        
        
        
        
        if showPlots&&bb==2
            subplot(1,2,2);
            pdeplot(transition_mode.meshd.Nodes,transition_mode.meshd.Elements,'XYData',P);
            colormap(redblue);
            axis equal off;
            caxis(fanMaxDp/2*[-1 1]);
            colorbar('off');
            drawnow;
        end
        
        
    end
    
    
    
    
    
    
    
    
    

end



