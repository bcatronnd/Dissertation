close all; clc; clearvars;

load('synthetic_wavefront.mat');

c = 340;
M = 0.6;
uBL_u = 0.8;

wf_f = WFfilter(wf.wf,'velocity-lowpass',[c*M*uBL_u*wf.sampleRate(1)/wf.sampleRate(3),0,0.0125]);

log_range = -14.5;
BlockSize = 2^10;
% Dispersion - Filtered WF
[WF,freq] = simpleDispersion(wf_f,'BlockSize',BlockSize,'SampleRate',wf.sampleRate);

out = [mean(nanrms(reshape(wf.wf,size(wf.wf,1)*size(wf.wf,2),[])),2)
mean(nanrms(reshape(wf.AO,size(wf.wf,1)*size(wf.wf,2),[])),2)
mean(nanrms(reshape(wf_f,size(wf.wf,1)*size(wf.wf,2),[])),2)];

disp(['Ratio Unfiltered: ' num2str(out(1)/out(2),'%0.3f')]);
disp(['Ratio Filtered:   ' num2str(out(3)/out(2),'%0.3f')]);

%%%%% Plot
views = [-125 25; -55 25; 180 0; 270 0];
f1 = figure(1);
for aa=1:4
    subplot(2,2,aa);
    scolor = parula(2);
    patch(isocaps(freq{1},freq{2},freq{3}(end/2+1:end),WF(:,:,end/2+1:end),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
    patch(isosurface(freq{1},freq{2},freq{3},WF,log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
    grid on;
    daspect([1 1 50]);
    xlim(wf.sampleRate(1)/2*[-1 1]);
    ylim(wf.sampleRate(2)/2*[-1 1]);
    zlim(wf.sampleRate(3)/2*[0 1]);
    xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
    ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
    zlabel('$f\ (Hz)$','Interpreter','Latex');
    view(views(aa,:));
    material dull;
    camlight;
    f1.Children(1).TickLabelInterpreter = 'latex';
end
f1.Units = 'inches';
f1.Position = [1 1 5.5 7.5];
saveas(f1,'filter_velocity.eps','epsc');

%%%%% Animation
nFrames = 150;
theta = (0:nFrames-1)/(nFrames)*4*pi;
az = rad2deg(theta);
el = 25*sin(theta/2);

f2 = figure(2);
scolor = parula(2);
patch(isocaps(freq{1},freq{2},freq{3}(end/2+1:end),WF(:,:,end/2+1:end),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq{1},freq{2},freq{3},WF,log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
grid on;
daspect([1 1 50]);
xlim(wf.sampleRate(1)/2*[-1 1]);
ylim(wf.sampleRate(2)/2*[-1 1]);
zlim(wf.sampleRate(3)/2*[0 1]);
xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
zlabel('$f\ (Hz)$','Interpreter','Latex');
f2.Children(1).TickLabelInterpreter = 'latex';
f2.Units = 'inches';
f2.Position = [1 1 5.5 6.25];
cl = camlight;

filename = 'filter_velocity.gif';
frameRate = 15;
for aa=1:nFrames-1
    view(az(aa),el(aa));
    camlight(cl);
    drawnow;
    frame = getframe(f2);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if aa==1
        imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',1/frameRate);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/frameRate);
    end
end



