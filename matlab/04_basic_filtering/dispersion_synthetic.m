close all; clc; clearvars;

load('synthetic_wavefront.mat');

%%%%% CODE 
BlockSize = 2^10;
[WF,freq] = simpleDispersion(wf.wf,'BlockSize',BlockSize,'SampleRate',wf.sampleRate);

log_range = -14.5;
disp(['Total Energy Inside Surface: ' num2str(sum(10.^WF(WF>log_range),'all')/sum(10.^WF,'all')*100,'%0.1f') '%']);

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
saveas(f1,'dispersion_synthetic.eps','epsc');

%%%%% Animation
% nFrames = 150;
% theta = (0:nFrames-1)/(nFrames)*4*pi;
% az = rad2deg(theta);
% el = 25*sin(theta/2);
% 
% f2 = figure(2);
% scolor = parula(2);
% patch(isocaps(freq{1},freq{2},freq{3}(end/2+1:end),WF(:,:,end/2+1:end),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
% patch(isosurface(freq{1},freq{2},freq{3},WF,log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
% grid on;
% daspect([1 1 50]);
% xlim(wf.sampleRate(1)/2*[-1 1]);
% ylim(wf.sampleRate(2)/2*[-1 1]);
% zlim(wf.sampleRate(3)/2*[0 1]);
% xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
% ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
% zlabel('$f\ (Hz)$','Interpreter','Latex');
% f2.Children(1).TickLabelInterpreter = 'latex';
% f2.Units = 'inches';
% f2.Position = [1 1 5.5 6.25];
% cl = camlight;
% 
% filename = 'dispersion_synthetic.gif';
% frameRate = 15;
% for aa=1:nFrames-1
%     view(az(aa),el(aa));
%     camlight(cl);
%     drawnow;
%     frame = getframe(f2);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if aa==1
%         imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',1/frameRate);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/frameRate);
%     end
% end

