close all; clc; clearvars;

directory = '/media/briancatron/ResearchData/2021-Spectral-Wavefront-Filtering/';
testPoint = '20210901003';

[WF,CineInfo,RunLog,WFInfo] = loadWF([directory testPoint '_WF.mat'],'Scale',1e6);

BlockSize = 2^10;
WF = WF(:,:,1:floor(size(WF,3)/BlockSize)*BlockSize);
WF = cat(1,WF,NaN*zeros(2^nextpow2(size(WF,1))-size(WF,1),size(WF,2),size(WF,3)));
WF = cat(2,WF,NaN*zeros(size(WF,1),2^nextpow2(size(WF,2))-size(WF,2),size(WF,3)));
[WF,freq] = simpleDispersion(WF,'BlockSize',BlockSize,'SampleRate',RunLog.samplerate);

log_range = -15;
disp(['Total Energy Inside Surface: ' num2str(sum(10.^WF(WF>log_range),'all')/sum(10.^WF,'all')*100,'%0.1f') '%']);

%%%%% Plot
f1 = figure(1);
scolor = parula(2);
patch(isocaps(freq{1},freq{2},freq{3}(end/2+1:end),WF(:,:,end/2+1:end),log_range(1),'all'),'facecolor','interp','edgecolor','none','facelighting','none');
patch(isosurface(freq{1},freq{2},freq{3},WF,log_range(1)),'edgecolor','none','facecolor',scolor(1,:),'facelighting','gouraud');%,'specularstrength',0.375);
grid on;
daspect([1 1 50]);
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[-1 1]);
zlim(RunLog.samplerate(3)/2*[0 1]);
xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
zlabel('$f\ (Hz)$','Interpreter','Latex');
material dull;
camlight;
f1.Children(1).TickLabelInterpreter = 'latex';
view(3);
% f1.Units = 'inches';
% f1.Position = [1 1 5.5 7.5];

%%%% Animation
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
xlim(RunLog.samplerate(1)/2*[-1 1]);
ylim(RunLog.samplerate(2)/2*[-1 1]);
zlim(RunLog.samplerate(3)/2*[0 1]);
xlabel('$\xi_x\ (m^{-1})$','Interpreter','Latex');
ylabel('$\xi_y\ (m^{-1})$','Interpreter','Latex');
zlabel('$f\ (Hz)$','Interpreter','Latex');
f2.Children(1).TickLabelInterpreter = 'latex';
f2.Units = 'inches';
f2.Position = [1 1 5.5 6.25];
cl = camlight;

filename = [testPoint '.gif'];
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

