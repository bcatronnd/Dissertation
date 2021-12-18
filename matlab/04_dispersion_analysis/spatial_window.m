close all; clc; clearvars;

n = 2^6;
c = [25 0 0;
    1 11 5;
    0.5 7 10];

c = permute(c,[3 2 1]);
[x,y] = meshgrid(-n/2:n/2-1);
r = sqrt(x.^2+y.^2);
t = atan2(y,x);
a = zeros(n);
R = sum(c(1,1,:).*cos(c(1,2,:).*t+c(1,3,:)),3);
a(r<=R) = 1;

[w,d] = createSpatialWindow(a);


f1 = figure(1);
colormap(gray)
subplot(1,2,1);
surf(a,'linestyle','none');
axis image off tight;
view(2);
subplot(1,2,2);
surf(w,'linestyle','none');
axis image off tight;
view(2);

f1.Units = 'inches';
f1.Position = [1 1 5.5 2.25];

saveas(f1,'spatial_window.eps','epsc');


function [window,maskDist] = createSpatialWindow(Mask_WF)
%CREATESPATIALWINDOW Summary of this function goes here
%   Detailed explanation goes here
maskSize = size(Mask_WF)+2;
mask = zeros(maskSize);
mask(2:end-1,2:end-1) = Mask_WF;
maskDist = zeros(maskSize);
[x,y] = meshgrid(1:maskSize(2),1:maskSize(1));
x = x(mask==0);
y = y(mask==0);
for aa=1:numel(mask)
    [row,col] = ind2sub(maskSize,aa);
    if mask(row,col)~=0
        maskDist(row,col) = min(sqrt((x-col).^2+(y-row).^2));
    end
end
maskDist = maskDist/max(maskDist,[],'all');
window = (cos((1-maskDist)*pi)+1)/2;
window = window(2:end-1,2:end-1);
end