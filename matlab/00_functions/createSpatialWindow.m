function [window] = createSpatialWindow(Mask_WF)
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

