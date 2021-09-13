function [sxx,freq] = simpleSXXn(x,fsamp,window)
switch nargin
    case 1
        window = 1;
        fsamp = ones(1,ndims(x));
    case 2
        window = 1;
end
if isempty(fsamp)
    fsamp = ones(1,ndims(x));
end
cw = 1/sqrt(sum(window.^2,'all')/numel(x));
sxx = cw*fftshift((abs(fftn(x.*window))).^2)/numel(x)/prod(fsamp);
for aa=1:ndims(x)
    freq{aa} = (-0.5:1/size(x,aa):0.5-1/size(x,aa))*fsamp(aa);
end
end

