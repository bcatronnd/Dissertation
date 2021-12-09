function [tf] = estimateTF(input,output,blocksize,dim)


switch nargin
    case 2
        blocksize = [];
        dim = [];
    case 3
        [~,dim] = max(size(input));
end

tf = output./input;
if ~isempty(blocksize)
    dimorder = 1:ndims(input);
    if dim~=1
        dimorder(dim) = dimorder(1);
        dimorder(1) = dim;
    end
    tf = permute(tf,dimorder);
    tf = [tf; tf(1,:)];
    step = blocksize/size(input,dim);
    tf = interp1((0:step:blocksize)',tf,(0:blocksize)','pchip');
    tf = tf(1:end-1,:);
    tf = ipermute(tf,dimorder)/10^nextpow2(step);
end
end

