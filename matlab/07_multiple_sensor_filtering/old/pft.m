function [X] = pft(x,k,n,dim)


switch nargin
    case 2
        dim = find(size(x)~=1,1,'first');
        n = size(x,dim);
    case 3
        dim = find(size(x)~=1,1,'first');
end
if isempty(n)
    n = size(x,dim);
end
xSize = size(x);
xSize(dim) = length(k);
dimorder = [dim 1:dim-1 dim+1:ndims(x)];
x = reshape(permute(x,dimorder),size(x,dim),prod(xSize(dimorder(2:end))));
X = complex(zeros(length(k),size(x,2)));

wn = exp(-2i*pi/n);
j = (0:n-1)';
for aa=1:length(k)
    X(aa,:) = sum(x.*wn.^j.^(k(aa)-1),1);
end
X = ipermute(reshape(X,length(k),xSize(dimorder(2:end))),dimorder);
end

