function [a,Phi,DataLocation,M,N,dim] = podCompute(u,dim)

if nargin==1
    dim = ndims(u);
end
if ~any(size(dim))
    dim = ndims(u);
end

% Permute 'u' if
if dim~=ndims(u)
    mpermute = [1:dim-1 dim+1:ndims(u) dim];
    u = permute(u,mpermute);
end
for aa=1:ndims(u)-1
    M(aa) = size(u,aa);
end
N = size(u,ndims(u));

% Reshape and Data Location
u = reshape(u,prod(M),N);
DataLocation = ~isnan(u);
u = reshape(u(DataLocation),[],N);

% Perform POD
[Phi,~] = eig(1/N*(u*u'),'vector');
Phi = flip(Phi,2);
a = (u'*Phi);
end