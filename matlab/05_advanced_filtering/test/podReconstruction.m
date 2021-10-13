function [u] = podReconstruction(a,Phi,DataLocation,M,N,dim)

u = NaN*ones(prod(M),N);
u(DataLocation) = Phi*a';

% Reshape and Inverse Permute
u = reshape(u,[M N]);
if dim~=ndims(u)
    u = ipermute(u,mpermute);
end
end

