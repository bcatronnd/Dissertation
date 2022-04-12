function [Xb,Yb,Zb,dl] = coordBEAM(lx,dl,gamma,ap,ap_n)
% Recalculate dl to fit
l = lx/sind(gamma);
dl = l/ceil(l/dl);
% Aperature Coordinates
[Xap,Yap] = coordAP(ap,ap_n);
% Generate Beam Coordinates
Xb = repmat(Xap,1,1,l/dl+1);
Yb = repmat(Yap,1,1,l/dl+1);
Zb = repmat(reshape(0:dl:l,1,1,[]),ap_n,ap_n)+Xb*cotd(gamma);
end