function [Xd,Yd,Zd,dl] = coordDUCT(yoff,lx,dl,gamma,ap,ap_n)
[Xb,Yb,Zb,dl] = coordBEAM(lx,dl,gamma,ap,ap_n);
% Transform to Duct Coordinates
Xd = Xb*cosd(gamma)-Zb*sind(gamma);
Yd = Yb+yoff;
Zd = Xb*sind(gamma)+Zb*cosd(gamma);
end