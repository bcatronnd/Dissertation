function [Xap,Yap,Rap,Tap,Rnorm] = coordAP(ap,ap_n)
ap_ratio = 0.975;
[Xap,Yap] = meshgrid(ap/2*(-1:2/(ap_n-1):1)*ap_ratio);
Rap = sqrt(Xap.^2+Yap.^2);
Tap = atan2(Yap,Xap);
Mask = ones(size(Rap));
Mask(Rap>ap/2) = NaN;
Xap = Xap.*Mask;
Yap = Yap.*Mask;
Rap = Rap.*Mask;
Tap = Tap.*Mask;
Rnorm = Rap/(ap/2);
end