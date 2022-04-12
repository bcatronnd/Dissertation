function x_zero = besseljzero(nu,k_max)
switch nargin
    case 1
        k_max = 1;
end

x_zero = NaN*ones(1,k_max);
tol = 1e-12;
int_max = 1e2;

mc = [+0.261036290 -0.258214500 -0.661353966 +0.329767986 -2.615458011 -0.659600335 +0.999988048];
bc = [+1.222492511 -4.133039658 +0.334123561 +0.069952674 -1.988359458 +0.214146406 +0.100341186];
m = mc(1)*(abs(nu)-mc(2))^mc(3)+mc(4)*(abs(nu)-mc(5))^mc(6)+mc(7);
b = bc(1)*(abs(nu)-bc(2))^bc(3)+bc(4)*(abs(nu)-bc(5))^bc(6)+bc(7);

for aa=1:k_max
    if aa==1
        x_zero(1) = m*abs(nu)+b;
    else
        x_zero(aa) = x_zero(aa-1)+pi;
    end
    for bb=1:int_max
        f = besselj(nu,x_zero(aa));
        if abs(f)<tol
            break
        end
        fp = (besselj(nu-1,x_zero(aa))-besselj(nu+1,x_zero(aa)))/2;
        fpp = (besselj(nu-2,x_zero(aa))-2*besselj(nu,x_zero(aa))+besselj(nu+2,x_zero(aa)))/4;
        x_zero(aa) = x_zero(aa)-(2*f*fp)/(2*fp^2-f*fpp);
    end
end
end