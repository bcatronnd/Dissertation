function d_zero = besseljdzero(nu,k_max,allow_x0)
switch nargin
    case 1
        k_max = 1;
        allow_x0 = 0;
    case 2
        allow_x0 = 0;
end

d_zero = NaN*ones(1,k_max);
tol = 1e-12;
int_max = 1e2;

if allow_x0 && abs((besselj(nu-1,0)-besselj(nu+1,0))/2)<tol
    d_zero(1) = 0;
    x_zero = besseljzero(nu,k_max);
    index = 2;
else
    x_zero = besseljzero(nu,k_max+1);    
    index = 1;
end

for aa=1:length(x_zero)
    if index>k_max
        break
    end
    d_zero(index) = x_zero(aa)-pi/2;
    for bb=1:int_max
        fp = (besselj(nu-1,d_zero(index))-besselj(nu+1,d_zero(index)))/2;
        if abs(fp)<tol
            if d_zero(index)<1e-6
                if allow_x0 && index==1
                    d_zero(index) = 0;
                else
                    break
                end
            end
            index = index+1;
            break
        end
        fpp = (besselj(nu-2,d_zero(index))-2*besselj(nu,d_zero(index))+besselj(nu+2,d_zero(index)))/4;
        fppp = (besselj(nu-3,d_zero(index))-3*besselj(nu-1,d_zero(index))+3*besselj(nu+1,d_zero(index))-besselj(nu-3,d_zero(index)))/8;
        d_zero(index) = d_zero(index)-(2*fp*fpp)/(2*fpp^2-fp*fppp);
    end
end
end