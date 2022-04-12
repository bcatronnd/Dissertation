function [zp] = ZernikePolynomial(j,rho,theta,scheme,epsilon)
%ZERNIKEPOLYNOMIAL Computes the jth Zernike Polynomial for a given scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Computes the jth Zernike Polynomial at the points rho and theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% j [scalar] = Index of Zernike polynomial
% rho [matrix] = Radial coordinates
% theta [matrix] = Azimuthal coordinates
% scheme [string] = Ordering scheme name:
%    'noll'    - Noll's sequential indices
%    'ansi'    - OSA/ANSI standard indices
%    'fringe'  - Fringe/University of Arizona indices
%    'zemax'   - Zemax standard indicies (Same as 'noll')
%    'osa'     - The Optical Society standard indicies (Same as 'ansi')
%    'annular' - Noll's sequential indices for annular pupils
% epsilon [scalar] = For annular pupils - epsilon=R_innner/R_outer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% zp [matrix] = Zernike Polynomial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCES:
% > Robert J. Noll. Zernike Polynomials and Atmospheric Turbulence.
% Journal of the Optical Society of America, 1976.
% > Larry N. Thibos, Raymond A. Applegate, James T. Schwiegerling, and
% Robert Webb. Standards for Reporting the Optical Aberrations
% of Eyes. Journal of Refractive Surgery, 2002.
% > James C. Wyant and Katherine Creath. Basic Wavefront Aberration Theory
% for Optical Metrology. Applied Optics and Optical Engineering, 1992.
% > Virendra N. Mahajan. Zernike Annular Polynomials for Imaging Systems
% with Annular Pupils. Journal of the Optical Society of America, 1981.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:
% Brian Catron - bcatron@nd.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERSION:
% 1.0 - 30 AUG 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<4
    warning('MATLAB:ZernikeScheme','No Scheme Given. Using Noll''s sequential indices.');
    scheme = 'noll';
end
if strcmp(scheme,'annular') && nargin<5
    error('MATLAB:ZernikeNoEpsilon','No Epsilon Provided')
end
if strcmp(scheme,'zemax')
    scheme = 'noll';
end
if strcmp(scheme,'osa')
    scheme = 'ansi';
end
if (~strcmp(scheme,'noll') && ~strcmp(scheme,'ansi')) && (~strcmp(scheme,'fringe') && ~strcmp(scheme,'annular'))
    warning('MATLAB:ZernikeName','Invalid Scheme Name. Using Noll''s sequential indices.');
    scheme = 'noll';
end
if sum(size(rho)~=size(theta),'all')
    error('MATLAB:ZernikeInvalidRT','Invalid Input.  rho and theta arrays are not the same size.')
end
[n,m] = ZernikeOrder(j,scheme);
mask = (rho<=1);
zp = Inf*zeros(size(rho));
if strcmp(scheme,'noll')
    m = abs(m);
    Rnm = zeros(size(rho(mask)));
    for ss=0:(n-m)/2
        Rnm = Rnm+((-1)^ss*factorial(n-ss))/(factorial(ss)*factorial((n+m)/2-ss)*factorial((n-m)/2-ss))*rho(mask).^(n-2*ss);
    end
    if m==0
        zp(mask) = sqrt(n+1)*Rnm;
    elseif mod(j,2)==0
        zp(mask) = sqrt(n+1)*Rnm*sqrt(2).*cos(m*theta(mask));
    elseif mod(j,2)==1
        zp(mask) = sqrt(n+1)*Rnm*sqrt(2).*sin(m*theta(mask));
    end
elseif strcmp(scheme,'ansi')
    Rnm = zeros(size(rho(mask)));
    for ss=0:(n-abs(m))/2
        Rnm = Rnm+((-1)^ss*factorial(n-ss))/(factorial(ss)*factorial((n+abs(m))/2-ss)*factorial((n-abs(m))/2-ss))*rho(mask).^(n-2*ss);
    end
    Nnm = sqrt(2*(n+1)/(1+(m==0)));
    if m>=0
        zp(mask) = Nnm*Rnm.*cos(m*theta(mask));
    else
        zp(mask) = -Nnm*Rnm.*sin(m*theta(mask));
    end
elseif strcmp(scheme,'fringe')
    Qnm = zeros(size(rho(mask)));
    for ss=0:(n-abs(m))
        Qnm = Qnm+((-1)^ss*factorial(2*n-abs(m)-ss))/(factorial(ss)*factorial(n-ss)*factorial(n-abs(m)-ss))*rho(mask).^(2*(n-abs(m)-ss));
    end
    if m==0
        zp(mask) = Qnm;
    elseif m>0
        zp(mask) = Qnm.*rho(mask).^m.*cos(m*theta(mask));
    else
        zp(mask) = Qnm.*rho(mask).^abs(m).*sin(abs(m)*theta(mask));
    end
elseif strcmp(scheme,'annular')
    m = abs(m);
    if epsilon<0 || epsilon>=1
        error('MATLAB:ZernikeWrongEpsilon','Epsilon Out Of Range.  0<=epsilon<1')
    end
    mask = ((epsilon<=rho) & (rho<=1));
    eta = 1e-9;
    num = [100 150 200];
    delta = zeros(1,100);
    for aa=1:100
        r = linspace(epsilon,1,num(aa));
        drho = r(2)-r(1);
        R = zeros(size(r));
        for ss=0:(n-m)/2
            R = R+((-1)^ss*factorial(n-ss))/(factorial(ss)*factorial((n+m)/2-ss)*factorial((n-m)/2-ss))*r.^(n-2*ss);
        end
        Nnm = sqrt((1-epsilon^2)/(1+(m==0))/trapz(drho,R.^2.*r));
        if aa>1
            delta(aa) = (Nnm-N_old)/Nnm;
        end
        N_old = Nnm;
        if aa>2
            if abs(delta(aa))<eta
                break
            end
            num(aa+1) = round(num(aa)-delta(aa)/(delta(aa)-delta(aa-1))*(num(aa)-num(aa-1)));
        end
    end
    Rnm = zeros(size(rho(mask)));
    for ss=0:(n-m)/2
        Rnm = Rnm+((-1)^ss*factorial(n-ss))/(factorial(ss)*factorial((n+m)/2-ss)*factorial((n-m)/2-ss))*rho(mask).^(n-2*ss);
    end
    if m==0
        zp(mask) = Nnm*Rnm;
    elseif mod(j,2)==0
        zp(mask) = Nnm*Rnm.*cos(m*theta(mask));
    elseif mod(j,2)==1
        zp(mask) = Nnm*Rnm.*sin(m*theta(mask));
    end
end
end
