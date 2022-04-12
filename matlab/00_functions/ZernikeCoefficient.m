function coeff = ZernikeCoefficient(j,phase_screen,rho,theta,scheme,epsilon)
%ZERNIKECOEFFICIENT Computes the Zernike coefficients of a phase screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Computes the Zernike coefficients of a phase screen array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% j [vector] = Vector of the indicies of Zernike polynomial
% phase_screen [array] = Wavefront/phase screen array
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
% coeff [matrix] = Coefficients for each j and phase screen
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


if nargin<5
    warning('MATLAB:ZernikeScheme','No Scheme Given. Using Noll''s sequential indices.');
    scheme = 'noll';
end
if strcmp(scheme,'annular') && nargin<6
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
mask = (rho<=1);
zp = zeros(length(j),length(rho(mask)));
for aa=1:length(j)
    if strcmp(scheme,'annular')
        mask = ((epsilon<=rho) & (rho<=1));
        if epsilon<0 || epsilon>=1
            error('MATLAB:ZernikeWrongEpsilon','Epsilon Out Of Range.  0<=epsilon<1')
        end
        zp(aa,:) = ZernikePolynomial(j(aa),rho(mask),theta(mask),scheme,epsilon);
    else
        zp(aa,:) = ZernikePolynomial(j(aa),rho(mask),theta(mask),scheme);
    end
end
coeff = zp'\reshape(phase_screen(repmat(mask,1,1,size(phase_screen,3))),size(zp,2),[]);
end