function [n,m] = ZernikeOrder(j,scheme)
%ZERNIKEORDER Computes Zernike n and m coefficients from a single-index
%scheme.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Computes the radial and azimuthal order of Zernike polynomial
% using various ordering scheme.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% j [scalar] = Index of Zernike polynomial
% scheme [string] = Ordering scheme name:
%    'noll'    - Noll's sequential indices
%    'ansi'    - OSA/ANSI standard indices
%    'fringe'  - Fringe/University of Arizona indices
%    'zemax'   - Zemax standard indicies (Same as 'noll')
%    'osa'     - The Optical Society standard indicies (Same as 'ansi')
%    'annular' - Noll's sequential indices for annular pupils
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% n [scalar] = Radial order of Zernike polynomial
% m [scalar] = Azimuthal order of Zernike polynomial
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


if nargin<2
    warning('MATLAB:ZernikeScheme','No Scheme Given. Using Noll''s sequential indices.');
    scheme = 'noll';
end
if strcmp(scheme,'zemax') || strcmp(scheme,'annular')
    scheme = 'noll';
end
if strcmp(scheme,'osa')
    scheme = 'ansi';
end
if strcmp(scheme,'ansi')
    if any(j<0)
        error('MATLAB:ZernikeInvalidIndex','Index must be >= 0');
    end
    n = ceil((-3+sqrt(9+8*j))/2);
    m = 2*j-n.*(n+2);
elseif strcmp(scheme,'fringe')
    if any(j<0)
        error('MATLAB:ZernikeInvalidIndex','Index must be >= 0');
    end
    if any(j>36)
        error('MATLAB:ZernikeInvalidIndex','Index must be <=36');
    end
    n_index = [0 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6];
    m_index = [0 1 -1 0 2 -2 1 -1 0 3 -3 2 -2 1 -1 0 4 -4 3 -3 2 -2 1 -1 0 5 -5 4 -4 3 -3 2 -2 1 -1 0 0];
    n = n_index(j+1);
    m = m_index(j+1);
else
    if ~strcmp(scheme,'noll')
        warning('MATLAB:ZernikeName','Invalid Scheme Name. Noll''s sequential indices.');
    end
    if any(j<=0)
        error('MATLAB:ZernikeInvalidIndex','Index must be > 0');
    end
    n = floor((sqrt(8*j-1)-1)/2);
    pD = j-n.*(n+1)/2-1;
    m = (pD+rem(rem(n,2)+pD,2)).*(mod(j-1,2)*2-1);
end
if nargout<2
    if size(j,1)==1
        n=[n; m];
    else
        n=[n m];
    end
end
end

