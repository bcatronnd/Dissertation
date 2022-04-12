function [name] = ZernikeName(j,scheme)
%ZERNIKENAME Returns the classical aberration name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Returns the classical aberration name for the j index of a scheme.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% j [vector] = Index of Zernike polynomial
% scheme [string] = Ordering scheme name:
%    'noll'    - Noll's sequential indices
%    'ansi'    - OSA/ANSI standard indices
%    'fringe'  - Fringe/University of Arizona indices
%    'zemax'   - Zemax standard indicies (Same as 'noll')
%    'osa'     - The Optical Society standard indicies (Same as 'ansi')
%    'annular' - Noll's sequential indices for annular pupils
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% name [cell] = Strings corresponding to the name of index j
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
if (~strcmp(scheme,'noll') && ~strcmp(scheme,'ansi')) && ~strcmp(scheme,'fringe')
    warning('MATLAB:ZernikeName','Invalid Scheme Name. Using Noll''s sequential indices.');
    scheme = 'noll';
end
s_name = {'piston','horizontal tilt','vertical tilt','defocus',...
    'oblique primary astigmatism','vertical primary astigmatism',...
    'vertical coma','horizontal coma','vertical trefoil',...
    'oblique trefoil','primary spherical',...
    'vertical secondary astigmatism','oblique secondary astigmatism',...
    'vertical quadrafoil','oblique quadrafoil',...
    'horizontal secondary coma','vertical secondary coma',...
    'oblique secondary trefoil','vertical secondary trefoil',...
    'oblique pentafoil','vertical pentafoil','-'};
if strcmp(scheme,'noll')
    j(j>21) = 22;
    name = s_name(j);
elseif strcmp(scheme,'ansi')
    [n,m] = ZernikeOrder(j,'noll');
    j = (n.*(n+2)+m)/2;
    j(j>20) = 21;
    name = s_name(j+1);
elseif strcmp(scheme,'fringe')
    s_name = {'piston','tilt x','tilt y','power','astigmatism x',...
        'astigmatism y','coma x','coma y','primary spherical',...
        'trifoil x','trifoil y','secondary astigmatism x'...
        'secondary astigmatism y','secondary coma x','secondary coma y'...
        'secondary spherical','tetrafoil x','tetrafoil y',...
        'secondary trifoil x','secondary trifoil y',...
        'tertiary astigmatism x','tertiary astigmatism y',...
        'tertiary coma x','tertiary coma y','tertiary spherical',...
        'pentafoil x','pentafoil y','secondary tetrafoil x',...
        'secondary tetrafoil y','tertiary trefoil x','tertiary trefoil y',...
        'quatenary astigmatism x','quatenary astigmatism y',...
        'quatenary coma x','quatenary coma y','quatenary spherical','-'};
    j(j>35) = 36;
    name = s_name(j+1);
end
end