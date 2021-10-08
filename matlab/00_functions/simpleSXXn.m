function [sxx,frequency] = simpleSXXn(x,sampleRate,window)
% Check Number of Inputs - Select Default Options
switch nargin
    case 1
        window = 1;
        sampleRate = ones(1,ndims(x));
    case 2
        window = 1;
end
% Check if 'sampleRate' is Empty
if isempty(sampleRate)
    sampleRate = ones(1,ndims(x));
end
% Check if 'window' is a Function Handle
N = size(x);
if isa(window,'function_handle')
    wfun = window;
    window = wfun(N(1));
    for aa=2:ndims(x)
        window = window.*permute(wfun(N(aa)),aa:-1:1);
    end
end
% Calculate Window Correction
if window==1
    cw = 1;
else
    cw = 1/sqrt(sum(window.^2,'all')/numel(x));
end
% Calculate Power Spectra
sxx = fftshift((abs(cw*fftn(x.*window))).^2)/numel(x)/prod(sampleRate);
% Calculate Frequency Ranges
frequency{1} = (-0.5:1/N(1):0.5-1/N(1))'*sampleRate(1);
for aa=2:ndims(x)
    frequency{aa} = permute((-0.5:1/N(aa):0.5-1/N(aa))'*sampleRate(aa),aa:-1:1);
end
end

