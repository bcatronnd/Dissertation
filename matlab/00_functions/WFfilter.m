function [wf] = WFfilter(wf,varargin)
%WFFILTER Summary of this function goes here
%   Detailed explanation goes here

% Check Number of Inputs
if mod(nargin,2)==0 || nargin<3
    error('Invalid Number of Inputs');
end
% Zero-Out NaN Values
mask = double(~isnan(wf(:,:,1)));
mask(mask==0) = NaN;
% disp(mask);
wf(isnan(wf)) = 0;
% 3D-FFT
wf = fftshift(fftn(wf));
% Calculate Frequency
BlockSize = size(wf);
Frequency.y = reshape(-1/2:1/BlockSize(1):1/2-1/BlockSize(1),[],1);
Frequency.x = reshape(-1/2:1/BlockSize(2):1/2-1/BlockSize(2),1,[]);
Frequency.t = reshape(-1/2:1/BlockSize(3):1/2-1/BlockSize(3),1,1,[]);
Frequency.rho = sqrt(Frequency.x.^2+Frequency.y.^2);
% Filters
for aa=1:length(varargin)/2
    switch lower(varargin{2*aa-1})
        case 'time-highpass'
            if length(varargin{2*aa})==1
                n = 1;
            else
                n = varargin{2*aa}(2);
            end
            [b,a] = butter(n,varargin{2*aa}(1),'high','s');
            gain = abs(reshape(freqs(b,a,reshape(Frequency.t,1,[])),1,1,[]));
            clear b a;
        case 'time-lowpass'
            if length(varargin{2*aa})==1
                n = 1;
            else
                n = varargin{2*aa}(2);
            end
            [b,a] = butter(n,varargin{2*aa}(1),'low','s');
            gain = abs(reshape(freqs(b,a,reshape(Frequency.t,1,[])),1,1,[]));
            clear b a;
        case {'time-passband' 'time-bandpass'}
            if length(varargin{2*aa})==2
                n = 1;
            else
                n = varargin{2*aa}(3);
            end
            [b,a] = butter(n,varargin{2*aa}(1:2),'bandpass','s');
            gain = abs(reshape(freqs(b,a,reshape(Frequency.t,1,[])),1,1,[]));
            clear b a;
        case {'time-bandstop' 'time-stop'}
            if length(varargin{2*aa})==2
                n = 1;
            else
                n = varargin{2*aa}(3);
            end
            [b,a] = butter(n,varargin{2*aa}(1:2),'stop','s');
            gain = abs(reshape(freqs(b,a,reshape(Frequency.t,1,[])),1,1,[])); 
            clear b a;
        case 'space-highpass'
            if length(varargin{2*aa})==1
                n = 1;
            else
                n = varargin{2*aa}(2);
            end
            gain = sqrt(1./(1+(Frequency.rho/varargin{2*aa}(1)).^(-2*n)));
        case 'space-lowpass'
            if length(varargin{2*aa})==1
                n = 1;
            else
                n = varargin{2*aa}(2);
            end
            gain = sqrt(1./(1+(Frequency.rho/varargin{2*aa}(1)).^(+2*n)));
        case 'velocity-highpass'
            if length(varargin{2*aa})==3
                n = 1;
            else
                n = varargin{2*aa}(4);
            end
            dist = abs(varargin{2*aa}(1)*Frequency.x+varargin{2*aa}(2)*Frequency.y-Frequency.t)/sqrt((varargin{2*aa}(1))^2+(varargin{2*aa}(2))^2+1);
            gain = sqrt(1./(1+(dist/varargin{2*aa}(3)).^(-2*n)));
        case 'velocity-lowpass'
            if length(varargin{2*aa})==3
                n = 1;
            else
                n = varargin{2*aa}(4);
            end
            dist = abs(varargin{2*aa}(1)*Frequency.x+varargin{2*aa}(2)*Frequency.y-Frequency.t)/sqrt((varargin{2*aa}(1))^2+(varargin{2*aa}(2))^2+1);
            gain = sqrt(1./(1+(dist/varargin{2*aa}(3)).^(+2*n)));
        case 'x-space'
            switch length(varargin{2*aa})==1
                case 1
                    n = 1;
                    threeDB = 0;
                case 2
                    n = varargin{2*aa}(2);
                    threeDB = 0;
                case 3
                    n = varargin{2*aa}(2);
                    threeDB = varargin{2*aa}(3);
            end
            if isempty(n)
                n = 1;
            end
            if threeDB
                gain = 1./(1+exp(-n*(Frequency.x-varargin{2*aa}(1)-log(sqrt(2)-1)/n)));
            else
                gain = 1./(1+exp(-n*(Frequency.x-varargin{2*aa}(1))));
            end
        case 'y-space'
            switch length(varargin{2*aa})==1
                case 1
                    n = 1;
                    threeDB = 0;
                case 2
                    n = varargin{2*aa}(2);
                    threeDB = 0;
                case 3
                    n = varargin{2*aa}(2);
                    threeDB = varargin{2*aa}(3);
            end
            if isempty(n)
                n = 1;
            end
            if threeDB
                gain = 1./(1+exp(-n*(Frequency.y-varargin{2*aa}(1)-log(sqrt(2)-1)/n)));
            else
                gain = 1./(1+exp(-n*(Frequency.y-varargin{2*aa}(1))));
            end
        case 'forward-moving'
            kx = BlockSize(2);
            kt = BlockSize(3);
            if ~isempty(varargin{2*aa})
                kx = kx*varargin{2*aa}(1);
                kt = kt*varargin{2*aa}(2);
            end
            gain = (2./(1+exp(-kx*Frequency.x))-1).*(2./(1+exp(-kt*Frequency.t))-1)/2+0.5;
        case 'backward-moving'
            kx = BlockSize(2);
            kt = BlockSize(3);
            if ~isempty(varargin{2*aa})
                kx = kx*varargin{2*aa}(1);
                kt = kt*varargin{2*aa}(2);
            end
            gain = (2./(1+exp(kx*Frequency.x))-1).*(2./(1+exp(-kt*Frequency.t))-1)/2+0.5;
        case 'forward-moving-ideal'
            gain = sign(Frequency.x.*Frequency.t)/2+0.5;
        case 'backward-moving-ideal'
            gain = -sign(Frequency.x.*Frequency.t)/2+0.5;
        case {'unity' 'no-filter'}
            gain = 1;
        otherwise
            warning('Invalid Filter Type.  Setting Gain to Unity.');
            gain = 1;
    end
%     disp(max(gain,[],'all'));
%     disp(min(gain,[],'all'));
    wf = wf.*gain;
    clear gain n;
end
% 3D-iFFT
wf = real(ifftn(ifftshift(wf))).*mask;
end

