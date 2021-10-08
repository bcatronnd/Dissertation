function [sxx,freq] = simpleSXX(x,fsamp,window)
N = length(x);
x = reshape(x,1,N);
switch nargin
    case 1
        window = reshape(hann(N),1,N);
        fsamp = 1;
    case 2
        window = reshape(hann(N),1,N);
end
if isempty(fsamp)
    fsamp = 1;
end
if isa(window,'function_handle')
    window = reshape(window(N),1,N);
end
cw = 1/sqrt(sum(window.^2,'all')/N);
sxx = fftshift((abs(cw*fft(x.*window))).^2)/N/fsamp;
freq = (-0.5:1/N:0.5-1/N)*fsamp;
end

