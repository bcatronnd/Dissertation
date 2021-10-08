close all; clc; clearvars;

nDim = 3;

x = randn(2.^(4+randi(6,1,nDim)));
w = hann(size(x,1));
% for aa=2:nDim
%     w = w.*permute(hann(size(x,aa)),aa:-1:1);
% end



[sxx,freq] = simpleSXXn(x,[],@hann);

for aa=1:length(freq)
    df(aa) = freq{aa}(2)-freq{aa}(1);
end

s1 = sum(x.^2,'all')/numel(x);
s2 = prod(df)*sum(sxx,'all');
