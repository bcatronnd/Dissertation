close all; clc; clearvars;

nPts = 2.^[6 14];
offset = 4;

cn = dsp.ColoredNoise(0.5,nPts(2),nPts(1));
x = cn()';

% g = normpdf(-nPts(1)/2:nPts(1)/2-1,0,1)'*2+1;
% x = x.*g;

w1 = hann(nPts(2))';
cw1 = 1/sqrt(sum(w1.^2,'all')/numel(w1));
X1 = mean((abs(cw1*fftshift(fft(x.*w1,[],2),2))).^2/nPts(2),1);
w2 = hann(nPts(1)).*hann(nPts(2))';
cw2 = 1/sqrt(sum(w2.^2,'all')/numel(w2));
X2 = (abs(cw2*fftshift(fft2(x.*w2)))).^2/numel(x);
F = -1/2:1/nPts(2):1/2-1/nPts(2);

dF = 1./nPts;

X1 = X1(end/2+1:end);
X2 = X2(:,end/2+1:end);
F = F(end/2+1:end);

f1 = figure(1);
loglog(F,dF(1)*sum(X2,1),F,X1);
legend('sum*df','1D')
grid on;


disp(['Error (sum): ' num2str(sum(((X1-dF(1)*sum(X2,1))./X1).^2))]);

% X1'\max(X2,[],1)'
% 
% X = 6:11;
% Y = [4.6883 5.3577 6.1286 6.6367 7.4252 8.1732];
% p = polyfit(X,Y,1);
% f2 = figure(2);
% plot(X,Y);
% grid on;


