function [WF] = filterLSE_POD(WF,y,range)

% Step 1 - FFT
WF = fft(WF,[],3);
y = fft(y,[],2);

% Step 2 - Compute Spectral POD
M = size(WF,1,2);
N = size(WF,3);
WF = reshape(WF,prod(M),N);
DataLocation = ~isnan(WF);
WF = reshape(WF(DataLocation),[],N);
[Phi,~] = eig(1/N*(WF*WF'),'vector');
Phi = flip(Phi,2);
a = (WF'*Phi);

% Step 3 - LSE
% L = (a'*conj(y'))/(y*conj(y'));
% a = a-(L*y)';
frequency = (-0.5:1/size(WF,3):0.5-1/size(WF,3))';
lPoints = and(abs(frequency)<=range(2),abs(frequency)>=range(1));
L = (a(lPoints',:)'*y(:,lPoints)')/(y(:,lPoints)*y(:,lPoints)');
a_lse = (L*y(:,lPoints))';
a_ao = a;
a_ao(lPoints',:) = a(lPoints',:)-a_lse;

% Step 4
WF = NaN*ones(prod(M),N);
WF(DataLocation) = Phi*a_ao';
WF = reshape(WF,[M N]);

% Step 5
WF = real(ifft(WF,[],3));





end

