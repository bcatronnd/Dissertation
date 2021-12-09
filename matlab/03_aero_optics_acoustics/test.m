close all; clc; clearvars;

A = ones(2,20);
A(2,1) = 1.25;

for aa=2:size(A,2)
    A(:,aa) = A(:,aa-1).*exp(-1i*[14; 16]*pi/180);
    
    
    
end




figure(1)
plot(A','-o')
axis equal;
grid on;
