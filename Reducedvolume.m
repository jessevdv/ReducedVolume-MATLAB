clear all;clc;

x1=0.15; % fiber radius
x2=1.584/2; % the maximum radius of the drop

d=36.447; % contact angle

C=cos(d*pi/180);
n=x2/x1;

a=(x2*C- x1)/(x2-x1*C);
k=(x2^2-(a*x1)^2);

disp(k)
disp(a)

[F,E]= ellipke(k^2); % the first and second kind elliptic integrals
L=2*(a*F*x1+x2*E);
L1=L/x1;

A = 2*a*a+3*a*n+2*n*2;
B = [(n^2-1)*(1-a^2)]^(1/2);

V1=(2*pi*n/3)*[A*E-(a^2)*F+B]-pi*L1;% the reduced volume

disp(V1)










