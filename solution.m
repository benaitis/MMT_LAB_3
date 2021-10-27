clc;
clear all;
close all;

var = mod(20162566,3);

k=1.38*10^-23;
e=1.602*10^-19;

i0=1*10^-6;
T=263.15;

U=linspace(-0.04, 0.16, 100);
I=i0*(exp((U.*e)/(k*T))-1);

plot(U, I);
hold on;

x0=[randn(1) randn(1)*300];
answer_newton=newtons('f46', x0);
answer_fsolve=fsolve('f46', x0);

I_newton = answer_newton(1)*(exp((U.*e)/(k*answer_newton(2)))-1);
I_fsolve = answer_fsolve(1)*(exp((U.*e)/(k*answer_fsolve(2)))-1);

% figure(2)
plot(U, I_newton);
plot(U, I_fsolve);
title('Diodo VACH');
legend('funckija', 'Newton method', 'FSolve');
hold off;

