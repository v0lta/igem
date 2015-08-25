close all;clear all;clc;

k1=5;
k2=15;

r0=5;

r=0:0.1:14;

F1=k2*(2*r0-r);
%F2=k1*(k2/k1*r0-r);
F2=k1*((k1+k2)/k1*r0-r);

figure(1);
hold on;
plot(r,F1);
plot(r,F2);
n=length(r);
plot(r,zeros(1,n))

legend('F1','F2');
