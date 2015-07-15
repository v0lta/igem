function [c,f,s] = pdechemopde(x,t,u,DuDx)

D1 = 0.1;
D2 = 0.1;
alfa = 0.7;

c = [1; 1];
f1 = D1.*DuDx(1) - u(1)./u(2).*DuDx(2).*2;
f2 = D2.*DuDx(2);
f = [f1; f2];
s = [0; -alfa.*u(1)];