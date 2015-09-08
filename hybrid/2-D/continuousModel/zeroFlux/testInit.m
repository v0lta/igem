%Take a look at initial conditions without simulating.

clear;
Xlen = 8.0; %cm
J = 100; %#

A = zeros(J,J);
B = zeros(J,J);

[x,y] = meshgrid(linspace(0,Xlen,J));
A = A + 1*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 7.5).^2 + (y - 7.5).^2 ));
%A = A + 2.5*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 4).^2 + (y - 4).^2 ));
%A = A + 8*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 4).^2 + (y - 2.5).^2 ));
%A = A + 8*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 2.5).^2 + (y - 4).^2 ));
%A = A + 1.25*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 5.5).^2 + (y - 4).^2 ));
%A = A + 1.25*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 4).^2   + (y - 5.5).^2 ));
%A = A + 3.75*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 5.5).^2 + (y - 2.5).^2 ));
%A = A + 0.6*(x - x.^2).*(y-y.^2).*exp(-20 .*((x  -5.5).^2   + (y - 5.5).^2 ));
%A = A + 24*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 2.5).^2 + (y - 2.5).^2 ));
%A = A + 4*(x - x.^2).*(y-y.^2).*exp(-20 .*((x  -2.5).^2   + (y - 5.5).^2 ));

B = B + 1*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 7.5).^2 + (y - 7.5).^2 ));
%B = B + 2.5*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 4).^2 + (y - 4).^2 ));
%B = B + 8*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 4).^2 + (y - 2.5).^2 ));
%B = B + 8*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 2.5).^2 + (y - 4).^2 ));
%B = B + 1.25*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 5.5).^2 + (y - 4).^2 ));
%B = B + 1.25*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 4).^2   + (y - 5.5).^2 ));
%B = B + 3.75*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 5.5).^2 + (y - 2.5).^2 ));
%B = B + 0.6*(x - x.^2).*(y-y.^2).*exp(-20 .*((x  -5.5).^2   + (y - 5.5).^2 ));
%B = B + 24*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 2.5).^2 + (y - 2.5).^2 ));
%B = B + 4*(x - x.^2).*(y-y.^2).*exp(-20 .*((x  -2.5).^2   + (y - 5.5).^2 ));

figure(1)
clf;
subplot(1,2,1)
surf(A);shading('flat');xlabel('x');ylabel('y');title('density type A');
view(0,90);
subplot(1,2,2)
surf(B);shading('flat');xlabel('x');ylabel('y');title('density type B');
view(0,90);

