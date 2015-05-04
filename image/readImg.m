clear all;

%Pattern recognition in matlab.

data = imread('pattern.jpeg');
data2 = rgb2gray(data);
data2 = data2(100:300,250:450);

posmin = data2 %< 200;

img = mat2gray(posmin);


figure
imshow(img)

%-------------------------------------------------------------------------

figure

tend = 0.2;
dt = 0.01;
J = 40;

%x in [0,1].
dx = 1/J;
%y in [0,1];
dy = 1/J;
%equal spacing in x and y direction.
mu = dt^2/dx^2

[x,y] = meshgrid(linspace(0,1,J));

%the boundary conditions are zero (of homogeneous diriclet type).
U1 = zeros(J);
U2 = zeros(J);

%the initial solution is u0(x,y) = sin(pi x) sin(pi y).
U = 15*(x - x.^2).*(y-y.^2).*exp(-50 .*((x - 0.5).^2 + (y - 0.5).^2 ));

Uold = U;
for t = 1:(tend/dt)
    elements = 2:J-1;  
    for i = 1:1:J
        %compute the columns where x is const.
        U1(elements,i) = mu*U(elements+1,i) + mu*U(elements-1,i); 
        %compute the columns where y is const.
        U2(i,elements) = mu*U(i,elements+1) + mu*U(i,elements-1);
    end
    Unew = (2 - 4*mu) .* U - Uold + U1 + U2;
    Uold = U;
    U = Unew;
end

surf(x,y,U); view(2);
        