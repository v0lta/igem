clear all;close all;clc;
%implicit adi method to solve the heat equation part for both equations.
%this method is brocken!!

%% load parameters.
Dn    = 0.1;
length = 3;
tend = 0.5;

dt = 0.005;
J = 100;
idx = 2:(J-1);

dx = length/J
dy = dx;
mu = Dn*(dt/(dx^2));

steps = floor(tend/dt)

%there are only cells at the center.
[x,y] = meshgrid(linspace(0,1,J));
N = 15*(x - x.^2).*(y-y.^2).*exp(-50 .*((x - 0.5).^2 + (y - 0.5).^2 ));
figure(1);
surf(N);shading('flat');

%N = 0.*x + 0.*y;
%N((J/2-10):(J/2+10),(J/2-10):(J/2+10)) = 1;
Nnew = N;

stepMat1 = toeplitz([(1 + 2*mu) -mu zeros(1,J-4)],...
             [(1 + 2*mu) -mu zeros(1,J-4)]);
stepMat2 = toeplitz([(1 - 2*mu) +mu zeros(1,J-4)],...
             [(1 - 2*mu) +mu zeros(1,J-4)]);
stepMat3 = toeplitz([(2*mu) -mu zeros(1,J-4)],...
             [(2*mu) -mu zeros(1,J-4)]);

%% start computations:
for t = 1:200
  for i = 2:(J-1)  
    %boundary condtiotions are always zero!
    Nnew(i,idx) = inv(stepMat1)*(stepMat2*N(idx,i));
    N(idx,i) = inv(stepMat1)*(Nnew(i,idx)' - stepMat3*N(idx,i));  
  end
    
end 
figure(2);
surf(N);shading('flat');
