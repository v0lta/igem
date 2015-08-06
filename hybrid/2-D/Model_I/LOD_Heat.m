
%implicit LOD method to solve the heat equation part for both equations.
%seems to work. Did not compare with the known analytic solution so far.

%% load parameters.
Dn    = 0.1;
length = 3;
tend = 6;

dt = 0.005;
J = 301;
idx = 2:(J-1);

dx = length/J
dy = dx;
mu = (dt/(dx^2));

steps = floor(tend/dt)

%there are only cells at the center.
[x,y] = meshgrid(linspace(0,1,J));
N = 15*(x - x.^2).*(y-y.^2).*exp(-50 .*((x - 0.5).^2 + (y - 0.5).^2 ));

%N = 0.*x + 0.*y;
%N((J/2-10):(J/2+10),(J/2-10):(J/2+10)) = 1;
Nint = N;

stepMat1 = toeplitz([(1 + mu) -0.5*mu zeros(1,J-4)],...
             [(1 + mu) -0.5*mu zeros(1,J-4)]);
stepMat2 = toeplitz([(1 - mu) 0.5*mu zeros(1,J-4)],...
             [(1 - mu) 0.5*mu zeros(1,J-4)]);

%% start computations:
for t = 1:steps
  %for i = 2:(J-1)  
    %boundary condtiotions are always zero!
    Nint(idx,idx) = inv(stepMat1)*(stepMat2*N(idx,idx));
    N(idx,idx) = inv(stepMat1)*(stepMat2*Nint(idx,idx)');  
  %end
    
end 
surf(N);shading('flat');