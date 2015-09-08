
%implicit adi method to solve the heat equation part for both equations.
%this method is brocken!!
%% load parameters.
Dn    = 0.1;
length = 3;
tend = 20;

dt = 0.005;
J = 100;
idx = 2:(J-1);

dx = length/J
dy = dx;
mux = Dn*(dt/(dx^2));
muy = Dn*(dt/(dx^2));


steps = floor(tend/dt)

%set initial condition:
%there are only cells at the center.
[x,y] = meshgrid(linspace(0,1,J));
%N = 15*(x - x.^2).*(y-y.^2).*exp(-50 .*((x - 0.5).^2 + (y - 0.5).^2 ));

N = 0.*x + 0.*y;
%N((J/2-10):(J/2+10),(J/2-10):(J/2+10)) = 1;

N = 0.*x + 0.*y;
%i
N(25:70,24:25) = 1;
N(73:78,23:26) = 1;
%G
N(25:75,34:35) = 1;
N(25:26,34:50) = 1;
N(74:75,34:50) = 1;
N(25:45,48:50) = 1;
N(40:45,40:50) = 1;

%E
N(25:70,54:55) = 1;
N(25:26,54:65) = 1;
N(50:51,54:65) = 1;
N(70:71,54:65) = 1;

%M
N(25:78,70:71) = 1;
N(25:78,85:86) = 1;
for i = 1:7
   N(78-i,71+i) = 1;
   N(79-i,71+i) = 1;
end
for i = 1:7
   N(79-i,86-i) = 1;
   N(78-i,86-i) = 1;
end
N(80:100,20:60) = 1;

NHalf = N;

%degradation
Kdeg = 2;

%% Define step matrices:
xStepRight = toeplitz([(1 - mux) 0.5*mux zeros(1,J-2)],...
             [(1- mux) 0.5*mux zeros(1,J-2)]);
yStepRight = toeplitz([(1 - muy) 0.5*muy zeros(1,J-2)],...
             [(1 - muy) 0.5*muy zeros(1,J-2)]);
xStepLeft  = toeplitz([(1 + mux) -0.5*mux zeros(1,J-2)],...
             [(1 + mux) -0.5*mux zeros(1,J-2)]);
yStepLeft  = toeplitz([(1 + muy) -0.5*muy zeros(1,J-2)],...
             [(1 + muy) -0.5*muy zeros(1,J-2)]);

%periodic boundary conditions:
xStepRight(1,J) = 0.5*mux; xStepRight(end,1) = 0.5*mux;
yStepRight(1,J) = 0.5*mux; yStepRight(end,1) = 0.5*mux;
xStepLeft(1,J) = -0.5*mux; xStepLeft(end,1) = -0.5*mux;
yStepLeft(1,J) = -0.5*mux; yStepLeft(end,1) = -0.5*mux;

vidObjRev=VideoWriter('ADISIM.avi');
set(vidObjRev,'FrameRate',24);
open(vidObjRev);

surf(N);shading('flat');%view(0,90);
writeVideo(vidObjRev,getframe(gcf));

%% start computations:
for t = 1:steps
    %boundary condtiotions are always zero!
    NHalf = xStepLeft\(yStepRight*N);
    N = (xStepRight*NHalf)/yStepLeft;
    surf(N);shading('flat');%view(0,90);
    zlim([0 1]);
    writeVideo(vidObjRev,getframe(gcf));
end

vidObjRev.close();
