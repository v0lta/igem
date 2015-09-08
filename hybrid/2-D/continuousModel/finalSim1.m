
%% Explicitly solve the Keller-Segel equations.

%% Load constants.
Da = 0.72*10^(-3);     %1*10^(-6)*3600; %(cm)^2/h;
Db = 2.376*10^(-3);     %5*10^(-6)*3600; %(cm)^2/h;
Dr = 60*10^(-3);        %(cm)^2/h;
Dh = 50*10^(-3);        %(cm)^2/h;

dt = 0.0015; %h
tend = 120; %h
Xlen = 8.0; %cm
J = 100; %#

Kc = 0.15;
gamma = 10^(-4);
As = 1.0 * 10^2; % (cl)^(-1);
Bs = 1.0 * 10^2; % (cl)^(-1);
kr = 1.584*10^(-4);     %nmol/h
kh = 1.5*1.584*10^(-4); %nmol/h
KAbsorbR = 10^(-5);
KAbsorbH = 10^(-5);


dx = Xlen/J;
steps = tend/dt
mu = dt/dx^2

A = zeros(J,J);
B = zeros(J,J);
R = zeros(J,J);
H = zeros(J,J);
preDif = zeros(1,J);
Astep = zeros(size(A));
Bstep = zeros(size(B));
Rstep = zeros(size(R));
Hstep = zeros(size(H));

idx = 2:1:(J-1);

%% Set intial conditions.
A = zeros(J,J);
B = zeros(J,J);

[x,y] = meshgrid(linspace(0,Xlen,J));
A = A + 2.6*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 4).^2 + (y - 4).^2 ));
A = A + 8*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 4).^2 + (y - 2.5).^2 ));
A = A + 8*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 2.5).^2 + (y - 4).^2 ));
A = A + 1.25*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 5.5).^2 + (y - 4).^2 ));
A = A + 1.25*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 4).^2   + (y - 5.5).^2 ));
%A = A + 3.75*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 5.5).^2 + (y - 2.5).^2 ));
%A = A + 0.6*(x - x.^2).*(y-y.^2).*exp(-20 .*((x  -5.5).^2   + (y - 5.5).^2 ));
%A = A + 24*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 2.5).^2 + (y - 2.5).^2 ));
%A = A + 4*(x - x.^2).*(y-y.^2).*exp(-20 .*((x  -2.5).^2   + (y - 5.5).^2 ));


B = B + 2.5*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 4).^2 + (y - 4).^2 ));
B = B + 8*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 4).^2 + (y - 2.5).^2 ));
B = B + 8*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 2.5).^2 + (y - 4).^2 ));
B = B + 1.25*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 5.5).^2 + (y - 4).^2 ));
B = B + 1.25*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 4).^2   + (y - 5.5).^2 ));
%B = B + 3.75*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 5.5).^2 + (y - 2.5).^2 ));
%B = B + 0.6*(x - x.^2).*(y-y.^2).*exp(-20 .*((x  -5.5).^2   + (y - 5.5).^2 ));
%B = B + 24*(x - x.^2).*(y-y.^2).*exp(-20 .*((x - 2.5).^2 + (y - 2.5).^2 ));
%B = B + 4*(x - x.^2).*(y-y.^2).*exp(-20 .*((x  -2.5).^2   + (y - 5.5).^2 ));
R(:,:) = 0.00001;
H(:,:) = 0;

%%set up movie obj.
vidObj=VideoWriter('finalSim9.avi');
set(vidObj,'FrameRate',24);
open(vidObj);


%% run simulation.
for t = 1:1:steps
   % Compute the inner grid values.
   %compute values related to cell A.
   i = idx;
   j = idx;
   [ Astep(i,j),Bstep(i,j),Rstep(i,j),Hstep(i,j) ] =  ...
     modKS( i,j,A,B,R,H,Da,Db,Dr,Dh,dt,J,Kc,gamma,As,Bs,kr,kh, ...
     KAbsorbR,KAbsorbH,dx);
   
    
   %% Do the boundary (i==1) Lower boundary
   i = 1;
   j = idx;
   [ Astep(i,j),Bstep(i,j),Rstep(i,j),Hstep(i,j) ] =  ...
    modKS( i,j,A,B,R,H,Da,Db,Dr,Dh,dt,J,Kc,gamma,As,Bs,kr,kh, ...
    KAbsorbR,KAbsorbH,dx);
                    
                
   %% Take care of (i==J) Top boundary
   i = J;
   j = idx;
   [ Astep(i,j),Bstep(i,j),Rstep(i,j),Hstep(i,j) ] =  ...
    modKS( i,j,A,B,R,H,Da,Db,Dr,Dh,dt,J,Kc,gamma,As,Bs,kr,kh, ...
    KAbsorbR,KAbsorbH,dx);
   

    %% Now (j==1) Left boundary   
   i = idx;
   j = 1;
   [ Astep(i,j),Bstep(i,j),Rstep(i,j),Hstep(i,j) ] =  ...
     modKS( i,j,A,B,R,H,Da,Db,Dr,Dh,dt,J,Kc,gamma,As,Bs,kr,kh, ...
     KAbsorbR,KAbsorbH,dx);
                
   %% (j==J) Right boundary
   i = idx;
   j = J;
   [ Astep(i,j),Bstep(i,j),Rstep(i,j),Hstep(i,j) ] =  ...
     modKS( i,j,A,B,R,H,Da,Db,Dr,Dh,dt,J,Kc,gamma,As,Bs,kr,kh, ...
     KAbsorbR,KAbsorbH,dx);
    
   %% (i==j==1) Bottom Left conrner
   i = 1;
   j = 1;
   [ Astep(i,j),Bstep(i,j),Rstep(i,j),Hstep(i,j) ] =  ...
     modKS( i,j,A,B,R,H,Da,Db,Dr,Dh,dt,J,Kc,gamma,As,Bs,kr,kh, ...
     KAbsorbR,KAbsorbH,dx);
                
   %% i=1,j==J Bottom right cornet
   i = 1;
   j = J;
   [ Astep(i,j),Bstep(i,j),Rstep(i,j),Hstep(i,j) ] =  ...
     modKS( i,j,A,B,R,H,Da,Db,Dr,Dh,dt,J,Kc,gamma,As,Bs,kr,kh, ...
     KAbsorbR,KAbsorbH,dx);
   
  %% i==J j==1 Top left corner
   i = J;
   j = 1;
   [ Astep(i,j),Bstep(i,j),Rstep(i,j),Hstep(i,j) ] =  ...
     modKS( i,j,A,B,R,H,Da,Db,Dr,Dh,dt,J,Kc,gamma,As,Bs,kr,kh,...
     KAbsorbR,KAbsorbH,dx);
                
   %% i==J, j==J
   i = J;
   j = J;
   [ Astep(i,j),Bstep(i,j),Rstep(i,j),Hstep(i,j) ] =  ...
     modKS( i,j,A,B,R,H,Da,Db,Dr,Dh,dt,J,Kc,gamma,As,Bs,kr,kh,...
     KAbsorbR,KAbsorbH,dx);
                
                
    %% Update the arrays!            
    A = Astep;
    B = Bstep;
    R = Rstep;
    H = Hstep;
    
    if (mod(t,250) == 0)
        
        subplot(1,4,1);
        surf(A);shading('flat');xlabel('x');ylabel('y');title('density type A');
        view(0,90);
        hours = dt*t;
        mFigure = gcf;
        % Create a uicontrol of type "text"
        mTextBox = uicontrol('style','text');
        string = [num2str(hours,3) ' h'];
        set(mTextBox,'string',string);
        mTextBox.Position = [700 270 80 20];
        mTextBox.BackgroundColor = [1 1 1];
        mTextBox.FontWeight = 'bold';
        subplot(1,4,2);
        surf(B);shading('flat');xlabel('x');ylabel('y');title('density type B');
        view(0,90);
        zlim([0 1300]);
        subplot(1,4,3);
        surf(R);shading('flat');xlabel('x');ylabel('y');title('repellent');
        view(0,90);
        subplot(1,4,4);
        surf(H);shading('flat');xlabel('x');ylabel('y');title('AHL');
        view(0,90);
        writeVideo(vidObj,getframe(gcf)); 
        
        hours
    end
    
end
vidObj.close();




