
%% Explicitly solve the Keller-Segel equations.

%% Load constants.
Da = 0.72*10^(-3);     %1*10^(-6)*3600; %(cm)^2/h;
Db = 2.376*10^(-3);     %5*10^(-6)*3600; %(cm)^2/h;
Dr = 60*10^(-3);        %(cm)^2/h;
Dh = 50*10^(-3);        %(cm)^2/h;

dt = 0.0015; %h
tend = 60; %h
Xlen = 8.0; %cm
J = 100; %#

Kc = 0.015;
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
[x,y] = meshgrid(linspace(0,Xlen,J));
B = 8*(x - x.^2).*(y-y.^2).*exp(-1 .*((x - 4).^2 + (y - 4).^2 ));
A = 8*(x - x.^2).*(y-y.^2).*exp(-1 .*((x - 4).^2 + (y - 4).^2 ));
R(:,:) = 0.00001;
H(:,:) = 0;

%%set up movie obj.
vidObj=VideoWriter('simulation_contAHLSel.avi');
set(vidObj,'FrameRate',60);
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
    
    if (mod(t,1000) == 0)
        figure(1);clf;
        %subplot(1,4,1);
        %plot(A(J/2,:));
        %hold on;
        %plot(B(J/2,:));
        %plot(R(J/2,:)*10);
        %plot(H(J/2,:)*10);
        %title('parameters in 1d at x=J/2');
        %legend('A','B','rep','AHL')
        %xlabel('y');ylabel('#cells and concentration*10');
        %ylim([0 1250]);
        
        
        
        subplot(1,4,1);
        surf(A);shading('flat');xlabel('x');ylabel('y');title('density type A');
        hours = dt*t;
        mFigure = gcf;
        % Create a uicontrol of type "text"
        mTextBox = uicontrol('style','text');
        string = [num2str(hours,2) ' h'];
        set(mTextBox,'string',string);
        mTextBox.Position = [525 270 80 20];
        mTextBox.BackgroundColor = [1 1 1];
        mTextBox.FontWeight = 'bold';
        zlim([0 1300]);
        subplot(1,4,2);
        surf(B);shading('flat');xlabel('x');ylabel('y');title('density type B');
        zlim([0 1300]);
        subplot(1,4,3);
        surf(R);shading('flat');xlabel('x');ylabel('y');title('repellent');
        zlim([0.0 3]);
        subplot(1,4,4);
        surf(H);shading('flat');xlabel('x');ylabel('y');title('AHL');
        zlim([0.0 3]);
        writeVideo(vidObj,getframe(gcf)); 
        
        hours
    end
    
end
vidObj.close();




