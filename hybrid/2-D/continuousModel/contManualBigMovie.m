

%A continuous model of the Keller-Segel system used within the hybrid 
%Model.

%% Load constants.
Dn = 1/10;
Ds = 1/10;

dt = 0.0001;
tend = 0.2;
Xlen = 1;
J = 30;

dx = Xlen/J;
steps = tend/dt
mu = dt/dx^2

B = zeros(J,J);
R = zeros(J,J);
preDif = zeros(1,J);

constK = -2;
ks = 0.8;
idx = 2:1:(J-1);

%% Set intial conditions.
[x,y] = meshgrid(linspace(0,Xlen,J));
B = 15*(x - x.^2).*(y-y.^2).*exp(-50 .*((x - 0.5).^2 + (y - 0.5).^2 ));
R(:,:) = 0.01;
Bstep = 0.1.*ones(size(B));
Rstep = 0.1.*ones(size(R));


%%set up movie obj.
vidObj=VideoWriter('simulation_contAHL2.avi');
set(vidObj,'FrameRate',24);
open(vidObj);


%% run simulation.
for t = 1:1:steps
    
   %% Compute the inner grid values.
   xN = B .* (constK./R);
   betaN = (xN(idx,idx+1) + xN(idx,idx))./2;
   betaS = (xN(idx,idx-1) + xN(idx,idx))./2;
   betaE = (xN(idx+1,idx) + xN(idx,idx))./2;
   betaW = (xN(idx-1,idx) + xN(idx,idx))./2;

   Bstep(idx,idx) = B(idx,idx) + (dt./dx^2) .* ( Dn .* ( B(idx+1,idx) + ...
                    B(idx-1,idx) + B(idx,idx+1) + B(idx,idx-1)- 4*B(idx,idx)) ...
                    -(+ betaE.*(R(idx+1,idx) - R(idx,idx))  ...
                    - betaW.*(R(idx,idx)   - R(idx-1,idx))  ...
                    + betaN.*(R(idx,idx+1) - R(idx,idx))    ... 
                    - betaS.*(R(idx,idx)   - R(idx,idx-1))));  
                
   Rstep(idx,idx) = R(idx,idx) + dt .* ( Ds .* (R(idx+1,idx) + R(idx-1,idx) ...
                    + R(idx,idx+1) + R(idx,idx-1) - 4.*R(idx,idx)) + ks.*B(idx,idx));
    
    %% Do the boundary (i==1)
    betaN = (xN(1,idx+1) + xN(1,idx))./2;
    betaS = (xN(1,idx-1) + xN(1,idx))./2;
    betaE = (xN(2,idx) + xN(1,idx))./2;
    betaW = (xN(J,idx) + xN(1,idx))./2;
    Bstep(1,idx) = B(1,idx) + (dt./dx^2) .* ( Dn .* ( B(2,idx) + ...
                   B(1,idx) + B(1,idx+1) + B(1,idx-1)- 4*B(1,idx)) ...
                   -(+ betaE.*(R(2,idx) - R(1,idx))  ...
                   - betaW.*(R(1,idx)   - R(J,idx))  ...
                   + betaN.*(R(1,idx+1) - R(1,idx))    ... 
                   - betaS.*(R(1,idx)   - R(1,idx-1))));  
    
     Rstep(1,idx) = R(1,idx) + dt .* ( Ds .* (R(2,idx) + R(J,idx) ...
                    + R(1,idx+1) + R(1,idx-1) - 4.*R(1,idx)) + ks.*B(1,idx));
                
    %% Take care of (i==J)
    betaN = (xN(J,idx+1) + xN(J,idx))./2;
    betaS = (xN(J,idx+1) + xN(J,idx))./2;
    betaE = (xN(1,idx)   + xN(J,idx))./2;
    betaW = (xN(J-1,idx) + xN(J,idx))./2;
    Bstep(J,idx) = B(J,idx) + (dt./dx^2) .* ( Dn .* ( B(1,idx) + ...
                     B(J-1,idx) + B(J,idx+1) + B(J,idx-1)- 4*B(J,idx)) ...
                     -(+ betaE.*(R(1,idx) - R(J,idx))  ...
                     - betaW.*(R(J,idx)   - R(J-1,idx))  ...
                     + betaN.*(R(J,idx+1) - R(J,idx))    ... 
                     - betaS.*(R(J,idx)   - R(J,idx-1))));  
                
    Rstep(J,idx) = R(J,idx) + dt .* ( Ds .* (R(1,idx) + R(J-1,idx) ...
                     + R(J,idx+1) + R(J,idx-1) - 4.*R(J,idx)) + ks.*B(J,idx));
    
    
    %% Now (j==1)            
    xN = B .* (constK./R);
    betaN = (xN(idx,2) + xN(idx,1))./2;
    betaS = (xN(idx,J) + xN(idx,1))./2;
    betaE = (xN(idx+1,1) + xN(idx,1))./2;
    betaW = (xN(idx-1,1) + xN(idx,1))./2;
    Bstep(idx,1) = B(idx,1) + (dt./dx^2) .* ( Dn .* ( B(idx+1,1) + ...
                     B(idx-1,1) + B(idx,2) + B(idx,J)- 4*B(idx,1)) ...
                     -(+ betaE.*(R(idx+1,1) - R(idx,1))  ...
                     - betaW.*(R(idx,1)     - R(idx-1,1))  ...
                     + betaN.*(R(idx,2)     - R(idx,1))    ... 
                     - betaS.*(R(idx,1)     - R(idx,J))));  
                
    Rstep(idx,1) = R(idx,1) + dt .* ( Ds .* (R(idx+1,1) + R(idx-1,1) ...
                     + R(idx,2) + R(idx,J) - 4.*R(idx,1)) + ks.*B(idx,1));
                 
   %% (j==J)
   xN = B .* (constK./R);
   betaN = (xN(idx,1) + xN(idx,J))./2;
   betaS = (xN(idx,J-1) + xN(idx,J))./2;
   betaE = (xN(idx+1,J) + xN(idx,J))./2;
   betaW = (xN(idx-1,J) + xN(idx,J))./2;

   Bstep(idx,J) =   B(idx,J) + (dt./dx^2) .* ( Dn .* ( B(idx+1,J) + ...
                    B(idx-1,J) + B(idx,1) + B(idx,J-1)- 4*B(idx,J)) ...
                    -(+ betaE.*(R(idx+1,J) - R(idx,J))  ...
                    - betaW.*(R(idx,J)   - R(idx-1,J))  ...
                    + betaN.*(R(idx,1) - R(idx,J))    ... 
                    - betaS.*(R(idx,J)   - R(idx,J-1))));  
                
   Rstep(idx,J) = R(idx,J) + dt .* ( Ds .* (R(idx+1,J) + R(idx-1,J) ...
                    + R(idx,1) + R(idx,J-1) - 4.*R(idx,J)) + ks.*B(idx,J));
   
    
   %% (i==j==1)
   betaN = (xN(1,2) + xN(1,1))./2;
   betaS = (xN(1,J) + xN(1,1))./2;
   betaE = (xN(2,1) + xN(1,1))./2;
   betaW = (xN(J,1) + xN(1,1))./2;
   Bstep(1,1) = B(1,1) + (dt./dx^2) .* ( Dn .* ( B(2,1) + ...
                    B(J,1) + B(1,2) + B(1,J)- 4*B(1,1)) ...
                    -(+ betaE.*(R(2,1) - R(1,1))  ...
                    - betaW.*(R(1,1)   - R(J,1))  ...
                    + betaN.*(R(1,2) - R(1,1))    ... 
                    - betaS.*(R(1,1)   - R(1,J))));                  
   Rstep(1,1) = R(1,1) + dt .* ( Ds .* (R(2,1) + R(J,1) ...
                    + R(1,2) + R(1,J) - 4.*R(1,1)) + ks.*B(1,1));
                
   %% i=1,j==J
   xN = B .* (constK./R);
   betaN = (xN(1,1) + xN(1,J))./2;
   betaS = (xN(1,J-1) + xN(1,J))./2;
   betaE = (xN(2,J) + xN(1,J))./2;
   betaW = (xN(J,J) + xN(1,J))./2;
   Bstep(1,J) = B(1,J) + (dt./dx^2) .* ( Dn .* ( B(2,J) + ...
                    B(J,J) + B(1,1) + B(1,J-1)- 4*B(1,J)) ...
                    -(+ betaE.*(R(2,J) - R(1,J))  ...
                    - betaW.*(R(1,J)   - R(J,J))  ...
                    + betaN.*(R(1,1) - R(1,J))    ... 
                    - betaS.*(R(1,J)   - R(1,J-1))));                  
   Rstep(1,J) = R(1,J) + dt .* ( Ds .* (R(2,J) + R(J,J) ...
                    + R(1,1) + R(1,J-1) - 4.*R(1,J)) + ks.*B(1,J));
   
  %% i==J j==1
   xN = B .* (constK./R);
   betaN = (xN(J,2) + xN(J,1))./2;
   betaS = (xN(J,J) + xN(J,1))./2;
   betaE = (xN(1,1) + xN(J,1))./2;
   betaW = (xN(J-1,1) + xN(J,1))./2;

   Bstep(J,1) = B(J,1) + (dt./dx^2) .* ( Dn .* ( B(1,1) + ...
                    B(J-1,1) + B(J,2) + B(J,J)- 4*B(J,1)) ...
                    -(+ betaE.*(R(1,1) - R(J,1))  ...
                    - betaW.*(R(J,1)   - R(J-1,1))  ...
                    + betaN.*(R(J,2)   - R(J,1))    ... 
                    - betaS.*(R(J,1)   - R(J,J))));  
                
   Rstep(J,1) = R(J,1) + dt .* ( Ds .* (R(1,1) + R(J-1,1) ...
                    + R(J,2) + R(J,J) - 4.*R(J,1)) + ks.*B(J,1));
                
                
   %% i==J, j==J
   xN = B .* (constK./R);
   betaN = (xN(J,1) + xN(J,J))./2;
   betaS = (xN(J,J-1) + xN(J,J))./2;
   betaE = (xN(1,J) + xN(J,J))./2;
   betaW = (xN(J-1,J) + xN(J,J))./2;

   Bstep(J,J) = B(J,J) + (dt./dx^2) .* ( Dn .* ( B(1,J) + ...
                    B(J-1,J) + B(J,1) + B(J,J-1)- 4*B(J,J)) ...
                    -(+ betaE.*(R(1,J) - R(J,J))    ...
                    - betaW.*(R(J,J)   - R(J-1,J))  ...
                    + betaN.*(R(J,1)   - R(J,J))    ... 
                    - betaS.*(R(J,J)   - R(J,J-1))));  
                
   Rstep(J,J) = R(J,J) + dt .* ( Ds .* (R(1,J) + R(J-1,J) ...
                    + R(J,1) + R(J,J-1) - 4.*R(J,J)) + ks.*B(J,J));
    
                
                
    %% Update the arrays!            
    R = Rstep;
    B = Bstep;
    
    %create big arrays to plot dots and ring patterns.
    P1 = [ B B B; B B B; B B B];
    P2 = [ R R R; R R R; R R R];
    
    figure(1);clf;
    subplot(1,2,1)
    surf(P1);shading('flat');xlabel('x');ylabel('y');title('density type B');
    zlim([0 1]);
    subplot(1,2,2)
    surf(P2);shading('flat');xlabel('x');ylabel('y');title('repellent');
    %zlim([0.075 0.1]);
    writeVideo(vidObj,getframe(gcf)); 

    (steps/t*100)

    
end
vidObj.close();




