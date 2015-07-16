

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

N = zeros(J,J);
S = zeros(J,J);
preDif = zeros(1,J);

constK = 2;
ks = 0.8;
idx = 2:1:(J-1);

%% Set intial conditions.
[x,y] = meshgrid(linspace(0,Xlen,J));
N = 15*(x - x.^2).*(y-y.^2).*exp(-50 .*((x - 0.5).^2 + (y - 0.5).^2 ));
S(:,:) = 0.1;
Nstep = 0.1.*ones(size(N));
Sstep = 0.1.*ones(size(S));


%%set up movie obj.
vidObj=VideoWriter('simulation_cont2.avi');
set(vidObj,'FrameRate',24);
open(vidObj);


%% run simulation.
for t = 1:1:steps
    
   %% Compute the inner grid values.
   xN = N .* (constK./S);
   betaN = (xN(idx,idx+1) + xN(idx,idx))./2;
   betaS = (xN(idx,idx-1) + xN(idx,idx))./2;
   betaE = (xN(idx+1,idx) + xN(idx,idx))./2;
   betaW = (xN(idx-1,idx) + xN(idx,idx))./2;

   Nstep(idx,idx) = N(idx,idx) + (dt./dx^2) .* ( Dn .* ( N(idx+1,idx) + ...
                    N(idx-1,idx) + N(idx,idx+1) + N(idx,idx-1)- 4*N(idx,idx)) ...
                    -(+ betaE.*(S(idx+1,idx) - S(idx,idx))  ...
                    - betaW.*(S(idx,idx)   - S(idx-1,idx))  ...
                    + betaN.*(S(idx,idx+1) - S(idx,idx))    ... 
                    - betaS.*(S(idx,idx)   - S(idx,idx-1))));  
                
   Sstep(idx,idx) = S(idx,idx) + dt .* ( Ds .* (S(idx+1,idx) + S(idx-1,idx) ...
                    + S(idx,idx+1) + S(idx,idx-1) - 4.*S(idx,idx)) - ks.*N(idx,idx));
    
    %% Do the boundary (i==1)
    betaN = (xN(1,idx+1) + xN(1,idx))./2;
    betaS = (xN(1,idx-1) + xN(1,idx))./2;
    betaE = (xN(2,idx) + xN(1,idx))./2;
    betaW = (xN(J,idx) + xN(1,idx))./2;
    Nstep(1,idx) = N(1,idx) + (dt./dx^2) .* ( Dn .* ( N(2,idx) + ...
                   N(1,idx) + N(1,idx+1) + N(1,idx-1)- 4*N(1,idx)) ...
                   -(+ betaE.*(S(2,idx) - S(1,idx))  ...
                   - betaW.*(S(1,idx)   - S(J,idx))  ...
                   + betaN.*(S(1,idx+1) - S(1,idx))    ... 
                   - betaS.*(S(1,idx)   - S(1,idx-1))));  
    
     Sstep(1,idx) = S(1,idx) + dt .* ( Ds .* (S(2,idx) + S(J,idx) ...
                    + S(1,idx+1) + S(1,idx-1) - 4.*S(1,idx)) - ks.*N(1,idx));
                
    %% Take care of (i==J)
    betaN = (xN(J,idx+1) + xN(J,idx))./2;
    betaS = (xN(J,idx+1) + xN(J,idx))./2;
    betaE = (xN(1,idx)   + xN(J,idx))./2;
    betaW = (xN(J-1,idx) + xN(J,idx))./2;
    Nstep(J,idx) = N(J,idx) + (dt./dx^2) .* ( Dn .* ( N(1,idx) + ...
                     N(J-1,idx) + N(J,idx+1) + N(J,idx-1)- 4*N(J,idx)) ...
                     -(+ betaE.*(S(1,idx) - S(J,idx))  ...
                     - betaW.*(S(J,idx)   - S(J-1,idx))  ...
                     + betaN.*(S(J,idx+1) - S(J,idx))    ... 
                     - betaS.*(S(J,idx)   - S(J,idx-1))));  
                
    Sstep(J,idx) = S(J,idx) + dt .* ( Ds .* (S(1,idx) + S(J-1,idx) ...
                     + S(J,idx+1) + S(J,idx-1) - 4.*S(J,idx)) - ks.*N(J,idx));
    
    
    %% Now (j==1)            
    xN = N .* (constK./S);
    betaN = (xN(idx,2) + xN(idx,1))./2;
    betaS = (xN(idx,J) + xN(idx,1))./2;
    betaE = (xN(idx+1,1) + xN(idx,1))./2;
    betaW = (xN(idx-1,1) + xN(idx,1))./2;
    Nstep(idx,1) = N(idx,1) + (dt./dx^2) .* ( Dn .* ( N(idx+1,1) + ...
                     N(idx-1,1) + N(idx,2) + N(idx,J)- 4*N(idx,1)) ...
                     -(+ betaE.*(S(idx+1,1) - S(idx,1))  ...
                     - betaW.*(S(idx,1)     - S(idx-1,1))  ...
                     + betaN.*(S(idx,2)     - S(idx,1))    ... 
                     - betaS.*(S(idx,1)     - S(idx,J))));  
                
    Sstep(idx,1) = S(idx,1) + dt .* ( Ds .* (S(idx+1,1) + S(idx-1,1) ...
                     + S(idx,2) + S(idx,J) - 4.*S(idx,1)) - ks.*N(idx,1));
                 
   %% (j==J)
   xN = N .* (constK./S);
   betaN = (xN(idx,1) + xN(idx,J))./2;
   betaS = (xN(idx,J-1) + xN(idx,J))./2;
   betaE = (xN(idx+1,J) + xN(idx,J))./2;
   betaW = (xN(idx-1,J) + xN(idx,J))./2;

   Nstep(idx,J) =   N(idx,J) + (dt./dx^2) .* ( Dn .* ( N(idx+1,J) + ...
                    N(idx-1,J) + N(idx,1) + N(idx,J-1)- 4*N(idx,J)) ...
                    -(+ betaE.*(S(idx+1,J) - S(idx,J))  ...
                    - betaW.*(S(idx,J)   - S(idx-1,J))  ...
                    + betaN.*(S(idx,1) - S(idx,J))    ... 
                    - betaS.*(S(idx,J)   - S(idx,J-1))));  
                
   Sstep(idx,J) = S(idx,J) + dt .* ( Ds .* (S(idx+1,J) + S(idx-1,J) ...
                    + S(idx,1) + S(idx,J-1) - 4.*S(idx,J)) - ks.*N(idx,J));
   
    
   %% (i==j==1)
   betaN = (xN(1,2) + xN(1,1))./2;
   betaS = (xN(1,J) + xN(1,1))./2;
   betaE = (xN(2,1) + xN(1,1))./2;
   betaW = (xN(J,1) + xN(1,1))./2;
   Nstep(1,1) = N(1,1) + (dt./dx^2) .* ( Dn .* ( N(2,1) + ...
                    N(J,1) + N(1,2) + N(1,J)- 4*N(1,1)) ...
                    -(+ betaE.*(S(2,1) - S(1,1))  ...
                    - betaW.*(S(1,1)   - S(J,1))  ...
                    + betaN.*(S(1,2) - S(1,1))    ... 
                    - betaS.*(S(1,1)   - S(1,J))));                  
   Sstep(1,1) = S(1,1) + dt .* ( Ds .* (S(2,1) + S(J,1) ...
                    + S(1,2) + S(1,J) - 4.*S(1,1)) - ks.*N(1,1));
                
   %% i=1,j==J
   xN = N .* (constK./S);
   betaN = (xN(1,1) + xN(1,J))./2;
   betaS = (xN(1,J-1) + xN(1,J))./2;
   betaE = (xN(2,J) + xN(1,J))./2;
   betaW = (xN(J,J) + xN(1,J))./2;
   Nstep(1,J) = N(1,J) + (dt./dx^2) .* ( Dn .* ( N(2,J) + ...
                    N(J,J) + N(1,1) + N(1,J-1)- 4*N(1,J)) ...
                    -(+ betaE.*(S(2,J) - S(1,J))  ...
                    - betaW.*(S(1,J)   - S(J,J))  ...
                    + betaN.*(S(1,1) - S(1,J))    ... 
                    - betaS.*(S(1,J)   - S(1,J-1))));                  
   Sstep(1,J) = S(1,J) + dt .* ( Ds .* (S(2,J) + S(J,J) ...
                    + S(1,1) + S(1,J-1) - 4.*S(1,J)) - ks.*N(1,J));
   
  %% i==J j==1
   xN = N .* (constK./S);
   betaN = (xN(J,2) + xN(J,1))./2;
   betaS = (xN(J,J) + xN(J,1))./2;
   betaE = (xN(1,1) + xN(J,1))./2;
   betaW = (xN(J-1,1) + xN(J,1))./2;

   Nstep(J,1) = N(J,1) + (dt./dx^2) .* ( Dn .* ( N(1,1) + ...
                    N(J-1,1) + N(J,2) + N(J,J)- 4*N(J,1)) ...
                    -(+ betaE.*(S(1,1) - S(J,1))  ...
                    - betaW.*(S(J,1)   - S(J-1,1))  ...
                    + betaN.*(S(J,2)   - S(J,1))    ... 
                    - betaS.*(S(J,1)   - S(J,J))));  
                
   Sstep(J,1) = S(J,1) + dt .* ( Ds .* (S(1,1) + S(J-1,1) ...
                    + S(J,2) + S(J,J) - 4.*S(J,1)) - ks.*N(J,1));
                
                
   %% i==J, j==J
   xN = N .* (constK./S);
   betaN = (xN(J,1) + xN(J,J))./2;
   betaS = (xN(J,J-1) + xN(J,J))./2;
   betaE = (xN(1,J) + xN(J,J))./2;
   betaW = (xN(J-1,J) + xN(J,J))./2;

   Nstep(J,J) = N(J,J) + (dt./dx^2) .* ( Dn .* ( N(1,J) + ...
                    N(J-1,J) + N(J,1) + N(J,J-1)- 4*N(J,J)) ...
                    -(+ betaE.*(S(1,J) - S(J,J))    ...
                    - betaW.*(S(J,J)   - S(J-1,J))  ...
                    + betaN.*(S(J,1)   - S(J,J))    ... 
                    - betaS.*(S(J,J)   - S(J,J-1))));  
                
   Sstep(J,J) = S(J,J) + dt .* ( Ds .* (S(1,J) + S(J-1,J) ...
                    + S(J,1) + S(J,J-1) - 4.*S(J,J)) - ks.*N(J,J));
    
                
                
    %% Update the arrays!            
    S = Sstep;
    N = Nstep;
    
    figure(1);clf;
    subplot(1,2,1)
    surf(N);shading('flat');xlabel('x');ylabel('y');title('cell density');
    %ylim([0 1]);
    subplot(1,2,2)
    surf(S);shading('flat');xlabel('x');ylabel('y');title('nutrient');
    %ylim([0.075 0.1]);
    writeVideo(vidObj,getframe(gcf)); 

    (steps/t*100)

    
end
vidObj.close();




