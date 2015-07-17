

%A continuous model of the Keller-Segel system used within the hybrid 
%Model.

%% Load constants.
Da = 1/20;
Db = 1/10;
Dr = 200;
Dh = 200;

dt = 0.0001;
tend = 0.08;
Xlen = 1;
J = 30;

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

constK = -2;
kR = 0.1;
kH = 0.2;
idx = 2:1:(J-1);

%% Set intial conditions.
[x,y] = meshgrid(linspace(0,Xlen,J));
B = 15*(x - x.^2).*(y-y.^2).*exp(-50 .*((x - 0.5).^2 + (y - 0.5).^2 ));
A = 15*(x - x.^2).*(y-y.^2).*exp(-50 .*((x - 0.5).^2 + (y - 0.5).^2 ));
R(:,:) = 0.00001;
H(:,:) = 0;

%%set up movie obj.
vidObj=VideoWriter('simulation_contAHLSel.avi');
set(vidObj,'FrameRate',24);
open(vidObj);


%% run simulation.
for t = 1:1:steps
   % Compute the inner grid values.
   %compute values related to cell A.
   Astep(idx,idx) = A(idx,idx) + (dt./dx^2) .* Da .* ( A(idx+1,idx) + ...
                    A(idx-1,idx) + A(idx,idx+1) + A(idx,idx-1)- 4*A(idx,idx));
                    
   %compute values related to cell B.
   xN = B .* 50.* (constK.*H./R);
   betaN = (xN(idx,idx+1) + xN(idx,idx))./2;
   betaS = (xN(idx,idx-1) + xN(idx,idx))./2;
   betaE = (xN(idx+1,idx) + xN(idx,idx))./2;
   betaW = (xN(idx-1,idx) + xN(idx,idx))./2;

   Bstep(idx,idx) = B(idx,idx) + (dt./dx^2) .* ( Db .* ( B(idx+1,idx) + ...
                    B(idx-1,idx) + B(idx,idx+1) + B(idx,idx-1)- 4*B(idx,idx)) ...
                    -(+ betaE.*(R(idx+1,idx) - R(idx,idx))  ...
                    - betaW.*(R(idx,idx)   - R(idx-1,idx))  ...
                    + betaN.*(R(idx,idx+1) - R(idx,idx))    ... 
                    - betaS.*(R(idx,idx)   - R(idx,idx-1))));  
                
   %Repellant related values. 
   Rstep(idx,idx) = R(idx,idx) + dt .* ( Dr .* (R(idx+1,idx) + R(idx-1,idx) ...
                    + R(idx,idx+1) + R(idx,idx-1) - 4.*R(idx,idx)) + kR.*A(idx,idx));
                
   %AHL related things.
   Hstep(idx,idx) = H(idx,idx) + dt .* ( Dh .* (H(idx+1,idx) + H(idx-1,idx) ...
                    + H(idx,idx+1) + H(idx,idx-1) - 4.*H(idx,idx)) + kH.*A(idx,idx));
   
    
    %% Do the boundary (i==1)
    Astep(1,idx) = A(1,idx) + (dt./dx^2) .* ( Da .* ( A(2,idx) + ...
                   A(1,idx) + A(1,idx+1) + A(1,idx-1)- 4*A(1,idx)));
    
    betaN = (xN(1,idx+1) + xN(1,idx))./2;
    betaS = (xN(1,idx-1) + xN(1,idx))./2;
    betaE = (xN(2,idx) + xN(1,idx))./2;
    betaW = (xN(J,idx) + xN(1,idx))./2;
    Bstep(1,idx) = B(1,idx) + (dt./dx^2) .* ( Db .* ( B(2,idx) + ...
                   B(1,idx) + B(1,idx+1) + B(1,idx-1)- 4*B(1,idx)) ...
                   -(+ betaE.*(R(2,idx) - R(1,idx))  ...
                   - betaW.*(R(1,idx)   - R(J,idx))  ...
                   + betaN.*(R(1,idx+1) - R(1,idx))    ... 
                   - betaS.*(R(1,idx)   - R(1,idx-1))));  
    
     Rstep(1,idx) = R(1,idx) + dt .* ( Dr .* (R(2,idx) + R(J,idx) ...
                    + R(1,idx+1) + R(1,idx-1) - 4.*R(1,idx)) + kR.*A(1,idx));
                
     Hstep(1,idx) = H(1,idx) + dt .* ( Dh .* (H(2,idx) + H(J,idx) ...
                    + H(1,idx+1) + H(1,idx-1) - 4.*H(1,idx)) + kH.*A(1,idx));           
                
    %% Take care of (i==J)
    Astep(J,idx) = A(J,idx) + (dt./dx^2) .* ( Da.* ( B(1,idx) + ...
               A(J-1,idx) + A(J,idx+1) + A(J,idx-1)- 4*B(J,idx)));
        
    betaN = (xN(J,idx+1) + xN(J,idx))./2;
    betaS = (xN(J,idx+1) + xN(J,idx))./2;
    betaE = (xN(1,idx)   + xN(J,idx))./2;
    betaW = (xN(J-1,idx) + xN(J,idx))./2;
    Bstep(J,idx) = B(J,idx) + (dt./dx^2) .* ( Db .* ( B(1,idx) + ...
                     B(J-1,idx) + B(J,idx+1) + B(J,idx-1)- 4*B(J,idx)) ...
                     -(+ betaE.*(R(1,idx) - R(J,idx))  ...
                     - betaW.*(R(J,idx)   - R(J-1,idx))  ...
                     + betaN.*(R(J,idx+1) - R(J,idx))    ... 
                     - betaS.*(R(J,idx)   - R(J,idx-1))));  
                
    Rstep(J,idx) = R(J,idx) + dt .* ( Dr .* (R(1,idx) + R(J-1,idx) ...
                     + R(J,idx+1) + R(J,idx-1) - 4.*R(J,idx)) + kR.*A(J,idx));
                 
    Hstep(J,idx) = H(J,idx) + dt .* ( Dh .* (H(1,idx) + H(J-1,idx) ...
                     + H(J,idx+1) + H(J,idx-1) - 4.*H(J,idx)) + kH.*A(J,idx));             
    
    
    %% Now (j==1)   
    Astep(idx,1) = A(idx,1) + (dt./dx^2) .* ( Da .* ( A(idx+1,1) + ...
                   A(idx-1,1) + A(idx,2) + A(idx,J)- 4*A(idx,1)));
    
    betaN = (xN(idx,2) + xN(idx,1))./2;
    betaS = (xN(idx,J) + xN(idx,1))./2;
    betaE = (xN(idx+1,1) + xN(idx,1))./2;
    betaW = (xN(idx-1,1) + xN(idx,1))./2;
    Bstep(idx,1) = B(idx,1) + (dt./dx^2) .* ( Db .* ( B(idx+1,1) + ...
                     B(idx-1,1) + B(idx,2) + B(idx,J)- 4*B(idx,1)) ...
                     -(+ betaE.*(R(idx+1,1) - R(idx,1))  ...
                     - betaW.*(R(idx,1)     - R(idx-1,1))  ...
                     + betaN.*(R(idx,2)     - R(idx,1))    ... 
                     - betaS.*(R(idx,1)     - R(idx,J))));  
                
    Rstep(idx,1) = R(idx,1) + dt .* ( Dr .* (R(idx+1,1) + R(idx-1,1) ...
                     + R(idx,2) + R(idx,J) - 4.*R(idx,1)) + kR.*A(idx,1));
                 
    Hstep(idx,1) = H(idx,1) + dt .* ( Dh .* (H(idx+1,1) + H(idx-1,1) ...
                     + H(idx,2) + H(idx,J) - 4.*H(idx,1)) + kH.*A(idx,1));
                 
   %% (j==J)
   Astep(idx,J) = A(idx,J) + (dt./dx^2) .* ( Da .* ( A(idx+1,J) + ...
                    A(idx-1,J) + A(idx,1) + A(idx,J-1)- 4*A(idx,J))); 
   
   betaN = (xN(idx,1) + xN(idx,J))./2;
   betaS = (xN(idx,J-1) + xN(idx,J))./2;
   betaE = (xN(idx+1,J) + xN(idx,J))./2;
   betaW = (xN(idx-1,J) + xN(idx,J))./2;

   Bstep(idx,J) =   B(idx,J) + (dt./dx^2) .* ( Db .* ( B(idx+1,J) + ...
                    B(idx-1,J) + B(idx,1) + B(idx,J-1)- 4*B(idx,J)) ...
                    -(+ betaE.*(R(idx+1,J) - R(idx,J))  ...
                    - betaW.*(R(idx,J)   - R(idx-1,J))  ...
                    + betaN.*(R(idx,1) - R(idx,J))    ... 
                    - betaS.*(R(idx,J)   - R(idx,J-1))));  
                
   Rstep(idx,J) = R(idx,J) + dt .* ( Dr .* (R(idx+1,J) + R(idx-1,J) ...
                    + R(idx,1) + R(idx,J-1) - 4.*R(idx,J)) + kR.*A(idx,J));
                
   Hstep(idx,J) = H(idx,J) + dt .* ( Dh .* (H(idx+1,J) + H(idx-1,J) ...
                    + H(idx,1) + H(idx,J-1) - 4.*H(idx,J)) + kH.*A(idx,J));             
   
    
   %% (i==j==1)
   Astep(1,1) = A(1,1) + (dt./dx^2) .* ( Da .* ( A(2,1) + ...
                    A(J,1) + A(1,2) + A(1,J)- 4*A(1,1)));
   
   betaN = (xN(1,2) + xN(1,1))./2;
   betaS = (xN(1,J) + xN(1,1))./2;
   betaE = (xN(2,1) + xN(1,1))./2;
   betaW = (xN(J,1) + xN(1,1))./2;
   Bstep(1,1) = B(1,1) + (dt./dx^2) .* ( Db .* ( B(2,1) + ...
                    B(J,1) + B(1,2) + B(1,J)- 4*B(1,1)) ...
                    -(+ betaE.*(R(2,1) - R(1,1))  ...
                    - betaW.*(R(1,1)   - R(J,1))  ...
                    + betaN.*(R(1,2) - R(1,1))    ... 
                    - betaS.*(R(1,1)   - R(1,J))));                  
   Rstep(1,1) = R(1,1) + dt .* ( Dr .* (R(2,1) + R(J,1) ...
                    + R(1,2) + R(1,J) - 4.*R(1,1)) + kR.*A(1,1));
                
   Hstep(1,1) = H(1,1) + dt .* ( Dh .* (H(2,1) + H(J,1) ...
                    + H(1,2) + H(1,J) - 4.*H(1,1)) + kH.*A(1,1));
                
   %% i=1,j==J
   Astep(1,J) = A(1,J) + (dt./dx^2) .* ( Da .* ( A(2,J) + ...
                    A(J,J) + A(1,1) + A(1,J-1)- 4*A(1,J)));

   betaN = (xN(1,1) + xN(1,J))./2;
   betaS = (xN(1,J-1) + xN(1,J))./2;
   betaE = (xN(2,J) + xN(1,J))./2;
   betaW = (xN(J,J) + xN(1,J))./2;
   Bstep(1,J) = B(1,J) + (dt./dx^2) .* ( Db .* ( B(2,J) + ...
                    B(J,J) + B(1,1) + B(1,J-1)- 4*B(1,J)) ...
                    -(+ betaE.*(R(2,J) - R(1,J))  ...
                    - betaW.*(R(1,J)   - R(J,J))  ...
                    + betaN.*(R(1,1) - R(1,J))    ... 
                    - betaS.*(R(1,J)   - R(1,J-1))));                  
   Rstep(1,J) = R(1,J) + dt .* ( Dr .* (R(2,J) + R(J,J) ...
                    + R(1,1) + R(1,J-1) - 4.*R(1,J)) + kR.*A(1,J));
   Hstep(1,J) = H(1,J) + dt .* ( Dh .* (H(2,J) + H(J,J) ...
                    + H(1,1) + H(1,J-1) - 4.*H(1,J)) + kH.*A(1,J));
   
  %% i==J j==1
   Astep(J,1) = A(J,1) + (dt./dx^2) .* ( Da .* ( A(1,1) + ...
                    A(J-1,1) + A(J,2) + A(J,J)- 4*A(J,1)));   

   betaN = (xN(J,2) + xN(J,1))./2;
   betaS = (xN(J,J) + xN(J,1))./2;
   betaE = (xN(1,1) + xN(J,1))./2;
   betaW = (xN(J-1,1) + xN(J,1))./2;

   Bstep(J,1) = B(J,1) + (dt./dx^2) .* ( Db .* ( B(1,1) + ...
                    B(J-1,1) + B(J,2) + B(J,J)- 4*B(J,1)) ...
                    -(+ betaE.*(R(1,1) - R(J,1))  ...
                    - betaW.*(R(J,1)   - R(J-1,1))  ...
                    + betaN.*(R(J,2)   - R(J,1))    ... 
                    - betaS.*(R(J,1)   - R(J,J))));  
                
   Rstep(J,1) = R(J,1) + dt .* ( Dr .* (R(1,1) + R(J-1,1) ...
                    + R(J,2) + R(J,J) - 4.*R(J,1)) + kR.*A(J,1));
   Hstep(J,1) = H(J,1) + dt .* ( Dh .* (H(1,1) + H(J-1,1) ...
                    + H(J,2) + H(J,J) - 4.*H(J,1)) + kH.*A(J,1));
                
                
   %% i==J, j==J
   Astep(J,J) = A(J,J) + (dt./dx^2) .* ( Da .* ( A(1,J) + ...
                A(J-1,J) + A(J,1) + A(J,J-1) - 4*A(J,J)));
   
   betaN = (xN(J,1) + xN(J,J))./2;
   betaS = (xN(J,J-1) + xN(J,J))./2;
   betaE = (xN(1,J) + xN(J,J))./2;
   betaW = (xN(J-1,J) + xN(J,J))./2;

   Bstep(J,J) = B(J,J) + (dt./dx^2) .* ( Db .* ( B(1,J) + ...
                    B(J-1,J) + B(J,1) + B(J,J-1)- 4*B(J,J)) ...
                    -(+ betaE.*(R(1,J) - R(J,J))    ...
                    - betaW.*(R(J,J)   - R(J-1,J))  ...
                    + betaN.*(R(J,1)   - R(J,J))    ... 
                    - betaS.*(R(J,J)   - R(J,J-1))));  
                
   Rstep(J,J) = R(J,J) + dt .* ( Dr .* (R(1,J) + R(J-1,J) ...
                    + R(J,1) + R(J,J-1) - 4.*R(J,J)) + kR.*A(J,J));
   
   Hstep(J,J) = H(J,J) + dt .* ( Dh .* (H(1,J) + H(J-1,J) ...
                    + H(J,1) + H(J,J-1) - 4.*H(J,J)) + kH.*A(J,J));
    
                
                
    %% Update the arrays!            
    R = Rstep;
    B = Bstep;
    R = Rstep;
    H = Hstep;
    
    figure(1);clf;
    subplot(1,4,1);
    plot(A(J/2,:));
    hold on;
    plot(B(J/2,:));
    plot(R(J/2,:)*100);
    plot(H(J/2,:)*100);
    title('parameters in 2d at x=J/2');
    legend('A','B','rep','AHL')
    xlabel('y');ylabel('#cells and concentration*100');
    subplot(1,4,2);
    %surf(A);shading('flat');xlabel('x');ylabel('y');title('density type A');
    %subplot(1,5,3);
    surf(B);shading('flat');xlabel('x');ylabel('y');title('density type B');
    zlim([-0.01 1]);
    subplot(1,4,3);
    surf(R);shading('flat');xlabel('x');ylabel('y');title('repellant');
    %zlim([0.075 0.1]);
    subplot(1,4,4);
    surf(H);shading('flat');xlabel('x');ylabel('y');title('AHL');
    writeVideo(vidObj,getframe(gcf)); 

    ((steps-t)/steps*100)

    
end
vidObj.close();




