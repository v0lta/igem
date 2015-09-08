function [ Astep,Bstep,Rstep,Hstep ] = modKS( i,j,A,B,R,H,Da,Db,Dr,Dh,dt, ...
            J,Kc,gamma,As,Bs,kr,kh,KAbsorbR,KAbsorbH,dx)
  %use the heat, chemotaxis and logostic growth functions as well as linear growth and
  %decay terms to assemble the Keller-Segel model.
                           
  Astep =  A(i,j) + heat(A,Da,dt,dx,i,j,J) + ...
                    logisticGrowth(gamma,As,A,i,j,dt);
                    
  %compute values related to cell B.
  Bstep = B(i,j) + heat(B,Db,dt,dx,i,j,J) + ...
                   chemotaxis(B,Kc,H,R,dx,dt,i,j,J) + ...
                   logisticGrowth(gamma,Bs,B,i,j,dt);
                
                    
                
                
  %Repellant related values. 
  Rstep = R(i,j) + heat(R,Dr,dt,dx,i,j,J) ...
                   + dt .* ( kr.*A(i,j)...
                   - KAbsorbR * R(i,j));
                
  %AHL related things.
  Hstep = H(i,j) + heat(H,Dh,dt,dx,i,j,J) ...
                   + dt .* ( kh.*A(i,j)...
                   - KAbsorbH * H(i,j));
end

