function res = chemotaxis(B,Kc,H,R,dx,dt,i,j,J)
%compute the chemotaxis term of the partial differential equation using a finite
%volume approach. (compare morton meyers numerical solution of partial differential
%equations page 198).


 %set the boundaries. 
 [ ip1,im1,jp1,jm1 ] = setBounds( i,j,J );
 
 %compute the values of the variable coefficient.
 Xn = -B .* Kc .* H./R;
 %compute the half steps for the finite volume scheme.
 betaN = (Xn(i,jp1) + Xn(i,j))./2;
 betaS = (Xn(i,jm1) + Xn(i,j))./2;
 betaE = (Xn(ip1,j) + Xn(i,j))./2;
 betaW = (Xn(im1,j) + Xn(i,j))./2;

%use the finite volume scheme to compute the chemotaxis term.
 res = dt./(dx^2) .*( -(+ betaE.*(R(ip1,j) - R(i,j))  ...
                        - betaW.*(R(i,j)   - R(im1,j))...
                        + betaN.*(R(i,jp1) - R(i,j))  ... 
                        - betaS.*(R(i,j)   - R(i,jm1))));
end
                    
