function res = chemotaxis(B,Kc,H,R,dx,dt,i,j,J)

 %set the boundaries. 
 [ ip1,im1,jp1,jm1 ] = setBounds( i,j,J );
 
 Xn = -B .* Kc .* H./R;
 betaN = (Xn(i,jp1) + Xn(i,j))./2;
 betaS = (Xn(i,jm1) + Xn(i,j))./2;
 betaE = (Xn(ip1,j) + Xn(i,j))./2;
 betaW = (Xn(im1,j) + Xn(i,j))./2;

 res = dt./(dx^2) .*( -(+ betaE.*(R(ip1,j) - R(i,j))  ...
                        - betaW.*(R(i,j)   - R(im1,j))...
                        + betaN.*(R(i,jp1) - R(i,j))  ... 
                        - betaS.*(R(i,j)   - R(i,jm1))));
end
                    