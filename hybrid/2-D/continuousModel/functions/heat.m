function res = heat(B,D,dt,dx,i,j,J)

 %set the boundaries:
 [ ip1,im1,jp1,jm1 ] = setBounds( i,j,J );
 res = (dt.* D)./dx^2  .* ( B(ip1,j) + B(im1,j) + B(i,jp1) ...
                            + B(i,jm1) - 4*B(i,j));
end