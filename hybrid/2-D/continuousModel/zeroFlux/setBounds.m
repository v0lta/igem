function [ ip1,im1,jp1,jm1 ] = setBounds( i,j,J )
%Zero flux boundary conditions nabla f(x) = 0. 

 % the center of the domain.
 if ((length(i) > 1) && (length(j) > 1))
    ip1 = i+1;
    im1 = i-1;
    jp1 = j+1;
    jm1 = j-1;
 % check if the i index is at a boundary?   
 elseif((length(i) == 1) && (length(j) > 1))
 % the left boundary   
    if (i == 1) 
        ip1 = i+1;
        im1 = i+1;
        jp1 = j+1 ;
        jm1 = j-1;
    elseif (i == J)
 % the right boundary.
        ip1 = J-1;
        im1 = J-1;
        jp1 = j+1 ;
        jm1 = j-1;
    end
 % check if the j index is @ a boundary.
 elseif((length(i) > 1) && (length(j) == 1))   
 % the bottom boundary.
    if (j == 1)
        ip1 = i+1;
        im1 = i-1;
        jp1 = j+1;
        jm1 = j+1;
 % the top boundary.
    elseif (j == J)
        ip1 = i+1;
        im1 = i-1;
        jp1 = J-1;
        jm1 = J-1;   
    end
 % possibly on of the corners.
 elseif (length(i) == 1) && (length(j) == 1)
     %left 
     if (i == 1)
        ip1 = i+1;
        im1 = i+1;
     end
     %right
     if (i == J)
        ip1 = J-1;
        im1 = J-1;
     end
     %bottom
     if (j == 1)
        jp1 = j+1;
        jm1 = j+1;
     end
     %top
     if (j == J)
        jp1 = J-1;
        jm1 = J-1;
     end  
 end
end

