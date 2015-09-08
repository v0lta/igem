function [ ip1,im1,jp1,jm1 ] = setBounds( i,j,J )
%Periodic boundary conditions, what disappears at the left reappears on the
%right, same with the top and bottom boundary.

 if ((length(i) > 1) && (length(j) > 1))
    ip1 = i+1;
    im1 = i-1;
    jp1 = j+1;
    jm1 = j-1;
 elseif((length(i) == 1) && (length(j) > 1))   
    if (i == 1) 
        ip1 = i+1;
        im1 = J;
        jp1 = j+1 ;
        jm1 = j-1;
    elseif (i == J)
        ip1 = 1;
        im1 = J-1;
        jp1 = j+1 ;
        jm1 = j-1;
    end
 elseif((length(i) > 1) && (length(j) == 1))   
    if (j == 1)
        ip1 = i+1;
        im1 = i-1;
        jp1 = j+1;
        jm1 = J;
    elseif (j == J)
        ip1 = i+1;
        im1 = i-1;
        jp1 = 1;
        jm1 = J-1;   
    end
 elseif (length(i) == 1) && (length(j) == 1)
     if (i == 1)
        ip1 = i+1;
        im1 = J;
     end
     if (i == J)
        ip1 = 1;
        im1 = J-1;
     end
     if (j == 1)
        jp1 = j+1;
        jm1 = J;
     end
     if (j == J)
        jp1 = 1;
        jm1 = J-1;
     end  
 end
end

