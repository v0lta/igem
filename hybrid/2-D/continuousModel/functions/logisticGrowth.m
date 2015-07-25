function res = logisticGrowth(gamma,As,A,i,j,dt)
 res = dt .* ( gamma.* A(i,j) .* (1 - A(i,j)./As));
end 