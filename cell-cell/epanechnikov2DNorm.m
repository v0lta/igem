function output=epanechnikov2DNorm(X,Y)

a=1.178097049578235;
%output=zeros(size(X));
normArray=sqrt(X.^2+Y.^2);
output=epanechnikov(normArray);
output=output/a;
