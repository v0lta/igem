function output=epanechnikov2DMult(X,Y)

%output=zeros(size(X));
output=epanechnikov(X).*epanechnikov(Y);
