function output=epanechnikov2DMultFun(x,y)

%output=zeros(size(X));
output=epanechnikov(x).*epanechnikov(y);
