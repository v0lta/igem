function output=epanechnikov2DNormFun(x,y)

a=1.178097049578235;
normal=sqrt(x.^2+y.^2);
output=epanechnikov(normal);
output=output/a;
