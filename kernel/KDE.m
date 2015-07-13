function KDEFunction=KDE(pointArray,kernelfun,bandwidth)

%define kernel function
%d=@uni;
%d=@tri;
%d=@(x) normpdf(x,0,1);
%d=@epanechnikov;

%turn into array function
%dArray=@(x) arrayfun(d,x);

%kernel density estimation
rhofixed=@(x) rho(x,kernelfun,pointArray,bandwidth);
%turn into array function
KDEFunction=@(x) arrayfun(rhofixed,x);

%KDEFunction=@(x) *K_h(x/bandwidth);
