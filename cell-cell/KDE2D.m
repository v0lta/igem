function rho=KDE2D(coordinateArray,kernelfun,X,Y,bandwidth)

rho=zeros(size(X));
[m,numCoordinate]=size(coordinateArray);

%parallel version
parfor i=1:numCoordinate
%serial version
%for i=1:numCoordinate
	x=coordinateArray(1,i);
	y=coordinateArray(2,i);

	rho=rho+1/bandwidth*kernelfun((X-x)/bandwidth,(Y-y)/bandwidth);
end
