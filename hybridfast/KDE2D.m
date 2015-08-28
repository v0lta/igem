function rho=KDE2D(coordinateArray,kernelfun,X,Y,bandwidth)

rho=zeros(size(X));
[~,numCoordinate]=size(coordinateArray);
bandwidth2=bandwidth^2;

%parallel version
parfor i=1:numCoordinate
%serial version
%for i=1:numCoordinate
	x=coordinateArray(1,i);
	y=coordinateArray(2,i);

	rho=rho+1/bandwidth2*kernelfun((X-x)/bandwidth,(Y-y)/bandwidth);
end
