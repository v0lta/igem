function output=interpol(xCoordinate,nutrientField)

x=nutrientField.domain;
C=nutrientField.concentration;

n=length(x);

%binary search for adjacent grid points
kMin=1;
kMax=n;
k=floor((kMin+kMax)/2);

while ~(x(k)<=xCoordinate && xCoordinate<=x(k+1))
	if xCoordinate<x(k)
		kMax=k;
		k=floor((kMin+kMax)/2);
	else
		kMin=k+1;
		k=floor((kMin+kMax)/2);
	end
end

%xCoordinate found between x(k) and x(k+1)
CLeft=C(k);
CRight=C(k+1);
xLeft=x(k);
xRight=x(k+1);

%calculate interpolated value
output=CLeft+(CRight-CLeft)/(xRight-xLeft)*(xCoordinate-xLeft);
