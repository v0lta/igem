function output=epanechnikov(xArray)

[m,n]=size(xArray);
output=zeros([m,n]);
for j=1:n
	for i=1:m
		x=xArray(i,j);
		if x<-1 || x>1
			output(i,j)=0;
		else
			output(i,j)=3/4*(1-x.^2);
		end
	end
end

end
