function output=epanechnikov(x)

if x<-1 || x>1
	output=0;
else
	output=3/4*(1-x^2);
end

end
