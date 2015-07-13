function output=rho(x,kernel,pointArray,bandwidth)

output=0;
for point=pointArray
	%output=output+kernel(x-point);
	output=output+1/bandwidth*kernel((x-point)/bandwidth);
end

end
