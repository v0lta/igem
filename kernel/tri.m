function output=tri(x)

if x<-1 || x>1
	output=0;
elseif x<=0
	output=1+x;
else
	output=1-x;
end

end
