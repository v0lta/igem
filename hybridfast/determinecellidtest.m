function cellid=determinecellidtest(x,y,XLength,rsearch)
	%determine cell id based on coordinates

	%Kx
	Kx=floor(XLength/rsearch)+1;

	%grid x coordinates
	xi=x-mod(x,rsearch);
	yi=y-mod(y,rsearch);

	%indices
	i=xi/rsearch+1;
	j=yi/rsearch+1;

	%cellid
	cellid=i+(j-1)*Kx;
end

