%Model II
classdef nutrient<handle
	properties
		domain;
		concentration;
		gradient;
		laplacian;
	end

	methods
		function obj=nutrient(domain,concentration)
		obj.domain=domain;
		obj.concentration=concentration;
		obj.calculategradient();
		obj.calculatelaplacian();
		end
		
		function domain=getdomain(obj)
		domain=obj.domain;
		end
		function concentration=getconcentration(obj)
		concentration=obj.concentration;
		end

		function calculategradient(obj)
		n=length(obj.concentration);
		obj.gradient=zeros(1,n);
		h=obj.domain(2)-obj.domain(1);

		%obj.gradient(1)=(obj.concentration(2)-obj.concentration(1))/h;
		%obj.gradient(n)=(obj.concentration(n)-obj.concentration(n-1))/h;

		for i=2:n-1
			obj.gradient(i)=(obj.concentration(i+1)-obj.concentration(i-1))/(2*h);
		end
		end

		function calculatelaplacian(obj)
		n=length(obj.concentration);
		obj.laplacian=zeros(1,n);
		h=obj.domain(2)-obj.domain(1);

		%obj.laplacian(1)=(obj.concentration(3)-2*obj.concentration(2)+obj.concentration(1))/(h^2);
		%obj.laplacian(n)=(obj.concentration(n)-2*obj.concentration(n-1)+obj.concentration(n-2))/(h^2);
		obj.laplacian(1)=0;
		obj.laplacian(n)=0;

		for i=2:n-1
			obj.laplacian(i)=(obj.concentration(i+1)-2*obj.concentration(i)+obj.concentration(i-1))/(h^2);
		end
		end

		function interpolatedValue=interpol(obj,field,xCoordinate)
		%interpolate nutrient concentration at bacterium position

		x=obj.domain;
		C=field;

		n=length(x);

		%binary search for adjacent grid points
		kMin=1;
		kMax=n;
		k=floor((kMin+kMax)/2);

		while ~(x(k)<=xCoordinate && xCoordinate<=x(k+1))
			%disp('New iteration');
			%disp('xCoordinate');
			%xCoordinate
			%disp('x(k)');
			%x(k)
			%disp('x(kMin)');
			%x(kMin)
			%disp('x(kMax)');
			%x(kMax)
			%%disp('x(1)');
			%%x(1)
			%pause(1);
			%coordinate smaller than x(k)
			if xCoordinate<x(k)
				kMax=k;
				k=floor((kMin+kMax)/2);
			%coordinate larger than x(k)
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
		interpolatedValue=CLeft+(CRight-CLeft)/(xRight-xLeft)*(xCoordinate-xLeft);
		end

		function interpolatedConcentration=interpolconc(obj,xCoordinate)
		field=obj.concentration;
		interpolatedConcentration=obj.interpol(field,xCoordinate);
		end

		function interpolatedGradient=interpolgrad(obj,xCoordinate)
		field=obj.gradient;
		interpolatedGradient=obj.interpol(field,xCoordinate);
		end

		function update(obj,bacteriaDensity,consumptionRate,Ds,timestep)
		%Update concentration based on bacteria density, consumption rate and diffusion constant
		n=length(obj.domain);

		%obj.concentration(1)=0;
		%obj.concentration(n)=0;
		for i=1:n
			x=obj.domain(i);
			rho=bacteriaDensity(x);
			obj.concentration(i)=obj.concentration(i)+(-consumptionRate*rho+Ds*obj.laplacian(i))*timestep;
			%obj.concentration(i)=obj.concentration(i)+(-consumptionRate*rho+Ds*obj.laplacian(i))*1e-6;
			%prevent negative concentration
			if obj.concentration(i)<1e-5
				obj.concentration(i)=1e-5;
			end
		end

		%update gradient
		obj.calculategradient();
		obj.calculatelaplacian();
		end
	end
end
