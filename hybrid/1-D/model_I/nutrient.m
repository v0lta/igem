classdef nutrient<handle
	properties
		domain;
		concentration;
		gradient;
	end

	methods
		function obj=nutrient(domain,concentration)
		obj.domain=domain;
		obj.concentration=concentration;
		end
		
		function domain=getdomain(obj)
		domain=obj.domain;
		end
		function concentration=getconcentration(obj)
		concentration=obj.concentration;
		end

		function calculategradient(obj)
		n=length(obj.concentration);
		obj.gradient=obj.concentration*0;
		h=obj.domain(2)-obj.domain(1);

		obj.gradient(1)=(obj.concentration(2)-obj.concentration(1))/h;
		obj.gradient(n)=(obj.concentration(n)-obj.concentration(n-1))/h;

		for i=2:n-1
			obj.gradient(i)=(obj.concentration(i+1)-obj.concentration(i-1))/(2*h);
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

		function update(obj,bacteriaDensity,consumptionRate)
		%Update concentration based on bacteria density and consumption rate
		n=length(obj.domain);
		for i=1:n
			x=obj.domain(i);
			rho=bacteriaDensity(x);
			obj.concentration(i)=obj.concentration(i)-consumptionRate*rho;
			%prevent negative concentration
			if obj.concentration(i)<0.01
				obj.concentration(i)=0.01;
			end
		end

		%update gradient
		obj.calculategradient();
		end
	end
end
