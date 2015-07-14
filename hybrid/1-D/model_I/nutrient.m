classdef nutrient<handle
	properties
		domain;
		concentration;
		gradient;
		leftBoundary;
		rightBoundary;
	end

	methods
		function obj=nutrient(domain,concentration,boundaries)
		obj.domain=domain;
		obj.concentration=concentration;
		obj.calculategradient();
		obj.leftBoundary=boundaries(1);
		obj.rightBoundary=boundaries(2);
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

		%xCoordinate
		%xmin=min(obj.domain)
		%xmax=max(obj.domain)

		x=obj.domain;
		C=field;
		%length(C)
		%C(1001)

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

		function update(obj,bacteriaDensity,consumptionRate,Ds,dt)
		%Update concentration based on bacteria density, consumption rate diffusion constant of nutrient and timestep
		%disp('length concentration');
		%disp(num2str(length(obj.concentration)));
		%disp('length domain');
		%disp(num2str(length(obj.domain)));
		n=length(obj.domain);
		dx=obj.domain(2)-obj.domain(1);
		r=dt/dx^2;

		%bacteria density
		rho=bacteriaDensity(obj.domain);
		%disp('length rho');
		%disp(num2str(length(rho)));

		%matrix
		A=diag(ones(n,1)*(1+2*Ds*r))+diag(ones(n-1,1)*(-Ds*r),1)+diag(ones(n-1,1)*(-Ds*r),-1);
		b=obj.concentration;
		b(1)=b(1)+Ds*r*obj.leftBoundary;
		b(n)=b(n)+Ds*r*obj.rightBoundary;
		%disp('length b');
		%disp(num2str(length(b)));
		%disp('length rho*consumptionRate*dt');
		%disp(num2str(length(rho*consumptionRate*dt)));
		b=b-rho*consumptionRate*dt;
		%[k,l]=size(b);
		%disp(['size of b: ' num2str(k) ' x ' num2str(l)]);
		%[k,l]=size(A);
		%disp(['size of A: ' num2str(k) ' x ' num2str(l)]);

		%Calculate new concentration field
		%disp('A');
		%A;
		%disp('b');
		%b;
		obj.concentration=(A\b')';

		%update gradient
		obj.calculategradient();
		end
	end
end
