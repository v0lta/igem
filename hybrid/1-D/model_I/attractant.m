classdef attractant<handle
	properties
		domain;
		concentration;
		gradient;
		leftBoundary;
		rightBoundary;
	end

	methods
		function obj=attractant(domain,concentration,boundaries)
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
		%interpolate attractant concentration at bacterium position

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

		function update(obj,bacteriaDensity,Da,eta,alpha,dt)
		%Update concentration based on bacteria density, production rate (eta),
		%diffusion constant of attractant, degradation rate (alpha) and timestep
		n=length(obj.domain);
		dx=obj.domain(2)-obj.domain(1);
		r=dt/dx^2;

		%bacteria density
		%zero flux
		rho=bacteriaDensity(obj.domain(1:end));
		%fixed concentration
		%rho=bacteriaDensity(obj.domain(2:end-1));

		%matrix
		%Zero flux boundary condition
		A=diag(ones(n,1)*(1+2*Da*r+alpha*dt));
		A=A+diag([-2*Da*r;ones(n-2,1)*(-Da*r)],1);
		A=A+diag([ones(n-2,1)*(-Da*r);-2*Da*r],-1);
		b=obj.concentration;

		%fixed concentration boundary condition
		%A=diag(ones(n-2,1)*(1+2*Da*r+alpha*dt));
		%A=A+diag(ones(n-3,1)*(-Da*r),1);
		%A=A+diag(ones(n-3,1)*(-Da*r),-1);
		%b=obj.concentration(2:end-1);
		%b(1)=b(1)+Da*r*obj.leftBoundary;
		%b(end)=b(end)+Da*r*obj.rightBoundary;

		b=b+rho*eta*dt;

		%Calculate new concentration field
		%zero flux
		obj.concentration=(A\b')';
		%fixed concentration
		%obj.concentration=[obj.leftBoundary,(A\b')',obj.rightBoundary];

		%correct for negative concentrations
		for i=1:n
			if obj.concentration(i)<1e-5
				obj.concentration(i)=1e-5;
			end
		end

		%update gradient
		obj.calculategradient();
		end
	end
end
