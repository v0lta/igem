classdef leucine<handle
	properties
		Dleucine;
		beta;
		k2;

		domain;
		domainGrid;
		concentration;
		gradient;
		westBoundary;
		eastBoundary;
		northBoundary;
		southBoundary;

		Jx;
		Jy;
	end

	methods
		function obj=leucine(paramleucine,concentration,boundaryArray,domain,domainGrid)
		%parameters
		obj.Dleucine=paramleucine.Dleucine;
		obj.beta=paramleucine.beta;
		obj.k2=paramleucine.k2;

		%domain, gradient, concentration and boundaries
		obj.domain=domain;
		obj.domainGrid=domainGrid;
		obj.Jx=length(domain.x);
		obj.Jy=length(domain.y);
		obj.concentration=concentration;
		obj.calculategradient();
		obj.westBoundary=boundaryArray(:,1);
		obj.eastBoundary=boundaryArray(:,2);
		obj.northBoundary=boundaryArray(:,3);
		obj.southBoundary=boundaryArray(:,4);
		end

		function domain=getdomain(obj)
		domain=obj.domain;
		end

		function domainGrid=getdomaingrid(obj)
		domainGrid=obj.domainGrid;
		end

		function concentration=getconcentration(obj)
		concentration=obj.concentration;
		end

		function calculategradient(obj)
		%domainx=obj.domain(:,1);
		%domainy=obj.domain(:,2);
		domain=obj.domain;
		concentration=obj.concentration;

		dx=domain.x(2)-domain.x(1);
		dy=domain.y(2)-domain.y(1);

		numX=length(domain.x);
		numY=length(domain.y);

		%x derivative
		idx=2:numX-1;
		%inner values
		xGradientInner=(concentration(:,idx+1)-concentration(:,idx-1))/(2*dx);
		%left boundary
		xGradientLeft=(concentration(:,2)-concentration(:,1))/dx;
		%right boundary
		xGradientRight=(concentration(:,numX)-concentration(:,numX-1))/dx;

		%combine
		xGradient=[xGradientLeft xGradientInner xGradientRight];

		%y derivative
		idy=2:numY-1;
		%inner values
		yGradientInner=(concentration(idy+1,:)-concentration(idy-1,:))/(2*dy);
		%upper boundary
		yGradientUpper=(concentration(2,:)-concentration(1,:))/dy;
		%lower boundary
		yGradientLower=(concentration(numY,:)-concentration(numY-1,:))/dy;

		%combine
		yGradient=[yGradientUpper;yGradientInner;yGradientLower];

		obj.gradient.x=xGradient;
		obj.gradient.y=yGradient;
		end

		function interpolatedValue=interpol(obj,field,coordinateVector)
		domain=obj.domain;
		xCoordinate=coordinateVector(1);
		yCoordinate=coordinateVector(2);
		m=length(domain.y);
		n=length(domain.x);

		dx=domain.x(2)-domain.x(1);
		dy=domain.y(2)-domain.y(1);

		%search for i
		yDiff=yCoordinate-domain.y(1);
		iDiff=yDiff/dy;
		i=floor(iDiff)+1;
		if i==m
			i=m-1;
		end

		%i
		%yCoordinate

		%search for j
		xDiff=xCoordinate-domain.x(1);
		jDiff=xDiff/dx;
		j=floor(jDiff)+1;
		if j==n
			j=n-1;
		end

		%delta x and delta y
		deltaY=yCoordinate-domain.y(i);
		deltaX=xCoordinate-domain.x(j);
		sigmaY=1-deltaY;
		sigmaX=1-deltaX;

		%error checking
		if deltaY>dy*1.01 || deltaX>dx*1.01
			disp('ALARM in interpolation');
			disp('i');
			i
			disp('deltaY');
			deltaY
			disp('dy');
			dy
			disp('yCoordinate');
			yCoordinate
			disp('j');
			j
			disp('deltaX');
			deltaX
			disp('dx');
			dx
			disp('xCoordinate');
			xCoordinate
		end
		
		interpolatedValue=1/(dx*dy)*(sigmaX*sigmaY*field(i,j)+...
			sigmaY*deltaX*field(i,j+1)+sigmaX*deltaY*field(i+1,j)+...
			deltaX*deltaY*field(i+1,j+1));
		end

		function interpolatedConcentration=interpolconc(obj,coordinateVector)
		field=obj.concentration;
		interpolatedConcentration=obj.interpol(field,coordinateVector);
		end

		function interpolatedGradient=interpolgrad(obj,coordinateVector)
		gradient=obj.gradient;
		interpolatedXDeriv=obj.interpol(gradient.x,coordinateVector);
		interpolatedYDeriv=obj.interpol(gradient.y,coordinateVector);
		interpolatedGradient=[interpolatedXDeriv interpolatedYDeriv];
        end

        function updateperiodic(obj,rhoANew,dt)
            %Update concentration based on bacterial density, diffusion
            %constant and consumption rate and timestep.
            L = obj.concentration;
            k2 = obj.k2;
            
            dx = obj.domain.x(2) - obj.domain.x(1);
		    dy = obj.domain.y(2) - obj.domain.y(1);
           
            muX = obj.Dleucine*dt/dx^2;
		    muY = obj.Dleucine*dt/dy^2;
            %% Define step matrices:
             %% Define step matrices:
            xStepRight = toeplitz([(1 - muX + 0.25*k2) 0.5*muX zeros(1,obj.Jx-2)],...
                 [(1 - muX + 0.25*k2) 0.5*muX zeros(1,obj.Jx-2)]);
            yStepRight = toeplitz([(1 - muY + 0.25*k2) 0.5*muY zeros(1,obj.Jy-2)],...
                 [(1 - muY + 0.25*k2 ) 0.5*muY zeros(1,obj.Jy-2)]);
            xStepLeft  = toeplitz([(1 + muX - 0.25*k2) -0.5*muX zeros(1,obj.Jx-2)],...
                 [(1 + muX - 0.25*k2) -0.5*muX zeros(1,obj.Jx-2)]);
            yStepLeft  = toeplitz([(1 + muY - 0.25*k2) -0.5*muY zeros(1,obj.Jy-2)],...
                 [(1 + muY - 0.25*k2) -0.5*muY zeros(1,obj.Jy-2)]);

            %periodic boundary conditions:
            xStepRight(1,obj.Jx) = 0.5*muX;  xStepRight(end,1) = 0.5*muX;
            yStepRight(1,obj.Jy) = 0.5*muY;  yStepRight(end,1) = 0.5*muY;
            xStepLeft(1,obj.Jx)  = -0.5*muX; xStepLeft(end,1) = -0.5*muX;
            yStepLeft(1,obj.Jy)  = -0.5*muY; yStepLeft(end,1) = -0.5*muY;
           
            Hhalf = xStepLeft\(yStepRight*L) + dt/2*obj.beta*rhoANew;
            L = (xStepRight*Hhalf)/yStepLeft + dt/2*obj.beta*rhoANew;
            
            obj.concentration = L;
            
        end  
                
		function updatedirichlet(obj,rhoANew,dt)
		%Update concentration based on bacterial density, diffusion constant,
	   	%and consumption rate and timestep
		Dleucine=obj.Dleucine;
		beta=obj.beta;
		domain=obj.domain;
		Jx=obj.Jx;
		Jy=obj.Jy;

		dx=domain.x(2)-domain.x(1);
		dy=domain.y(2)-domain.y(1);

		idx=2:Jx-1;
		idy=2:Jy-1;

		muX=Dleucine*dt/dx^2;
		muY=Dleucine*dt/dy^2;

		west=obj.westBoundary';
		east=obj.eastBoundary';
		north=obj.northBoundary';
		south=obj.southBoundary';

		%ADI fixed boundary conditions, assuming initial conditions and boundary conditions are compatible
		%First half step
		%Calculate RHS 1
		A=obj.concentration;
		RHS1=(1-muY)*A(idy,idx)+...
			1/2*muY*A(idy-1,idx)+...
			1/2*muY*A(idy+1,idx)+...
			dt/2*beta*rhoANew(idy,idx);

		%add known boundary conditions at t+1/2
		RHS1=RHS1';
		RHS1(1,:)=RHS1(1,:)+1/2*muX*west(idy);
		RHS1(end,:)=RHS1(end,:)+1/2*muX*east(idy);

		%coefficient matrix
		M=toeplitz([1+muX -1/2*muX zeros(1,Jx-4)],[1+muX -1/2*muX zeros(1,Jx-4)]);

		%calculate concentration at half timestep
		B=M\RHS1;
		%B=inv(M)*RHS1;
		A(idy,idx)=B';

		%Second half step
		%Calculate RHS 2
		RHS2=(1-muX)*A(idy,idx)+...
			1/2*muX*A(idy,idx-1)+...
			1/2*muX*A(idy,idx+1)+...
			dt/2*beta*rhoANew(idy,idx);

		%add known boundary conditions at t+1
		RHS2(1,:)=RHS2(1,:)+1/2*muY*south(idx);
		RHS2(end,:)=RHS2(end,:)+1/2*muY*north(idx);

		%coefficient matrix
		M=toeplitz([1+muY -1/2*muY zeros(1,Jy-4)],[1+muX -1/2*muX zeros(1,Jy-4)]);

		%calculate concentration at full timestep
		B=M\RHS2;
		%B=inv(M)*RHS2;

		obj.concentration(idy,idx)=B;
		end

		function updatezeroflux(obj,rhoANew,dt)
		%Update concentration based on bacterial density, diffusion constant,
	   	%and production rate, linear degradation constant and timestep
		Dleucine=obj.Dleucine;
		beta=obj.beta;
		k2=obj.k2;
		domain=obj.domain;
		Jx=obj.Jx;
		Jy=obj.Jy;

		dx=domain.x(2)-domain.x(1);
		dy=domain.y(2)-domain.y(1);

		idx=2:Jx-1;
		idy=2:Jy-1;

		muX=Dleucine*dt/dx^2;
		muY=Dleucine*dt/dy^2;

		%ADI zero flux boundary conditions
		%Define matrices
		MRHSx=toeplitz([1-muX 1/2*muX zeros(1,Jx-2)],[1-muX 1/2*muX zeros(1,Jx-2)]);
		MRHSx(1,2)=muX;
		MRHSx(end,end-1)=muX;

		MRHSy=toeplitz([1-muY 1/2*muY zeros(1,Jy-2)],[1-muY 1/2*muY zeros(1,Jy-2)]);
		MRHSy(1,2)=muY;
		MRHSy(end,end-1)=muY;

		MLHSx=toeplitz([1+muX+1/2*k2*dt -1/2*muX zeros(1,Jx-2)],[1+muX+1/2*k2*dt -1/2*muX zeros(1,Jx-2)]);
		MLHSx(1,2)=-muX;
		MLHSx(end,end-1)=-muX;

		MLHSy=toeplitz([1+muY+1/2*k2*dt -1/2*muY zeros(1,Jy-2)],[1+muY+1/2*k2*dt -1/2*muY zeros(1,Jy-2)]);
		MLHSy(1,2)=-muY;
		MLHSy(end,end-1)=-muY;

		%First half step
		%Calculate RHS 1
		A=obj.concentration;
		RHS1=MRHSy*A;
		RHS1=RHS1+dt/2*beta*rhoANew;

		%calculate concentration at half timestep
		A=(MLHSx\RHS1')';

		%Second half step
		%Calculate RHS 2
		RHS2=(MRHSx*A')';
		RHS2=RHS2+dt/2*beta*rhoANew;

		%calculate concentration at full timestep
		A=MLHSy\RHS2;

		obj.concentration=A;
		end

		function update(obj,rhoANew,dt)
		%Update concentration based on bacterial density, diffusion constant,
	   	%and consumption rate and timestep

		Jx=obj.Jx;
		Jy=obj.Jy;

		%Dirichlet boundary conditions
		%obj.updatedirichlet(rhoANew,dt);

		%Zero flux boundary conditions
		%obj.updatezeroflux(rhoANew,dt);
        
        %Periodic boundary conditions:
        obj.updateperiodic(rhoANew,dt);

		concentration=obj.concentration;
		%Correct for negative concentration
		for j=1:Jx
			for i=1:Jy
				if concentration(i,j)<1e-5
					concentration(i,j)=1e-5;
				end
			end
		end
		obj.concentration=concentration;

		%update gradient
		obj.calculategradient();
		end
	end
end
