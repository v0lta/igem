classdef AHL<handle
	properties
		domain;
		domainGrid;
		concentration;
		gradient;
		westBoundary;
		eastBoundary;
		northBoundary;
		southBoundary;
	end

	methods
		function obj=AHL(domain,concentration,boundaryArray)
		obj.domain=domain;
		obj.concentration=concentration;
		obj.calculategradient();
		obj.westBoundary=boundaryArray(:,1);
		obj.eastBoundary=boundaryArray(:,2);
		obj.northBoundary=boundaryArray(:,3);
		obj.southBoundary=boundaryArray(:,4);

		[X,Y]=meshgrid(domain.x,domain.y);
		domainGrid.X=X;
		domainGrid.Y=Y;
		obj.domainGrid=domainGrid;
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
		interpolatedXDeriv=obj.interpol(obj.gradient.x,coordinateVector);
		interpolatedYDeriv=obj.interpol(obj.gradient.y,coordinateVector);
		interpolatedGradient=[interpolatedXDeriv interpolatedYDeriv];
		end

		function updatedirichlet(obj,rhoAOld,rhoANew,DAHL,alpha,dt)
		%Update concentration based on bacterial density, diffusion constant,
	   	%and consumption rate and timestep
		domain=obj.domain;
		Jx=length(domain.x);
		Jy=length(domain.y);

		dx=domain.x(2)-domain.x(1);
		dy=domain.y(2)-domain.y(1);

		idx=2:Jx-1;
		idy=2:Jy-1;

		muX=DAHL*dt/dx^2;
		muY=DAHL*dt/dy^2;

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
			dt/2*alpha*rhoAOld(idy,idx);

		%disp('RHS1');
		%pause(1);
		%RHS1
		%pause(1);

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

		%disp('B');
		%pause(1);
		%B
		%pause(1);

		%Second half step
		%Calculate RHS 2
		RHS2=(1-muX)*A(idy,idx)+...
			1/2*muX*A(idy,idx-1)+...
			1/2*muX*A(idy,idx+1)+...
			dt/2*alpha*rhoANew(idy,idx);

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

		function updatezeroflux(obj,rhoAOld,rhoANew,DAHL,alpha,k1,dt)
		%Update concentration based on bacterial density, diffusion constant,
	   	%and consumption rate and timestep
		domain=obj.domain;
		Jx=length(domain.x);
		Jy=length(domain.y);

		dx=domain.x(2)-domain.x(1);
		dy=domain.y(2)-domain.y(1);

		idx=2:Jx-1;
		idy=2:Jy-1;

		muX=DAHL*dt/dx^2;
		muY=DAHL*dt/dy^2;

		%west=obj.westBoundary';
		%east=obj.eastBoundary';
		%north=obj.northBoundary';
		%south=obj.southBoundary';

		%ADI zero flux boundary conditions
		%Define matrices
		MRHSx=toeplitz([1-muX 1/2*muX zeros(1,Jx-2)],[1-muX 1/2*muX zeros(1,Jx-2)]);
		MRHSx(1,2)=muX;
		MRHSx(end,end-1)=muX;

		MRHSy=toeplitz([1-muY 1/2*muY zeros(1,Jy-2)],[1-muY 1/2*muY zeros(1,Jy-2)]);
		MRHSy(1,2)=muY;
		MRHSy(end,end-1)=muY;

		MLHSx=toeplitz([1+muX+1/2*k1*dt -1/2*muX zeros(1,Jx-2)],[1+muX+1/2*k1*dt -1/2*muX zeros(1,Jx-2)]);
		MLHSx(1,2)=-muX;
		MLHSx(end,end-1)=-muX;

		MLHSy=toeplitz([1+muY+1/2*k1*dt -1/2*muY zeros(1,Jy-2)],[1+muY+1/2*k1*dt -1/2*muY zeros(1,Jy-2)]);
		MLHSy(1,2)=-muY;
		MLHSy(end,end-1)=-muY;

		%First half step
		%Calculate RHS 1

		A=obj.concentration;
		RHS1=MRHSy*A;
		RHS1=RHS1+dt/2*alpha*rhoAOld;

		%A=obj.concentration;
		%RHS1alt=zeros(size(A));
		%RHS1alt(idy,:)=(1-muY)*A(idy,:)+...
		%			1/2*muY*A(idy-1,:)+...
		%			1/2*muY*A(idy+1,:);

		%RHS1alt(1,:)=(1-muY)*A(1,:)+muY*A(2,:);
		%RHS1alt(end,:)=(1-muY)*A(end,:)+muY*A(end-1,:);
		%RHS1alt=RHS1alt+dt/2*alpha*rhoAOld;

		%RHS1;
		%pause(1);
		%RHS1alt;
		%pause(1);
		%if sum(sum(RHS1==RHS1alt))~=Jx*Jy
		%	disp('RHS matrices do not match!!');
		%end

		%disp('RHS1');
		%pause(1);
		%RHS1
		%pause(1);


		%calculate concentration at half timestep

		%RHS1=RHS1';
		%B=MLHSx\RHS1;
		%A=B';

		A=(MLHSx\RHS1')';
		%B=inv(M)*RHS1;

		%disp('B');
		%pause(1);
		%B
		%pause(1);

		%Second half step
		%Calculate RHS 2

		%RHS2=(1-muX)*A(idy,idx)+...
		%	1/2*muX*A(idy,idx-1)+...
		%	1/2*muX*A(idy,idx+1)+...
		%	dt/2*alpha*rhoANew(idy,idx);

		RHS2=(MRHSx*A')';
		RHS2=RHS2+dt/2*alpha*rhoANew;

		%calculate concentration at full timestep
		A=MLHSy\RHS2;
		%B=inv(M)*RHS2;

		obj.concentration=A;
		end

		function update(obj,rhoAOld,rhoANew,DAHL,alpha,k1,dt)
		%Update concentration based on bacterial density, diffusion constant,
	   	%and consumption rate and timestep

		domain=obj.domain;
		Jx=length(domain.x);
		Jy=length(domain.y);

		%Dirichlet boundary conditions
		obj.updatedirichlet(rhoAOld,rhoANew,DAHL,alpha,dt);

		%Zero flux boundary conditions
		%obj.updatezeroflux(rhoAOld,rhoANew,DAHL,alpha,k1,dt);

		%Correct for negative concentration
		for j=1:Jx
			for i=1:Jy
				if obj.concentration(i,j)<1e-5
					obj.concentration(i,j)=1e-5;
				end
			end
		end

		%update gradient
		obj.calculategradient();
		end
	end
end
