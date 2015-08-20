classdef bacteriaPopulationA < handle
	properties
		kernelfun;
		bandwidth;
		muA;
		VthA;
		kappaA;
		r0;	%cell radius
		k;	%spring constant
		gamma;	%friction parameter

		bacteria;
		domain;
		domainGrid;

		currentIsBlack;
		counter;
		modulo;
	end

	methods
		function obj=bacteriaPopulationA(paramA,bacteriaA,domain,domainGrid)
		%parameters
		obj.kernelfun=paramA.kernelfun;
		obj.bandwidth=paramA.bandwidth;
		obj.muA=paramA.muA;
		obj.VthA=paramA.VthA;
		obj.kappaA=paramA.kappaA;
		obj.r0=paramA.r0;
		obj.k=paramA.k;
		obj.gamma=paramA.gamma;
		obj.modulo=paramA.modulo;

		%bacteria list and domain
		obj.bacteria=bacteriaA;
		obj.domain=domain;
		obj.domainGrid=domainGrid;

		%neighbor variables
		obj.currentIsBlack=0;
		obj.counter=0;
		end

		function addBacterium(obj,bacterium)
		obj.bacteria=[obj.bacteria bacterium];
		end

		function coordinateArray=coordinates(obj);
		coordinateArray=[];
		for bacterium=obj.bacteria
			x=bacterium.getxcoordinate();
			y=bacterium.getycoordinate();
			coordinateArray(:,end+1)=[x;y];
		end
		end
	
		function rho=bacteriadensity(obj)
		%return bacteria density array evaluated on grid points of domain
		kernelfun=obj.kernelfun;
		bandwidth=obj.bandwidth;
		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;
		coordinateArray=obj.coordinates();
		rho=KDE2D(coordinateArray,kernelfun,X,Y,bandwidth);
		end

		function update(obj,AHLField,dt)
		%update bacteria position based on AHL, AHL gradient, diffusion constant,
		%Threshold concentration, chemotactic sensitivity constant and timestep
		muA=obj.muA;
		VthA=obj.VthA;
		kappaA=obj.kappaA;
		domain=obj.domain;
		
		%parallel version
		bacteria=obj.bacteria;
		n=length(bacteria);
		%parallel
		for i=1:n
		%	bacterium=bacteria(i);
		%serial version
		%for bacterium=obj.bacteria
		%for i=1:n
			bacterium=bacteria(i);
			%disp('test');
			x=bacterium.getxcoordinate();
			y=bacterium.getycoordinate();
			
			AHL=AHLField.interpolconc([x y]);
			dAHL=AHLField.interpolgrad([x y]);
			dAHLx=dAHL(1);
			dAHLy=dAHL(2);

			if AHL>VthA	%AHL above threshold => low diffusion
				currentMuA=muA.low;
			else		%AHL below threshold => high diffusion
				currentMuA=muA.high;
			end

			%calculate new position
			%new x
			xNew=x+currentMuA*kappaA/AHL*dAHLx+sqrt(2*currentMuA*dt)*normrnd(0,1);
			%new y
			yNew=y+currentMuA*kappaA/AHL*dAHLy+sqrt(2*currentMuA*dt)*normrnd(0,1);

			if AHL==0
				disp('ALARM');
				xNew
				yNew
			end

			%correct for going out of boundary
			%x
			if xNew < domain.x(1);
				%wall boundary condition
				xNew=domain.x(1);
			elseif xNew > domain.x(end)
				%wall boundary condition
				xNew=domain.x(end);
			end

			%y
			if yNew < domain.y(1);
				%wall boundary condition
				yNew=domain.y(1);
			elseif yNew > domain.y(end)
				%wall boundary condition
				yNew=domain.y(end);
			end

			%set new position
			%bacterium.setxcoordinate(xNew);
			%bacterium.setycoordinate(yNew);
			bacterium.xCoordinate=xNew;
			bacterium.yCoordinate=yNew;

			bacteria(i)=bacterium;
		end
		obj.bacteria=bacteria;
		end

		function r=calculater(obj,bacterium1,bacterium2)
		x1=bacterium1.getxcoordinate();
		y1=bacterium1.getycoordinate();

		x2=bacterium2.getxcoordinate();
		y2=bacterium2.getycoordinate();

		r=sqrt((x1-x2)^2+(y1-y2)^2);

		end

		function refreshneighbors(obj)
		n=length(obj.bacteria);
		obj.currentIsBlack=~obj.currentIsBlack;

		for i=1:n
			bacterium=obj.bacteria(i);
			bacterium.purge(obj.currentIsBlack);

			newNbArray=[];
			newNbXArray=[];
			newNbYArray=[];
			newNbIsBlackArray=[];
			newNbRArray=[];

			l=i+1;
			%parallel
			for j=l:n
				otherBacterium=obj.bacteria(j);
				r=obj.calculater(bacterium,otherBacterium);

				if r<=obj.r0
					newNbArray=[newNbArray otherBacterium];
					x=otherBacterium.getxcoordinate();
					newNbXArray=[newNbXArray x];
					y=otherBacterium.getycoordinate();
					newNbYArray=[newNbYArray y];
					newNbIsBlackArray=[newNbIsBlackArray obj.currentIsBlack];
					newNbRArray=[newNbRArray r];
				end
			end

			bacterium.addneighborarray(newNbArray,newNbXArray,newNbYArray,newNbIsBlackArray,newNbRArray);
		end

		end

		function updateneighbors(obj)
		bacteria=obj.bacteria;
		n=length(bacteria);

		%parallel
		for i=1:n
			bacterium=bacteria(i);

			nbArray=bacterium.getneighborarray();

			newNbArray=[];
			newNbXArray=[];
			newNbYArray=[];
			newNbIsBlackArray=[];
			newNbRArray=[];

			for neighbor=nbArray
				r=obj.calculater(bacterium,neighbor);

				%newNbArray=[newNbArray neighbor];
				newNbXArray=[newNbXArray neighbor.getxcoordinate()];
				newNbYArray=[newNbYArray neighbor.getycoordinate()];
				%newNbIsBlackArray=[newNbIsBlackArray currentIsBlack];
				newNbRArray=[newNbRArray r];
			end

			bacterium.updateneighborscoordinates(newNbXArray,newNbYArray,newNbRArray);

			bacteria(i)=bacterium;
		end
		obj.bacteria=bacteria;
		end
		
		function [vx,vy]=cellforce(obj,x1,y1,x2,y2,r)
		%Calculates displacement of bacterium 1 due to bacterium 2

		if r>=obj.r0
			vx=0;
			vy=0;
			return
		end

		if r==0
			theta=rand*2*pi;
			ex=cos(theta);
			ey=sin(theta);
		else
			ex=(x1-x2)/r;
			ey=(y1-y2)/r;
		end

		Fr=obj.k*(obj.r0-r);

		vx=obj.gamma*Fr*ex;
		vy=obj.gamma*Fr*ey;

		end

		function updateraw(obj,dt)
		%update bacteria position based on AHL, AHL gradient, diffusion constant,
		%Threshold concentration, chemotactic sensitivity constant and timestep
		currentMuA=obj.muA;
		VthA=obj.VthA;
		kappaA=obj.kappaA;
		domain=obj.domain;
		
		%parallel version
		bacteria=obj.bacteria;
		n=length(bacteria);
		%parallel
		for i=1:n
		%	bacterium=bacteria(i);
		%serial version
		%for bacterium=obj.bacteria
		%for i=1:n
			bacterium=bacteria(i);
			%disp('test');
			x=bacterium.getxcoordinate();
			y=bacterium.getycoordinate();
			
			%calculate new position
			%new x
			xNew=x+sqrt(2*currentMuA*dt)*normrnd(0,1);
			%new y
			yNew=y+sqrt(2*currentMuA*dt)*normrnd(0,1);

			%correct for going out of boundary
			%x
			if xNew < domain.x(1);
				%wall boundary condition
				xNew=domain.x(1);
			elseif xNew > domain.x(end)
				%wall boundary condition
				xNew=domain.x(end);
			end

			%y
			if yNew < domain.y(1);
				%wall boundary condition
				yNew=domain.y(1);
			elseif yNew > domain.y(end)
				%wall boundary condition
				yNew=domain.y(end);
			end

			%set new position
			%bacterium.setxcoordinate(xNew);
			%bacterium.setycoordinate(yNew);
			bacterium.xCoordinate=xNew;
			bacterium.yCoordinate=yNew;

			bacteria(i)=bacterium;
		end
		obj.bacteria=bacteria;
		end

		function updateinteraction(obj,dt)
		%update bacteria position based on AHL, AHL gradient, diffusion constant,
		%Threshold concentration, chemotactic sensitivity constant and timestep
		currentMuA=obj.muA;
		VthA=obj.VthA;
		kappaA=obj.kappaA;
		domain=obj.domain;

		obj.counter=obj.counter+1;

		%update neighbors
		if mod(obj.counter-1,obj.modulo)==0
			obj.refreshneighbors();
		else
			obj.updateneighbors();
			%obj.refreshneighbors();
		end
		
		%parallel version
		bacteria=obj.bacteria;
		n=length(bacteria);
		%parallel
		for i=1:n
		%	bacterium=bacteria(i);
		%serial version
		%for bacterium=obj.bacteria
		%for i=1:n
			bacterium=bacteria(i);
			%disp('test');
			x=bacterium.getxcoordinate();
			y=bacterium.getycoordinate();

			[nbXArray,nbYArray,nbRArray]=bacterium.getneighborscoordinates();
			numNeighbors=length(nbXArray);
			%[cellvx,cellvy]=cellforce

			cellvxArray=[];
			cellvyArray=[];
			for k=1:numNeighbors
				x2=nbXArray(k);
				y2=nbYArray(k);
				r=nbRArray(k);

				[vx,vy]=obj.cellforce(x,y,x2,y2,r);
				cellvxArray=[cellvxArray vx];
				cellvyArray=[cellvyArray vy];
			end

			%disp('test');
			celldx=sum(cellvxArray*dt*currentMuA);
			celldy=sum(cellvyArray*dt*currentMuA);
			
			%calculate new position
			%new x
			xNew=x;
			xNew=xNew+sqrt(2*currentMuA*dt)*normrnd(0,1);
			xNew=xNew+celldx;
			%new y
			yNew=y;
			yNew=y+sqrt(2*currentMuA*dt)*normrnd(0,1);
			yNew=y+celldy;

			%correct for going out of boundary
			%x
			if xNew < domain.x(1);
				%wall boundary condition
				xNew=domain.x(1);
			elseif xNew > domain.x(end)
				%wall boundary condition
				xNew=domain.x(end);
			end

			%y
			if yNew < domain.y(1);
				%wall boundary condition
				yNew=domain.y(1);
			elseif yNew > domain.y(end)
				%wall boundary condition
				yNew=domain.y(end);
			end

			%set new position
			%bacterium.setxcoordinate(xNew);
			%bacterium.setycoordinate(yNew);
			bacterium.xCoordinate=xNew;
			bacterium.yCoordinate=yNew;

			bacteria(i)=bacterium;
		end
		obj.bacteria=bacteria;
		end
	end
end
