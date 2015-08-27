classdef bacteriaPopulationB < handle
	properties
		kernelfun;
		bandwidth;
		muB;
		VthB;
		kappaB;
		r0;	%cell radius
		k1;	%spring constant
		k2;	%spring constant
		gamma;	%friction parameter

		bacteria;
		domain;
		domainGrid;

		currentIsBlack;
		counter;
		modulo;
	end

	methods
		function obj=bacteriaPopulationB(paramB,bacteriaB,domain,domainGrid)
		%parameters
		obj.kernelfun=paramB.kernelfun;
		obj.bandwidth=paramB.bandwidth;
		obj.muB=paramB.muB;
		obj.VthB=paramB.VthB;
		obj.kappaB=paramB.kappaB;
		obj.r0=paramB.r0;
		obj.k1=paramB.k1;
		obj.k2=paramB.k2;
		obj.gamma=paramB.gamma;
		obj.modulo=paramB.modulo;

		%bacteria list and domain
		obj.bacteria=bacteriaB;
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

		function update(obj,AHLField,leucineField,dt)
		%update bacteria position based on AHL and leucine fields, diffusion constant,
		%threshold concentration, chemotactic sensitivity and timestep
		muB=obj.muB;
		VthB=obj.VthB;
		kappaB=obj.kappaB;
		domain=obj.domain;
		
		%parallel version
		bacteria=obj.bacteria;
		n=length(bacteria);
		parfor i=1:n
		%	bacterium=bacteria(i);
		%serial version
		%for bacterium=obj.bacteria
		%for i=1:n
			bacterium=bacteria(i);
			%disp('test');
			x=bacterium.getxcoordinate();
			y=bacterium.getycoordinate();
			
			AHL=AHLField.interpolconc([x y]);

			%disp('class of leucineField');
			%class(leucineField)
			%disp('leucineField');
			%leucineField
			leucine=leucineField.interpolconc([x y]);
			dleucine=leucineField.interpolgrad([x y]);
			dleucinex=dleucine(1);
			dleuciney=dleucine(2);

			if AHL>VthB	%AHL above threshold => high diffusion
				currentMuB=muB.high;
			else		%AHL below threshold => low diffusion
				currentMuB=muB.low;
			end

			%calculate new position
			%repelled by AHL
			%new x
			xNew=x-currentMuB*kappaB/leucine*dleucinex+sqrt(2*currentMuB*dt)*normrnd(0,1);
			%new y
			yNew=y-currentMuB*kappaB/leucine*dleuciney+sqrt(2*currentMuB*dt)*normrnd(0,1);

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
			bacterium.setxcoordinate(xNew);
			bacterium.setycoordinate(yNew);
			%bacterium.xCoordinate=xNew;
			%bacterium.yCoordinate=yNew;

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

		function refreshorupdateneighbors(obj)

		%update neighbors
		if mod(obj.counter,obj.modulo)==0
			obj.refreshneighbors();
		else
			obj.updateneighbors();
			%obj.refreshneighbors();
		end

		obj.counter=obj.counter+1;
		end

		function refreshneighbors(obj)
		n=length(obj.bacteria);
		obj.currentIsBlack=~obj.currentIsBlack;

		for i=1:n
			bacterium=obj.bacteria(i);
			%delete all bacteria, except for bacteria previously added in current timestep
			%saves all bacteria A too
			bacterium.purge(obj.currentIsBlack);

			l=i+1;
			for j=l:n
				otherBacterium=obj.bacteria(j);
				r=obj.calculater(bacterium,otherBacterium);

				%if r<=2*obj.rcut
				if r<=2*obj.r0
					
					otherX=otherBacterium.getxcoordinate();
					otherY=otherBacterium.getycoordinate();

					currentX=bacterium.getxcoordinate();
					currentY=bacterium.getycoordinate();

					bacterium.addneighborb(otherBacterium,otherX,otherY,obj.currentIsBlack,r)
					otherBacterium.addneighborb(bacterium,currentX,currentY,obj.currentIsBlack,r)

				end
			end
		end

		end

		function updateneighbors(obj,bacteriaPopB)
		bacteriaB=obj.bacteria;
		n=length(bacteriaB);

		%bacteriaB=bacteriaPopB.bacteria;
		%m=length(bacteriaB);

		for i=1:n
			bacterium=bacteriaB(i);

			%neighboring cells A
			nbAArray=bacterium.getneighboraarray();

			%newNbArray=[];
			newNbXArray=[];
			newNbYArray=[];
			%newNbIsBlackArray=[];
			newNbRArray=[];

			%loop over neighbors A
			for neighbor=nbAArray
				%calculate and update new variables
				r=obj.calculater(bacterium,neighbor);

				%newNbArray=[newNbArray neighbor];
				newNbXArray=[newNbXArray neighbor.getxcoordinate()];
				newNbYArray=[newNbYArray neighbor.getycoordinate()];
				%newNbIsBlackArray=[newNbIsBlackArray currentIsBlack];
				newNbRArray=[newNbRArray r];
			end

			bacterium.updateneighborsacoordinates(newNbXArray,newNbYArray,newNbRArray);

			%neighboring cells B
			nbBArray=bacterium.getneighborbarray();

			%newNbArray=[];
			newNbXArray=[];
			newNbYArray=[];
			%newNbIsBlackArray=[];
			newNbRArray=[];

			%loop over neighbors B
			for neighbor=nbBArray
				%calculate and update new variables
				r=obj.calculater(bacterium,neighbor);

				%newNbArray=[newNbArray neighbor];
				newNbXArray=[newNbXArray neighbor.getxcoordinate()];
				newNbYArray=[newNbYArray neighbor.getycoordinate()];
				%newNbIsBlackArray=[newNbIsBlackArray currentIsBlack];
				newNbRArray=[newNbRArray r];
			end

			bacterium.updateneighborsbcoordinates(newNbXArray,newNbYArray,newNbRArray);

			bacteriaB(i)=bacterium;
		end
		obj.bacteria=bacteriaB;
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

		function [vx,vy]=cellforce(obj,x1,y1,x2,y2,r)
		%Calculates displacement of bacterium A1 due to bacterium B2

		%if r>=2*obj.r0
		r0=obj.r0;
		%rcut=obj.rcut;

		k1=obj.k1;
		k2=obj.k2;
		%k3=obj.k3;
		if r>=r0*2
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

		if r>=r0
			Fr=k2*(2*r0-r);
		else
			Fr=k1*((k1+k2)/k1*r0-r);
		end

		vx=1/obj.gamma*Fr*ex;
		%vx=obj.gamma*Fr*ex;
		vy=1/obj.gamma*Fr*ey;
		%vy=obj.gamma*Fr*ey;

		end

		function [cellvx,cellvy]=totalcellforce(obj,bacterium)
		%Calculate total force due to neighboring cells
		x=bacterium.getxcoordinate();
		y=bacterium.getycoordinate();

		%force due to neighboring cells A
		[nbAXArray,nbAYArray,nbARArray]=bacterium.getneighborsacoordinates();
		numNeighborsA=length(nbAXArray);
		%[cellvx,cellvy]=cellforce

		cellvxAArray=[];
		cellvyAArray=[];
		for k=1:numNeighborsA
			x2=nbAXArray(k);
			y2=nbAYArray(k);
			r=nbARArray(k);

			[vx,vy]=obj.cellforce(x,y,x2,y2,r);
			cellvxAArray=[cellvxAArray vx];
			cellvyAArray=[cellvyAArray vy];
		end

		%force due to neighboring cells B
		[nbBXArray,nbBYArray,nbBRArray]=bacterium.getneighborsbcoordinates();
		numNeighborsB=length(nbBXArray);
		%[cellvx,cellvy]=cellforce

		cellvxBArray=[];
		cellvyBArray=[];
		for k=1:numNeighborsB
			x2=nbBXArray(k);
			y2=nbBYArray(k);
			r=nbBRArray(k);

			[vx,vy]=obj.cellforce(x,y,x2,y2,r);
			cellvxBArray=[cellvxBArray vx];
			cellvyBArray=[cellvyBArray vy];
		end

		%sum of all neighboring cells
		cellvx=sum(cellvxAArray)+sum(cellvxBArray);
		cellvy=sum(cellvyAArray)+sum(cellvyBArray);
		end

		function updateinteraction(obj,AHLField,leucineField,dt)
		%update bacteria position based on AHL, AHL gradient, diffusion constant,
		%Threshold concentration, chemotactic sensitivity constant and timestep
		muB=obj.muB;
		VthB=obj.VthB;
		kappaB=obj.kappaB;
		domain=obj.domain;

		%neighbors are already updated
		
		bacteria=obj.bacteria;
		newBacteria=[];
		n=length(bacteria);
		%parfor i=1:n
		%	bacterium=bacteria(i);
		%serial version
		%for bacterium=obj.bacteria
		for i=1:n
			bacterium=bacteria(i);
			%disp('test');
			x=bacterium.getxcoordinate();
			y=bacterium.getycoordinate();
			
			AHL=AHLField.interpolconc([x y]);

			leucine=leucineField.interpolconc([x y]);
			dleucine=leucineField.interpolgrad([x y]);
			dleucinex=dleucine(1);
			dleuciney=dleucine(2);

			if AHL>VthB	%AHL above threshold => high diffusion
				currentMuB=muB.high;
			else		%AHL below threshold => low diffusion
				currentMuB=muB.low;
			end

			%calculate displacement due to neighboring cells
			[cellvx,cellvy]=obj.totalcellforce(bacterium);
			celldx=cellvx*dt;
			%celldx=cellvx*dt*currentMuB;
			celldy=cellvy*dt;
			%celldy=cellvy*dt*currentMuB;
			
			%calculate new position
			%new x
			xNew=x;
			xNew=xNew+sqrt(2*currentMuB*dt)*normrnd(0,1);
			%xNew=x-currentMuB*kappaB/leucine*dleucinex;
			xNew=x-currentMuB*kappaB*dleucinex;
			xNew=xNew+celldx;
			%new y
			yNew=y;
			yNew=yNew+sqrt(2*currentMuB*dt)*normrnd(0,1);
			%yNew=y-currentMuB*kappaB/leucine*dleuciney;
			yNew=y-currentMuB*kappaB*dleuciney;
			yNew=yNew+celldy;

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
			bacterium.setxcoordinate(xNew);
			bacterium.setycoordinate(yNew);
			%bacterium.xCoordinate=xNew;
			%bacterium.yCoordinate=yNew;

			%newBacteria(i)=bacterium;
			newBacteria=[newBacteria bacterium];
		end
		obj.bacteria=newBacteria;
		end
	end
end
