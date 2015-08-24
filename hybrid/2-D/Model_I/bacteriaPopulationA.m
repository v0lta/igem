classdef bacteriaPopulationA < handle
	properties
		kernelfun;
		bandwidth;
		muA;
		VthA;
		kappaA;

		bacteria;
		domain;
		domainGrid;
	end

	methods
		function obj=bacteriaPopulationA(paramA,bacteriaA,domain,domainGrid)
		%parameters
		obj.kernelfun=paramA.kernelfun;
		obj.bandwidth=paramA.bandwidth;
		obj.muA=paramA.muA;
		obj.VthA=paramA.VthA;
		obj.kappaA=paramA.kappaA;

		%bacteria list and domain
		obj.bacteria=bacteriaA;
		obj.domain=domain;
		obj.domainGrid=domainGrid;
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
			bacterium.setxcoordinate(xNew);
			bacterium.setycoordinate(yNew);
			%bacterium.xCoordinate=xNew;
			%bacterium.yCoordinate=yNew;

			bacteria(i)=bacterium;
		end
		obj.bacteria=bacteria;
		end
	end
end
