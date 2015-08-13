classdef bacteriaPopulationB < handle
	properties
		domain;
		domainGrid;
		bacteria;
	end

	methods
		function obj=bacteriaPopulationB(bacteria,domain)
		%initialize bacteria list
		obj.bacteria=bacteria;
		obj.domain=domain;

		%domainx=domain(:,1);
		%domainy=domain(:,2);
		[X,Y]=meshgrid(domain.x,domain.y);
		%domainGrid=[];
		domainGrid.X=X;
		domainGrid.Y=Y;
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
	
		function rho=bacteriadensity(obj,kernelfun,bandwidth)
		%return bacteria density array evaluated on grid points of domain
		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;
		coordinateArray=obj.coordinates();
		rho=KDE2D(coordinateArray,kernelfun,X,Y,bandwidth);
		end

		function update(obj,AHLField,muB,VthB,kappaB,dt)
		%only random diffusion based on diffusion constant
		%domainx=obj.domain(:,1);
		%domainy=obj.domain(:,2);
		domain=obj.domain;
		
		for bacterium=obj.bacteria
			x=bacterium.getxcoordinate();
			y=bacterium.getycoordinate();
			
			AHL=AHLField.interpolconc([x y]);
			dAHL=AHLField.interpolgrad([x y]);
			dAHLx=dAHL(1);
			dAHLy=dAHL(2);

			if AHL>VthB	%AHL above threshold => high diffusion
				currentMuB=muB.high;
			else		%AHL below threshold => low diffusion
				currentMuB=muB.low;
			end

			%calculate new position
			%repelled by AHL
			%new x
			xNew=x-currentMuB*kappaB/AHL*dAHLx+sqrt(2*currentMuB*dt)*normrnd(0,1);
			%new y
			yNew=y-currentMuB*kappaB/AHL*dAHLy+sqrt(2*currentMuB*dt)*normrnd(0,1);

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
		end
		end
	end
end
