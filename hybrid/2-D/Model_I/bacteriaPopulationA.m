classdef bacteriaPopulationA < handle
	properties
		domain;
		domainGrid;
		bacteria;
	end

	methods
		function obj=bacteriaPopulationA(bacteria,domain)
		%initialize bacteria list
		obj.bacteria=bacteria;
		obj.domain=domain;

		[X,Y]=meshgrid(domain(:,1),domain(:,2));
		domainGrid=[];
		domainGrid(:,:,1)=X;
		domainGrid(:,:,2)=Y;
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
		X=obj.domainGrid(:,:,1);
		Y=obj.domainGrid(:,:,2);
		coordinateArray=obj.coordinates();
		rho=KDE2D(coordinateArray,kernelfun,X,Y,bandwidth);
		end

		function update(obj,muA,dt)
		%only random diffusion based on diffusion constant
		domainx=obj.domain(:,1);
		domainy=obj.domain(:,2);
		
		for bacterium=obj.bacteria
			x=bacterium.getxcoordinate();
			y=bacterium.getycoordinate();

			%calculate new position
			%new x
			xNew=x+sqrt(2*muA*dt)*normrnd(0,1);

			%correct for going out of boundary
			if xNew < domainx(1);
				%wall boundary condition
				xNew=domainx(1);
			elseif xNew > domainx(end)
				%wall boundary condition
				xNew=domainx(end);
			end

			%new y
			yNew=y+sqrt(2*muA)*normrnd(0,1);

			%correct for going out of boundary
			if yNew < domainy(1);
				%wall boundary condition
				yNew=domainy(1);
			elseif yNew > domainy(end)
				%wall boundary condition
				yNew=domainy(end);
			end

			%set new position
			bacterium.setxcoordinate(xNew);
			bacterium.setycoordinate(yNew);
		end
		end
	end
end
