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

		function newDirection=turn(obj,direction,lambda0B,kappaB,speedB,leucine,dleucine,timestep)
		dleucinex=dleucine(1);
		dleuciney=dleucine(2);

		lambda=lambda0B+kappaB/(2*leucine)*speedB*(cos(direction)*dleucinex+sin(direction)*dleuciney);

		%disp('dS');
		%dS

		%disp('lambda');
		%lambda

		f=rand;
		if f<lambda*timestep	%turn
			newDirection=rand*2*pi;
		else		%don't turn
			newDirection=direction;
		end
		end

		function update(obj,AHLField,leucineField,lambda0B,speedB,kappaB,VthB,dt)
		%update bacterie position based on AHL and leucine fields, diffusion constant,
		%threshold concentration, chemotactic sensitivity and timestep
		%domainx=obj.domain(:,1);
		%domainy=obj.domain(:,2);
		domain=obj.domain;
		
		for bacterium=obj.bacteria
			x=bacterium.getxcoordinate();
			y=bacterium.getycoordinate();
			direction=bacterium.getdirection();%direction
			
			AHL=AHLField.interpolconc([x y]);

			%disp('class of leucineField');
			%class(leucineField)
			%disp('leucineField');
			%leucineField
			leucine=leucineField.interpolconc([x y]);
			dleucine=leucineField.interpolgrad([x y]);
			dleucinex=dleucine(1);
			dleuciney=dleucine(2);

			if AHL>VthB	%AHL is above threshold => high diffusion
				%variable speed, fixed turning frequency
				%currentSpeedB=speedB(2);
				%currentLambda0B=lambda0B;
				%variable turning frequency, fixed speed
				currentSpeedB=speedB;
				currentLambda0B=lambda0B.low;
			else		%AHL is below threshold => low diffusion
				%variable speed, fixed turning frequency
				%currentSpeedB=speedB(1);
				%currentLambda0B=lambda0B;
				%variable turning frequency, fixed speed
				currentSpeedB=speedB;
				currentLambda0B=lambda0B.high;
			end

			newDirection=obj.turn(direction,currentLambda0B,kappaB,currentSpeedB,leucine,dleucine,dt);%new direction

			%calculate new position
			%repelled by leucine
			%new x
			xNew=x+currentSpeedB*cos(newDirection)*dt;
			%new y
			yNew=y+currentSpeedB*sin(newDirection)*dt;

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

			%set new direction
			bacterium.setdirection(newDirection);

			%set new position
			bacterium.setxcoordinate(xNew);
			bacterium.setycoordinate(yNew);
		end
		end
	end
end
