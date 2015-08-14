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

		function newDirection=turn(obj,direction,lambda0A,kappaA,speedA,AHL,dAHL,timestep)
		dAHLx=dAHL(1);
		dAHLy=dAHL(2);

		%lambda0A
		%kappaA
		%AHL
		%currentSpeedA
		%direction
		%dAHLx
		%dAHLy
		lambda=lambda0A-kappaA/(2*AHL)*speedA*(cos(direction)*dAHLx+sin(direction)*dAHLy);

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

		function update(obj,AHLField,lambda0A,speedA,kappaA,VthA,dt)
		%only random diffusion based on diffusion constant
		%domainx=obj.domain(:,1);
		%domainy=obj.domain(:,2);
		domain=obj.domain;
		
		for bacterium=obj.bacteria
			x=bacterium.getxcoordinate();
			y=bacterium.getycoordinate();
			direction=bacterium.getdirection();
			
			AHL=AHLField.interpolconc([x y]);
			dAHL=AHLField.interpolgrad([x y]);

			if AHL>VthA	%AHL above threshold => low diffusion
				%variable speed, fixed turning frequency
				%currentSpeedA=speedA(1);
				%currentLambda0A=lambda0A;
				%variable turning frequency, fixed speed
				currentSpeedA=speedA;
				currentLambda0A=lambda0A.high;
			else		%AHL below threshold => high diffusion
				%variable speed, fixed turning frequency
				%currentSpeedA=speedA(2);
				%currentLambda0A=lambda0A;
				%variable turning frequency, fixed speed
				currentSpeedA=speedA;
				currentLambda0A=lambda0A.high;
			end

			newDirection=obj.turn(direction,currentLambda0A,kappaA,currentSpeedA,AHL,dAHL,dt);%new direction
			%calculate new position
			%new x
			xNew=x+currentSpeedA*cos(newDirection)*dt;
			%new y
			yNew=y+currentSpeedA*sin(newDirection)*dt;

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

			%set new direction
			bacterium.setdirection(newDirection);

			%set new position
			bacterium.setxcoordinate(xNew);
			bacterium.setycoordinate(yNew);
		end
		end
	end
end
