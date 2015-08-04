%Model II
classdef bacteriaPopulationA < handle
	properties
		domain;
		bacteria;
	end

	methods
		function obj=bacteriaPopulationA(bacteria,domain)
		%initialize bacteria list
		obj.bacteria=bacteria;
		obj.domain=domain;
		end

		function addBacterium(obj,bacterium)
		obj.bacteria=[obj.bacteria bacterium];
		end

		function coordinateArray=coordinates(obj);
		coordinateArray=[];
		for bacterium=obj.bacteria
			coordinateArray(end+1)=bacterium.getxcoordinate();
		end
		end

		function densityfun=bacteriadensity(obj,kernelfun,bandwidth)
		%Return continuous density function of bacteria population

		coordinateArray=obj.coordinates();
		densityfun=KDE(coordinateArray,kernelfun,bandwidth);
		end

		function newDirection=turn(obj,direction,lambda0A,kappaA,currentSpeedA,AHL,dAHL,timestep)
		foo=kappaA*currentSpeedA*dAHL/(2*AHL);
		%disp('dS');
		%dS

		%if direction==0 && dS >0
		%	disp('against the gradient!');
		%end

		if direction==0	%moving to the left
			lambda=lambda0A+foo;
		else			%moving to the right
			lambda=lambda0A-foo;
		end

		%disp('lambda');
		%lambda

		f=rand;
		if f<lambda*timestep	%turn
			newDirection=~direction;
		else		%don't turn
			newDirection=direction;
		end
		end

		function update(obj,AHLField,lambda0A,speedA,kappaA,VthA,timestep)
		%Update bacteria position based on AHL field, base turning frequency, constant speed,
		%threshold AHL concentration and chemotactic sensitivity constant
		for bacterium=obj.bacteria
			x=bacterium.getxcoordinate();%x coordinate
			direction=bacterium.getdirection();%direction
			AHL=AHLField.interpolconc(x);%local AHL concentration
			dAHL=AHLField.interpolgrad(x);%local AHL gradient

			if AHL>VthA	%AHL below threshold => low diffusion
				currentSpeedA=speedA(1);
			else		%AHL above threshold => high diffusion
				currentSpeedA=speedA(2);
			end

			newDirection=obj.turn(direction,lambda0A,kappaA,currentSpeedA,AHL,dAHL,timestep);%new direction

			%calculate new position
			if newDirection==0
				xNew=x-currentSpeedA*timestep;
			else
				xNew=x+currentSpeedA*timestep;
			end

			dx=obj.domain(2)-obj.domain(1);
			%correct for going out of boundary
			if xNew < obj.domain(1);
				%wall boundary condition
				%xNew=obj.domain(1);
				%periodic boundary condition
				%xNew=xNew+obj.domain(end);
				xNew=xNew+(obj.domain(end)-obj.domain(1));
			elseif xNew > obj.domain(end)+dx
				%wall boundary condition
				%xNew=obj.domain(end);
				%periodic boundary condition
				%xNew=xNew-obj.domain(end);
				xNew=xNew-(obj.domain(end)-obj.domain(1)+dx);
			end

			%set new position
			bacterium.setxcoordinate(xNew);
			bacterium.setdirection(newDirection);
		end
		end
	end
end
