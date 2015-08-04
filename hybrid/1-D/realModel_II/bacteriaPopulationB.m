%Model II
classdef bacteriaPopulationB < handle
	properties
		domain;
		bacteria;
	end

	methods
		function obj=bacteriaPopulationB(bacteria,domain)
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

		function newDirection=turn(obj,direction,lambda0B,kappaB,currentSpeedB,leucine,dLeucine,timestep)
		foo=kappaB*currentSpeedB*dLeucine/(2*leucine);
		%disp('dS');
		%dS

		%if direction==0 && dS >0
		%	disp('against the gradient!');
		%end

		%Repelled by leucine
		if direction==0	%moving to the left
			lambda=lambda0B-foo;
		else			%moving to the right
			lambda=lambda0B+foo;
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

		function update(obj,AHLField,leucineField,lambda0B,speedB,kappaB,VthB,timestep)
		%Update bacteria position based on AHL field, leucine field, base turning frequency, constant speed,
		%threshold AHL concentration and chemotactic sensitivity constant
		for bacterium=obj.bacteria
			x=bacterium.getxcoordinate();%x coordinate
			direction=bacterium.getdirection();%direction
			AHL=AHLField.interpolconc(x);%local AHL concentration
			%dAHL=AHLField.interpolgrad(x);%local AHL gradient
			leucine=leucineField.interpolconc(x);%local leucine concentration
			dLeucine=leucineField.interpolgrad(x);%local leucine gradient

			if AHL>VthB	%AHL is above threshold => high diffusion
				currentSpeedB=speedB(2);
			else		%AHL is below threshold => low diffusion
				currentSpeedB=speedB(1);
			end

			newDirection=obj.turn(direction,lambda0B,kappaB,currentSpeedB,leucine,dLeucine,timestep);%new direction

			%calculate new position
			if newDirection==0
				xNew=x-currentSpeedB*timestep;
			else
				xNew=x+currentSpeedB*timestep;
			end

			dx=obj.domain(2)-obj.domain(1);
			%correct for going out of boundary
			if xNew < obj.domain(1);
				%wall boundary condition
				%xNew=obj.domain(1);
				%periodic boundary condition
				xNew=xNew+(obj.domain(end)-obj.domain(1));
			elseif xNew > obj.domain(end)+dx
				%wall boundary condition
				%xNew=obj.domain(end);
				%periodic boundary condition
				xNew=xNew-(obj.domain(end)-obj.domain(1)+dx);
			end

			%set new position
			bacterium.setxcoordinate(xNew);
			bacterium.setdirection(newDirection);
		end
		end
	end
end
