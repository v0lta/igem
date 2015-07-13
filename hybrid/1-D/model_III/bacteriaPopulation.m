%Model II
classdef bacteriaPopulation < handle
	properties
		domain;
		bacteria;
	end

	methods
		function obj=bacteriaPopulation(bacteria,domain)
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

		function newDirection=turn(obj,direction,lambda0,kappa,speed,S,dS,timestep)
		foo=kappa*speed*dS/(2*S);
		%disp('dS');
		%dS

		if direction==0	%moving to the left
			lambda=lambda0+foo;
		else			%moving to the right
			lambda=lambda0-foo;
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

		function update(obj,nutrientField,lambda0,kappa,speed,timestep)
		%Update bacteria position based on nutrient field, base turning frequency, chemotactic sensitivity constant & constant speed
		for bacterium=obj.bacteria
			x=bacterium.getxcoordinate();%x coordinate
			direction=bacterium.getdirection();%direction
			S=nutrientField.interpolconc(x);%local nutrient concentration
			dS=nutrientField.interpolgrad(x);%local nutrient gradient
			newDirection=obj.turn(direction,lambda0,kappa,speed,S,dS,timestep);%new direction

			%calculate new position
			if newDirection==0
				xNew=x-speed*timestep;
			else
				xNew=x+speed*timestep;
			end

			%correct for going out of boundary
			%stick boundary condition
			if xNew < obj.domain(1);
				xNew=obj.domain(1);
			elseif xNew > obj.domain(end)
				xNew=obj.domain(end);
			end

			%set new position
			bacterium.setxcoordinate(xNew);
			bacterium.setdirection(newDirection);
		end
		end
	end
end
