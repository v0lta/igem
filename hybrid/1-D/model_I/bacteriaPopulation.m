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

		function update(obj,nutrientField,mu,kappa)
		%Update bacteria position based on nutrient field, diffusion constant and chemotactic sensitivity constant
		for bacterium=obj.bacteria
			x=bacterium.getxcoordinate();%x coordinate
			S=nutrientField.interpolconc(x);%local nutrient concentration
			dS=nutrientField.interpolgrad(x);%local nutrient gradient

			%calculate new position
			xNew=x+mu*kappa/S*dS+sqrt(2*mu)*normrnd(0,1);
			%correct for going out of boundary
			if xNew < obj.domain(1);
				xNew=obj.domain(1);
			elseif xNew > obj.domain(end)
				xNew=obj.domain(end);
			end

			%set new position
			bacterium.setxcoordinate(xNew);
		end
		end
	end
end
