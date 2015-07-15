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

		function update(obj,attractantField,mu,Vth,kappa)
		%Update bacteria position based on attractant field, diffusion constant and chemotactic sensitivity constant
		for bacterium=obj.bacteria
			x=bacterium.getxcoordinate();%x coordinate
			S=attractantField.interpolconc(x);%local attractant concentration
			dS=attractantField.interpolgrad(x);%local attractant gradient

			if S>Vth %low diffusion
				currentMu=mu(1);
			else	%high diffusion
				currentMu=mu(2);
			end

			%calculate new position
			xNew=x+currentMu*kappa/S*dS+sqrt(2*currentMu)*normrnd(0,1);
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
