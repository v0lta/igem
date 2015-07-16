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

		function update(obj,AHLField,leucineField,muB,VthB,kappaB)
		%Update bacteria position based on AHL field, leucine field, diffusion constant,
		%threshold AHL concentration and chemotactic sensitivity constant
		for bacterium=obj.bacteria
			x=bacterium.getxcoordinate();%x coordinate
			AHL=AHLField.interpolconc(x);%local AHL concentration
			%dAHL=AHLField.interpolgrad(x);%local AHL gradient

			leucine=leucineField.interpolconc(x);%local leucine concentration
			dLeucine=leucineField.interpolgrad(x);%local leucine gradient

			if AHL>VthB	%AHL is above threshold => high diffusion
				currentMuB=muB(2);
			else		%AHL is below threshold => low diffusion
				currentMuB=muB(1);
			end

			%calculate new position
			%repelled by leucine
			%x
			%currentMuB
			%kappaB
			%leucine
			%dLeucine
			xNew=x-currentMuB*kappaB/leucine*dLeucine+sqrt(2*currentMuB)*normrnd(0,1);
			%correct for going out of boundary
			if xNew < obj.domain(1);
				%wall boundary condition
				%xNew=obj.domain(1);
				%periodic boundary condition
				xNew=xNew+obj.domain(end);
			elseif xNew > obj.domain(end)
				%wall boundary condition
				%xNew=obj.domain(end);
				%periodic boundary condition
				xNew=xNew-obj.domain(end);
			end

			%set new position
			bacterium.setxcoordinate(xNew);
		end
		end
	end
end
