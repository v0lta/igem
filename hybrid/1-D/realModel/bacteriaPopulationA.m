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

		function update(obj,AHLField,muA,VthA,kappaA)
		%Update bacteria position based on AHL field, diffusion constant,
		%threshold AHL concentration and chemotactic sensitivity constant
		for bacterium=obj.bacteria
			x=bacterium.getxcoordinate();%x coordinate
			AHL=AHLField.interpolconc(x);%local AHL concentration
			dAHL=AHLField.interpolgrad(x);%local AHL gradient

			if AHL>VthA	%AHL below threshold => low diffusion
				currentMuA=muA(1);
			else		%AHL above threshold => high diffusion
				currentMuA=muA(2);
			end

			%calculate new position
			%attracted by AHL
			xNew=x+currentMuA*kappaA/AHL*dAHL+sqrt(2*currentMuA)*normrnd(0,1);
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
