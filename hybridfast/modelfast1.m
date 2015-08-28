classdef modelfast1 < handle
	properties
		%bacterial population A, B and AHL, leucine fields
		bacteriaPopAB;
		AHLField;
		leucineField;

		%domain and domain grid
		domain;
		domainGrid;

		%array of density arrays & coordinate arrays
		AHLArray;
		leucineArray;
		rhoAArray;
		rhoBArray;
		coordinateAMatrix;
		coordinateBMatrix;
	end

	methods
		function obj=modelfast1(paramAB,paramAHL,paramleucine,...
		bacteriaA,bacteriaB,...
		AHLconcentration,AHLboundaries,leucineconcentration,leucineboundaries,...
		domain)
		%domain
		obj.domain=domain;
		[X,Y]=meshgrid(domain.x,domain.y);
		domainGrid.X=X;
		domainGrid.Y=Y;
		obj.domainGrid=domainGrid;

		%create objects
		obj.bacteriaPopAB=bacteriaPopulationAB(paramAB,bacteriaA,bacteriaB,domain,domainGrid);
		obj.AHLField=AHL(paramAHL,AHLconcentration,AHLboundaries,domain,domainGrid);
		obj.leucineField=leucine(paramleucine,leucineconcentration,leucineboundaries,domain,domainGrid);

		%record initial density functions, coordinates & AHL field
		%bacteria A
		obj.rhoAArray=obj.bacteriaPopAB.bacteriadensityA();		%density
		obj.coordinateAMatrix=obj.bacteriaPopAB.coordinatesA();	%coordinates

		%bacteria B
		obj.rhoBArray=obj.bacteriaPopAB.bacteriadensityB();		%density
		obj.coordinateBMatrix=obj.bacteriaPopAB.coordinatesB();	%coordinates

		%AHL field
		obj.AHLArray=obj.AHLField.getconcentration();			%AHL Field

		%leucine field
		obj.leucineArray=obj.leucineField.getconcentration();	%leucine Field
		end

		function update(obj,dt)
		%Current rho of bacteria A
		rhoAOld=obj.rhoAArray(:,:,end);

		%update bacteria positions
		obj.bacteriaPopAB.updatefast(obj.AHLField,obj.leucineField,dt);
		%calculate bacteria density
		rhoA=obj.bacteriaPopAB.bacteriadensityA();
		rhoB=obj.bacteriaPopAB.bacteriadensityB();

		%update AHL field
		obj.AHLField.update(rhoAOld,rhoA,dt);

		%update leucine field
		obj.leucineField.update(rhoAOld,rhoA,dt);

		%record rho
		obj.rhoAArray(:,:,end+1)=rhoA;
		obj.rhoBArray(:,:,end+1)=rhoB;
		%record coordinates
		obj.coordinateAMatrix(:,:,end+1)=obj.bacteriaPopAB.coordinatesA();
		obj.coordinateBMatrix(:,:,end+1)=obj.bacteriaPopAB.coordinatesB();
		%record AHL field
		obj.AHLArray(:,:,end+1)=obj.AHLField.getconcentration();
		%record leucine field
		obj.leucineArray(:,:,end+1)=obj.leucineField.getconcentration();
		end

		function k=getlength(obj)
		[n,m,k]=size(obj.rhoAArray);
		end

		function domain=getdomain(obj)
		domain=obj.domain;
		end

		function domainGrid=getdomainGrid(obj)
		domainGrid=obj.domainGrid;
		end

		function AHLArray=getAHLArray(obj)
		AHLArray=obj.AHLArray;
		end

		function leucineArray=getleucineArray(obj)
		leucineArray=obj.leucineArray;
		end

		function rhoAArray=getrhoAArray(obj)
		rhoAArray=obj.rhoAArray;
		end

		function rhoBArray=getrhoBArray(obj)
		rhoBArray=obj.rhoBArray;
		end

		function coordinateAMatrix=getcoordinateAMatrix(obj)
		coordinateAMatrix=obj.coordinateAMatrix;
		end

		function coordinateBMatrix=getcoordinateBMatrix(obj)
		coordinateBMatrix=obj.coordinateBMatrix;
		end
	end
end
