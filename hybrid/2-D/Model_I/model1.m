classdef model1 < handle
	properties
		%bacterial population A and AHL field
		bacteriaPopA;
		bacteriaPopB;
		AHLField;
		leucineField;

		%constants
		alpha;		%production rate of AHL
		beta;		%production rate of leucine
		k1;			%degradation rate of AHL
		k2;			%degradation rate of leucine
		muA;		%diffusion constant of bacteria A
		muB;		%diffusion constant of bacteria B
		DAHL;		%diffusion constant of AHL
		Dleucine;	%diffusion constant of leucine
		kappaA;		%chemotactic sensitivity of bacteria A
		kappaB;		%chemotactic sensitivity of bacteria B
		VthA;		%threshold concentration of AHL of bacteria A
		VthB;		%threshold concentration of AHL of bacteria B

		%kernel function and bandwidth
		kernelfun;
		bandwidth;
		
		%domain and domain grid
		domain;
		domainGrid;

		%timestep
		dt;

		%array of density arrays & coordinate arrays
		AHLArray;
		leucineArray;
		rhoAArray;
		rhoBArray;
		coordinateAMatrix;
		coordinateBMatrix;
	end

	methods
		function obj=model1(paramA,paramB,paramAHL,paramleucine,...
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
		obj.bacteriaPopA=bacteriaPopulationA(paramA,bacteriaA,domain,domainGrid);
		obj.bacteriaPopB=bacteriaPopulationB(paramB,bacteriaB,domain,domainGrid);
		obj.AHLField=AHL(paramAHL,AHLconcentration,AHLboundaries,domain,domainGrid);
		obj.leucineField=leucine(paramleucine,leucineconcentration,leucineboundaries,domain,domainGrid);

		%record initial density functions, coordinates & AHL field
		%bacteria A
		rhoA=obj.bacteriaPopA.bacteriadensity();
		obj.rhoAArray=rhoA;									%density
		obj.coordinateAMatrix=obj.bacteriaPopA.coordinates();	%coordinates

		%bacteria B
		rhoB=obj.bacteriaPopB.bacteriadensity();
		obj.rhoBArray=rhoB;									%density
		obj.coordinateBMatrix=obj.bacteriaPopB.coordinates();	%coordinates

		%AHL field
		obj.AHLArray=obj.AHLField.getconcentration();			%AHL Field

		%leucine field
		obj.leucineArray=obj.leucineField.getconcentration();	%leucine Field
		end

		function update(obj,dt)
		%Current rho of bacteria A
		rhoAOld=obj.rhoAArray(:,:,end);

		%bacteria A
		%update bacteria positions
		obj.bacteriaPopA.update(obj.AHLField,dt);
		%calculate bacteria density
		rhoA=obj.bacteriaPopA.bacteriadensity();

		%bacteria B
		%update bacteria positions
		obj.bacteriaPopB.update(obj.AHLField,obj.leucineField,dt);
		%calculate bacteria density
		rhoB=obj.bacteriaPopB.bacteriadensity();

		%update AHL field
		obj.AHLField.update(rhoAOld,rhoA,dt);

		%update leucine field
		obj.leucineField.update(rhoAOld,rhoA,dt);

		%record rho
		obj.rhoAArray(:,:,end+1)=rhoA;
		obj.rhoBArray(:,:,end+1)=rhoB;
		%record coordinates
		obj.coordinateAMatrix(:,:,end+1)=obj.bacteriaPopA.coordinates();
		obj.coordinateBMatrix(:,:,end+1)=obj.bacteriaPopB.coordinates();
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
