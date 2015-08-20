classdef cellmodel < handle
	properties
		%bacterial population A and AHL field
		bacteriaPopA;

		%constants
		muA;		%diffusion constant of bacteria A

		%kernel function and bandwidth
		kernelfun;
		bandwidth;
		
		%domain and domain grid
		domain;
		domainGrid;

		%timestep
		dt;

		%array of density arrays & coordinate arrays
		rhoAArray;
		coordinateAMatrix;
	end

	methods
		function obj=cellmodel(paramA,bacteriaA,domain)
		%domain
		obj.domain=domain;
		[X,Y]=meshgrid(domain.x,domain.y);
		domainGrid.X=X;
		domainGrid.Y=Y;
		obj.domainGrid=domainGrid;

		%create objects
		obj.bacteriaPopA=bacteriaPopulationA(paramA,bacteriaA,domain,domainGrid);

		%record initial density functions, coordinates & AHL field
		%bacteria A
		rhoA=obj.bacteriaPopA.bacteriadensity();
		obj.rhoAArray=rhoA;									%density
		obj.coordinateAMatrix=obj.bacteriaPopA.coordinates();	%coordinates
		end

		function update(obj,dt)
		%Current rho of bacteria A
		%rhoAOld=obj.rhoAArray(:,:,end);

		%bacteria A
		%update bacteria positions
		%obj.bacteriaPopA.updateraw(obj.AHLField,dt);
		%obj.bacteriaPopA.updateraw(dt);
		obj.bacteriaPopA.updateinteraction(dt);
		%calculate bacteria density
		rhoA=obj.bacteriaPopA.bacteriadensity();

		%record rho
		obj.rhoAArray(:,:,end+1)=rhoA;
		%record coordinates
		obj.coordinateAMatrix(:,:,end+1)=obj.bacteriaPopA.coordinates();
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

		function rhoAArray=getrhoAArray(obj)
		rhoAArray=obj.rhoAArray;
		end

		function coordinateAMatrix=getcoordinateAMatrix(obj)
		coordinateAMatrix=obj.coordinateAMatrix;
		end
	end
end
