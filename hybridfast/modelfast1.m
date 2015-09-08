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

		%counter
		N;
		counter;

		%memory
		mFile;
	end

	methods
		function obj=modelfast1(filename,paramAB,paramAHL,paramleucine,...
		bacteriaA,bacteriaB,...
		AHLconcentration,AHLboundaries,leucineconcentration,leucineboundaries,...
		domain,N)
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

		%%initialize initial density functions, coordinates & AHL field
		Jx=numel(domain.x);
		Jy=numel(domain.y);
		Nend=uint16(N+1);
		rhoAArray=zeros(Jy,Jx,Nend);
		rhoBArray=zeros(Jy,Jx,Nend);
		AHLArray=zeros(Jy,Jx,Nend);
		leucineArray=zeros(Jy,Jx,Nend);
		[nBacteriaA,~]=size(bacteriaA);
		[nBacteriaB,~]=size(bacteriaB);
		coordinateAMatrix=zeros(nBacteriaA,2,Nend);
		coordinateBMatrix=zeros(nBacteriaB,2,Nend);
		counter=1;

		%record initial density functions, coordinates & AHL field
		%bacteria A
		obj.rhoAArray(:,:,counter)=obj.bacteriaPopAB.bacteriadensityA();		%density
		%obj.rhoAArray=rhoAArray;
		obj.coordinateAMatrix(:,:,counter)=obj.bacteriaPopAB.coordinatesA();	%coordinates
		%obj.coordinateAMatrix=coordinateAMatrix;

		%bacteria B
		obj.rhoBArray(:,:,counter)=obj.bacteriaPopAB.bacteriadensityB();		%density
		%obj.rhoBArray=rhoBArray;
		obj.coordinateBMatrix(:,:,counter)=obj.bacteriaPopAB.coordinatesB();	%coordinates
		%obj.coordinateBMatrix=coordinateBMatrix;

		%AHL field
		obj.AHLArray(:,:,counter)=obj.AHLField.getconcentration();			%AHL Field
		%obj.AHLArray=AHLArray;

		%leucine field
		obj.leucineArray(:,:,counter)=obj.leucineField.getconcentration();	%leucine Field
		%obj.leucineArray=leucineArray;

		obj.counter=counter;
		end

		function update(obj,dt)

		%update bacteria positions
		obj.bacteriaPopAB.updatefast(obj.AHLField,obj.leucineField,dt);
		%calculate bacteria density
		rhoA=obj.bacteriaPopAB.bacteriadensityA();
		rhoB=obj.bacteriaPopAB.bacteriadensityB();

		%update AHL field
		obj.AHLField.update(rhoA,dt);

		%update leucine field
		obj.leucineField.update(rhoA,dt);

		%increment counter
		counter=obj.counter;
		counter=counter+1;

		%record rho
		obj.rhoAArray(:,:,counter)=rhoA;
		obj.rhoBArray(:,:,counter)=rhoB;
		%record coordinates
		obj.coordinateAMatrix(:,:,counter)=obj.bacteriaPopAB.coordinatesA();
		obj.coordinateBMatrix(:,:,counter)=obj.bacteriaPopAB.coordinatesB();
		%record AHL field
		obj.AHLArray(:,:,counter)=obj.AHLField.getconcentration();
		%record leucine field
		obj.leucineArray(:,:,counter)=obj.leucineField.getconcentration();

		obj.counter=counter;
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
