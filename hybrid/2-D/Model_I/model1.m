classdef model1 < handle
	properties
		%bacterial population A and AHL field
		bacteriaPopA;
		AHLField;

		%constants
		alpha;	%consumption rate of AHL
		muA;	%diffusion constant of bacteria A
		DAHL;	%diffusion constant of AHL
		kappaA;	%chemotactic sensitivity of bacteria A

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
		rhoAArray;
		coordinateAMatrix;

		%scaling for plotting concentrations
		scaling;
	end

	methods
		function obj=model1(bacteriaPopA,AHLField,alpha,muA,DAHL,kappaA,kernelfun,bandwidth,dt,scaling)
		%process arguments
		obj.bacteriaPopA=bacteriaPopA;
		obj.AHLField=AHLField;
		obj.alpha=alpha;
		obj.muA=muA;
		obj.DAHL=DAHL;
		obj.kappaA=kappaA;
		obj.kernelfun=kernelfun;
		obj.bandwidth=bandwidth;
		obj.dt=dt;
		obj.scaling=scaling;

		%domain
		obj.domain=bacteriaPopA.domain;
		obj.domainGrid=bacteriaPopA.domainGrid;

		%record initial density functions, coordinates & AHL field
		%bacteria A
		rhoA=obj.bacteriaPopA.bacteriadensity(obj.kernelfun,obj.bandwidth);
		obj.rhoAArray=rhoA;									%density
		obj.coordinateAMatrix=bacteriaPopA.coordinates();	%coordinates

		%AHL field
		obj.AHLArray=AHLField.getconcentration();			%AHL Field
		end

		function update(obj)
		%Current rho
		rhoAOld=obj.rhoAArray(:,:,end);

		%bacteria A
		%update bacteria positions
		obj.bacteriaPopA.update(obj.AHLField,obj.muA,obj.dt,obj.kappaA);
		%calculate bacteria density
		rhoANew=obj.bacteriaPopA.bacteriadensity(obj.kernelfun,obj.bandwidth);

		%update AHL field
		obj.AHLField.update(rhoAOld,rhoANew,obj.DAHL,obj.alpha,obj.dt);

		%record rho
		obj.rhoAArray(:,:,end+1)=rhoANew;
		%record coordinates
		obj.coordinateAMatrix(:,:,end+1)=obj.bacteriaPopA.coordinates();
		%record AHL field
		obj.AHLArray(:,:,end+1)=obj.AHLField.getconcentration();
		end

		function k=getlength(obj)
		[n,m,k]=size(obj.rhoAArray);
		end

		function plotrhos3D(obj,k,fig)
		figure(fig);
		hold on;

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		rhoA=obj.rhoAArray(:,:,k);
		%surf(X,Y,rhoA);
		mesh(X,Y,rhoA,'facecolor','none');
		%surf(X,Y,rhoA);
		view(3);
		
		hold off;
		end

		function plotAHL3D(obj,k,fig,scaling)
		figure(fig);
		hold on;

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		AHL=obj.AHLArray(:,:,k);
		%multiply for scaling
		mesh(X,Y,AHL*scaling,'facecolor','none');
		%surf(X,Y,AHL*scaling);
		view(3);

		%%double plot
		%plot(domain+domain(end),AHL*obj.scaling);

		hold off;
		end

		function plotbacteria3D(obj,k,fig)
		figure(fig);
		hold on;

		coordinateAArray=obj.coordinateAMatrix(:,:,k);

		xCoordinateArray=coordinateAArray(1,:);
		yCoordinateArray=coordinateAArray(2,:);
		n=length(xCoordinateArray);

		plot3(xCoordinateArray,yCoordinateArray,zeros(1,n),'k.','MarkerSize',20);

		hold off;
		end

		function plot3D(obj,k,fig)

		scaling=obj.scaling;
		maxRho=max(max(max(obj.rhoAArray)));
		currentMaxRho=max(max(obj.rhoAArray(:,:,k)));
		maxAHL=max(max(max(obj.AHLArray)));
		currentMaxAHL=max(max(obj.AHLArray(:,:,k)));
		zlim([0 max([maxRho maxAHL*scaling])]);

		%if maxRho>maxAHL*obj.scaling
		%	if currentMaxRho>maxRho/2
		%		zlim([0 maxRho]);
		%		scaling=obj.scaling;
		%	elseif currentMaxRho >maxRho/4
		%		zlim([0 maxRho/2]);
		%		scaling=obj.scaling/2;
		%	else
		%		zlim([0 maxRho/4]);
		%		scaling=obj.scaling/4;
		%	end
		%else
		%	if currentMaxAHL>maxAHL/2
		%		scaling=obj.scaling;
		%		zlim([0 maxAHL*scaling]);
		%	elseif currentMaxAHL >maxAHL/4
		%		scaling=obj.scaling/2;
		%		zlim([0 maxAHL*scaling]);
		%	else
		%		scaling=obj.scaling/4;
		%		zlim([0 maxAHL*scaling]);
		%	end
		%end


		obj.plotrhos3D(k,fig);
		obj.plotAHL3D(k,fig,scaling);
		obj.plotbacteria3D(k,fig);
		legend('Density bacteria A','AHL');
		end

		function plotrhos2D(obj,k,fig)
		figure(fig);
		hold on;

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		rhoA=obj.rhoAArray(:,:,k);
		%surf(X,Y,rhoA);
		%mesh(X,Y,rhoA,'facecolor','none');
		surf(X,Y,rhoA);
		view(2);
		
		hold off;
		end

		function plotAHL2D(obj,k,fig,scaling)
		figure(fig);
		hold on;

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		AHL=obj.AHLArray(:,:,k);
		%multiply for scaling
		%mesh(X,Y,AHL*scaling,'facecolor','none');
		surf(X,Y,AHL*scaling);
		view(2);

		%%double plot
		%plot(domain+domain(end),AHL*obj.scaling);

		hold off;
		end

		function plotbacteria2D(obj,k,fig)
		figure(fig);
		hold on;

		coordinateAArray=obj.coordinateAMatrix(:,:,k);

		xCoordinateArray=coordinateAArray(1,:);
		yCoordinateArray=coordinateAArray(2,:);
		n=length(xCoordinateArray);

		plot3(xCoordinateArray,yCoordinateArray,ones(1,n),'k.','MarkerSize',20);

		hold off;
		end

		function plot2D(obj,k,fig)

		scaling=1;
		%maxRho=max(max(max(obj.rhoAArray)));
		%currentMaxRho=max(max(obj.rhoAArray(:,:,k)));
		%if currentMaxRho>maxRho/2
		%	zlim([0 maxRho]);
		%	scaling=obj.scaling;
		%elseif currentMaxRho >maxRho/4
		%	zlim([0 maxRho/2]);
		%	scaling=obj.scaling/2;
		%else
		%	zlim([0 maxRho/4]);
		%	scaling=obj.scaling/4;
		%end

		subplot(1,2,1);
		obj.plotrhos2D(k,fig);
		subplot(1,2,2);
		obj.plotAHL2D(k,fig,scaling);
		legend('AHL');
		subplot(1,2,1);
		obj.plotbacteria2D(k,fig);
		legend('Density bacteria A');
		end
	end
end
