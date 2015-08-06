classdef model1 < handle
	properties
		%bacterial population A
		bacteriaPopA;

		%constants
		muA;	%diffusion constant of bacteria A

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
		function obj=model1(bacteriaPopA,muA,kernelfun,bandwidth,dt)
		%process arguments
		obj.bacteriaPopA=bacteriaPopA;
		obj.muA=muA;
		obj.kernelfun=kernelfun;
		obj.bandwidth=bandwidth;
		obj.dt=dt;

		%domain
		obj.domain=bacteriaPopA.domain;
		obj.domainGrid=bacteriaPopA.domainGrid;

		%record initial density functions, coordinates
		%bacteria A
		rhoA=obj.bacteriaPopA.bacteriadensity(obj.kernelfun,obj.bandwidth);
		obj.rhoAArray=rhoA;									%density
		obj.coordinateAMatrix=bacteriaPopA.coordinates();	%coordinates
		end

		function update(obj)

		%bacteria A
		%update bacteria positions
		obj.bacteriaPopA.update(obj.muA,obj.dt);
		%calculate bacteria density
		rhoA=obj.bacteriaPopA.bacteriadensity(obj.kernelfun,obj.bandwidth);

		%record rho
		obj.rhoAArray(:,:,end+1)=rhoA;
		%record coordinates
		obj.coordinateAMatrix(:,:,end+1)=obj.bacteriaPopA.coordinates();
		end

		function k=getlength(obj)
		[n,m,k]=size(obj.rhoAArray);
		end

		function plotrhos(obj,k,fig)
		figure(fig);
		hold on;

		X=obj.domainGrid(:,:,1);
		Y=obj.domainGrid(:,:,2);

		rhoA=obj.rhoAArray(:,:,k);
		%surf(X,Y,rhoA);
		mesh(X,Y,rhoA,'facecolor','none');
		view(3);
		
		hold off;
		end

		function plotbacteria(obj,k,fig)
		figure(fig);
		hold on;

		coordinateAArray=obj.coordinateAMatrix(:,:,k);

		xCoordinateArray=coordinateAArray(1,:);
		yCoordinateArray=coordinateAArray(2,:);
		n=length(xCoordinateArray);

		plot3(xCoordinateArray,yCoordinateArray,zeros(1,n),'k.','MarkerSize',20);

		hold off;
		end

		function plot(obj,k,fig)
		obj.plotrhos(k,fig);
		obj.plotbacteria(k,fig);
		maxRho=max(max(max(obj.rhoAArray)));
		currentMaxRho=max(max(obj.rhoAArray(:,:,k)));
		if currentMaxRho>maxRho/2
			zlim([0 maxRho]);
		elseif currentMaxRho >maxRho/4
			zlim([0 maxRho/2]);
		else
			zlim([0 maxRho/4]);
		end
		legend('Density bacteria A');
		end
	end
end
