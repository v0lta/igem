classdef model1 < handle
	properties
		%bacterial population A and AHL field
		bacteriaPopA;
		bacteriaPopB;
		AHLField;

		%constants
		alpha;	%consumption rate of AHL
		muA;	%diffusion constant of bacteria A
		muB;	%diffusion constant of bacteria B
		DAHL;	%diffusion constant of AHL
		kappaA;	%chemotactic sensitivity of bacteria A
		kappaB;	%chemotactic sensitivity of bacteria B
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
		rhoAArray;
		rhoBArray;
		coordinateAMatrix;
		coordinateBMatrix;

		%scaling for plotting concentrations
		scaling;
	end

	methods
		function obj=model1(bacteriaPopA,bacteriaPopB,AHLField,...
		alpha,muA,muB,DAHL,kappaA,kappaB,VthA,VthB,...
		kernelfun,bandwidth,dt,scaling)
		%process arguments
		obj.bacteriaPopA=bacteriaPopA;
		obj.bacteriaPopB=bacteriaPopB;
		obj.AHLField=AHLField;
		obj.alpha=alpha;
		obj.muA=muA;
		obj.muB=muB;
		obj.DAHL=DAHL;
		obj.kappaA=kappaA;
		obj.kappaB=kappaB;
		obj.VthA=VthA;
		obj.VthB=VthB;
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

		%bacteria B
		rhoB=obj.bacteriaPopB.bacteriadensity(obj.kernelfun,obj.bandwidth);
		obj.rhoBArray=rhoB;									%density
		obj.coordinateBMatrix=bacteriaPopB.coordinates();	%coordinates

		%AHL field
		obj.AHLArray=AHLField.getconcentration();			%AHL Field
		end

		function update(obj)
		%Current rho of bacteria A
		rhoAOld=obj.rhoAArray(:,:,end);

		%bacteria A
		%update bacteria positions
		obj.bacteriaPopA.update(obj.AHLField,obj.muA,obj.VthA,obj.kappaA,obj.dt);
		%calculate bacteria density
		rhoA=obj.bacteriaPopA.bacteriadensity(obj.kernelfun,obj.bandwidth);

		%bacteria B
		%update bacteria positions
		obj.bacteriaPopB.update(obj.AHLField,obj.muB,obj.VthB,obj.kappaB,obj.dt);
		%calculate bacteria density
		rhoB=obj.bacteriaPopB.bacteriadensity(obj.kernelfun,obj.bandwidth);

		%update AHL field
		obj.AHLField.update(rhoAOld,rhoA,obj.DAHL,obj.alpha,obj.dt);

		%record rho
		obj.rhoAArray(:,:,end+1)=rhoA;
		obj.rhoBArray(:,:,end+1)=rhoB;
		%record coordinates
		obj.coordinateAMatrix(:,:,end+1)=obj.bacteriaPopA.coordinates();
		obj.coordinateBMatrix(:,:,end+1)=obj.bacteriaPopB.coordinates();
		%record AHL field
		obj.AHLArray(:,:,end+1)=obj.AHLField.getconcentration();
		end

		function k=getlength(obj)
		[n,m,k]=size(obj.rhoAArray);
		end

		function plotrhoA3D(obj,k,fig)
		figure(fig);
		hold on;

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		rhoA=obj.rhoAArray(:,:,k);
		%surf(X,Y,rhoA);
		mesh(X,Y,rhoA,'facecolor','none');
		%surf(X,Y,rhoA);
		view(3);
		%view(-37.5,0);
		
		hold off;
		end

		function plotrhoB3D(obj,k,fig)
		figure(fig);
		hold on;

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		rhoB=obj.rhoBArray(:,:,k);
		%surf(X,Y,rhoA);
		mesh(X,Y,rhoB,'facecolor','none');
		%surf(X,Y,rhoA);
		view(3);
		%view(-37.5,0);
		
		hold off;
		end

		function plotAHL3D(obj,k,fig,scaling)
		figure(fig);
		hold on;

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		AHLMax=max(max(max(obj.AHLArray)));

		AHL=obj.AHLArray(:,:,k);
		%multiply for scaling
		%mesh(X,Y,(AHL-AHLMax)*scaling,'facecolor','none');
		mesh(X,Y,AHL*scaling,'facecolor','none');
		%surf(X,Y,AHL*scaling);
		view(3);
		%view(-37.5,10);
		%view(-40,0);

		%%double plot
		%plot(domain+domain(end),AHL*obj.scaling);

		hold off;
		end

		function plotbacteriaA3D(obj,k,fig)
		figure(fig);
		hold on;

		coordinateAArray=obj.coordinateAMatrix(:,:,k);

		xCoordinateAArray=coordinateAArray(1,:);
		yCoordinateAArray=coordinateAArray(2,:);
		n=length(xCoordinateAArray);

		plot3(xCoordinateAArray,yCoordinateAArray,zeros(1,n),'k.','MarkerSize',20);

		hold off;
		end

		function plotbacteriaB3D(obj,k,fig)
		figure(fig);
		hold on;

		coordinateBArray=obj.coordinateBMatrix(:,:,k);

		xCoordinateBArray=coordinateBArray(1,:);
		yCoordinateBArray=coordinateBArray(2,:);
		n=length(xCoordinateBArray);

		plot3(xCoordinateBArray,yCoordinateBArray,zeros(1,n),'k.','MarkerSize',20);

		hold off;
		end

		function limit=calculatelimit(obj,maximum,current)
		if current>maximum/2
			limit=maximum;
		elseif current>maximum/4
			limit=maximum/2;
		else
			limit=maximum/4;
		end
		end


		function [rhoALimit,rhoBLimit,AHLLimit]=limitoptimizer(obj,k)
		maxRhoA=max(max(max(obj.rhoAArray)));
		currentMaxRhoA=max(max(obj.rhoAArray(:,:,k)));

		maxRhoB=max(max(max(obj.rhoBArray)));
		currentMaxRhoB=max(max(obj.rhoBArray(:,:,k)));

		maxAHL=max(max(max(obj.AHLArray)));
		currentMaxAHL=max(max(obj.AHLArray(:,:,k)));

		rhoALimit=obj.calculatelimit(maxRhoA,currentMaxRhoA);
		rhoBLimit=obj.calculatelimit(maxRhoB,currentMaxRhoB);
		AHLLimit=obj.calculatelimit(maxAHL,currentMaxAHL);

		end

		function plot3D(obj,k,fig)

		scaling=obj.scaling;
		[rhoALimit,rhoBLimit,AHLLimit]=obj.limitoptimizer(k);

		subplot(1,3,1);
		obj.plotrhoA3D(k,fig);
		obj.plotbacteriaA3D(k,fig);
		title('Bacteria A');
		%legend('Density','Bacterium');
		zlim([0 rhoALimit]);
		%foo = get(gca,'dataaspectratio');
		%set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);

		subplot(1,3,2);
		obj.plotrhoB3D(k,fig);
		obj.plotbacteriaB3D(k,fig);
		title('Bacteria B');
		%legend('Density','Bacterium');
		zlim([0 rhoBLimit]);
		%foo = get(gca,'dataaspectratio');
		%set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);

		subplot(1,3,3);
		obj.plotAHL3D(k,fig,scaling);
		title('AHL');
		%legend('Concentration');
		zlim([0 AHLLimit*scaling]);
		%foo = get(gca,'dataaspectratio');
		%set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
		end

		function plotrhoA2D(obj,k,fig)
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

		function plotrhoB2D(obj,k,fig)
		figure(fig);
		hold on;

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		rhoB=obj.rhoBArray(:,:,k);
		%surf(X,Y,rhoB);
		%mesh(X,Y,rhoB,'facecolor','none');
		surf(X,Y,rhoB);
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

		function plotbacteriaA2D(obj,k,fig)
		figure(fig);
		hold on;

		coordinateAArray=obj.coordinateAMatrix(:,:,k);

		xCoordinateArray=coordinateAArray(1,:);
		yCoordinateArray=coordinateAArray(2,:);
		n=length(xCoordinateArray);

		zl=zlim;
		%marksize=20;
		marksize=1;
		plot3(xCoordinateArray,yCoordinateArray,ones(1,n)*zl(2),'k.','MarkerSize',marksize);

		hold off;
		end

		function plotbacteriaB2D(obj,k,fig)
		figure(fig);
		hold on;

		coordinateBArray=obj.coordinateBMatrix(:,:,k);

		xCoordinateArray=coordinateBArray(1,:);
		yCoordinateArray=coordinateBArray(2,:);
		n=length(xCoordinateArray);

		zl=zlim;
		%marksize=20;
		marksize=1;
		plot3(xCoordinateArray,yCoordinateArray,ones(1,n)*zl(2),'k.','MarkerSize',marksize);

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

		subplot(1,3,1);
		obj.plotrhoA2D(k,fig);
		%obj.plotbacteriaA2D(k,fig);
		title('Bacteria A');
		axis equal;

		subplot(1,3,2);
		obj.plotrhoB2D(k,fig);
		%obj.plotbacteriaB2D(k,fig);
		title('Bacteria B');
		axis equal;

		subplot(1,3,3);
		obj.plotAHL2D(k,fig,scaling);
		title('AHL');
		axis equal;

		end
	end
end
