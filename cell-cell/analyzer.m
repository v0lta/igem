classdef analyzer < handle
	properties
		%array of density arrays & coordinate arrays
		rhoAArray;
		coordinateAMatrix;

		%analyzer parameters
		framerate;
		nFrames;
		tPause;

		%domain
		domain;
		domainGrid;
	end
	
	methods
		function obj=analyzer(paramAnal,model)
		%parameters
		obj.framerate=paramAnal.framerate;
		obj.tPause=paramAnal.tPause;

		%extract history of model
		obj.rhoAArray=model.getrhoAArray();
		obj.coordinateAMatrix=model.getcoordinateAMatrix();
		obj.nFrames=model.getlength();

		%domain
		obj.domain=model.getdomain();
		obj.domainGrid=model.getdomainGrid();
		end

		function preview(obj)
		close all;
		nFrames=obj.nFrames;
		fig=figure();
		set(fig,'units','normalized','outerposition',[0 0 1 1]);
		fig=figure(1);
		for i=1:nFrames
			%3D
			obj.plot3D(i,fig);
			%ylim([-5 75]);
			pause(obj.tPause);
			clf;
		end

		for i=1:nFrames
			%2D
			obj.plot2D(i,fig);
			pause(obj.tPause);
			%ylim([-5 75]);
			clf;
		end
		end

		function makevideo(obj,filename)
		framerate=obj.framerate;
		nFrames=obj.nFrames;

		save(filename);
		disp('Saving workspace and videos');
		vidObj3D=VideoWriter([filename '_3D.avi']);
		set(vidObj3D,'FrameRate',framerate);
		open(vidObj3D);

		vidObj2D=VideoWriter([filename '_2D.avi']);
		set(vidObj2D,'FrameRate',framerate);
		open(vidObj2D);

		fig=figure();
		set(fig,'units','normalized','outerposition',[0 0 1 1]);
		%fig=figure(1);
		for i=1:nFrames
			%3D
			obj.plot3D(i,fig);
			%ylim([-5 75]);
			writeVideo(vidObj3D,getframe(fig));
			clf;
		end

		for i=1:nFrames
			%2D
			obj.plot2D(i,fig);
			%ylim([-5 75]);
			writeVideo(vidObj2D,getframe(fig));
			clf;
		end

		close(fig);
		vidObj3D.close();
		vidObj2D.close();
		end
			
		function plotrhoA3D(obj,k,fig)
		figure(fig);
		hold on;

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		rhoA=obj.rhoAArray(:,:,k);
		%surf(X,Y,rhoA);
		mesh(X,Y,rhoA,'facecolor','none');
		%m=meshc(X,Y,rhoA);
		%set(m,'facecolor','none');
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
		%m=meshc(X,Y,rhoB);
		%set(m,'facecolor','none');
		%surf(X,Y,rhoA);
		view(3);
		%view(-37.5,0);
		
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


		function [rhoALimit]=limitoptimizer(obj,k)
		maxRhoA=max(max(max(obj.rhoAArray)));
		currentMaxRhoA=max(max(obj.rhoAArray(:,:,k)));

		%maxRhoB=max(max(max(obj.rhoBArray)));
		%currentMaxRhoB=max(max(obj.rhoBArray(:,:,k)));

		%rhoALimit=obj.calculatelimit(maxRhoA,currentMaxRhoA);
		%rhoBLimit=obj.calculatelimit(maxRhoB,currentMaxRhoB);
		%AHLLimit=obj.calculatelimit(maxAHL,currentMaxAHL);
		%leucineLimit=obj.calculatelimit(maxleucine,currentMaxleucine);

		rhoALimit=maxRhoA;
		%rhoBLimit=maxRhoB;
		end

		function plot3D(obj,k,fig)

		%scaling=obj.scaling;
		[rhoALimit]=obj.limitoptimizer(k);

		%Bacteria A
		%subplot(2,2,1);
		obj.plotrhoA3D(k,fig);
		obj.plotbacteriaA3D(k,fig);
		title('Bacteria A');
		%legend('Density','Bacterium');
		zlim([0 rhoALimit]);
		foo = get(gca,'dataaspectratio');
		set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
		xlabel('x');
		ylabel('y');
		zlabel('Density');

		%Bacteria B
		%subplot(2,2,2);
		%obj.plotrhoB3D(k,fig);
		%obj.plotbacteriaB3D(k,fig);
		%title('Bacteria B');
		%%legend('Density','Bacterium');
		%zlim([0 rhoBLimit]);
		%foo = get(gca,'dataaspectratio');
		%set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
		%xlabel('x');
		%ylabel('y');
		%zlabel('Density');

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
		shading('flat');
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
		shading('flat');
		view(2);
		
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

		%subplot(2,1,1);
		obj.plotrhoA2D(k,fig);
		%obj.plotbacteriaA2D(k,fig);
		title('Bacteria A');
		axis equal;
		xlabel('x');
		ylabel('y');
		zlabel('Density');

		%subplot(2,1,2);
		%obj.plotrhoB2D(k,fig);
		%%obj.plotbacteriaB2D(k,fig);
		%title('Bacteria B');
		%axis equal;
		%xlabel('x');
		%ylabel('y');
		%zlabel('Density');
		end

	end
end
