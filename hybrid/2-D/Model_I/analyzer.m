classdef analyzer < handle
	properties
		%array of density arrays & coordinate arrays
		AHLArray;
		leucineArray;
		rhoAArray;
		rhoBArray;
		coordinateAMatrix;
		coordinateBMatrix;

		%scaling for plotting concentrations
		scaling;
		framerate;

		%analyzer parameters
		nFrames;

		%domain
		domain;
		domainGrid;
	end
	
	methods
		function obj=analyzer(paramAnal,model)
		%parameters
		obj.scaling=paramAnal.scaling;
		obj.framerate=paramAnal.framerate;

		%extract history of model
		obj.AHLArray=model.getAHLArray();
		obj.leucineArray=model.getleucineArray();
		obj.rhoAArray=model.getrhoAArray();
		obj.rhoBArray=model.getrhoBArray();
		obj.coordinateAMatrix=model.getcoordinateAMatrix();
		obj.coordinateBMatrix=model.getcoordinateBMatrix();
		obj.nFrames=model.getlength();

		%domain
		obj.domain=model.getdomain();
		obj.domainGrid=model.getdomainGrid();
		end

		function interactionpreview(obj)
		close all;
		nFrames=obj.nFrames;
		fig=figure();
		set(fig,'units','normalized','outerposition',[0 0 1 1]);

		for i=1:nFrames
			%2D
			obj.plotinteraction(i,fig);
			pause(0.1);
			%ylim([-5 75]);
			clf;
		end
		end

		function plotinteraction(obj,k,fig)

		figure(fig);
		hold on;

		%cell A
		coordinateAArray=obj.coordinateAMatrix(:,:,k);

		xCoordinateArray=coordinateAArray(1,:);
		yCoordinateArray=coordinateAArray(2,:);
		n=length(xCoordinateArray);

		zl=zlim;
		marksize=20;
		%marksize=1;
		%plot3(xCoordinateArray,yCoordinateArray,ones(1,n)*zl(2),'k.','MarkerSize',marksize);
		plot(xCoordinateArray,yCoordinateArray,'k.','MarkerSize',marksize);

		viscircles(coordinateAArray',1*ones(n,1));
		viscircles(coordinateAArray',1.25*ones(n,1));

		%cell B
		coordinateBArray=obj.coordinateBMatrix(:,:,k);

		xCoordinateArray=coordinateBArray(1,:);
		yCoordinateArray=coordinateBArray(2,:);
		n=length(xCoordinateArray);

		zl=zlim;
		marksize=20;
		%marksize=1;
		%plot3(xCoordinateArray,yCoordinateArray,ones(1,n)*zl(2),'k.','MarkerSize',marksize);
		plot(xCoordinateArray,yCoordinateArray,'k.','MarkerSize',marksize);

		viscircles(coordinateBArray',1*ones(n,1));
		%viscircles(coordinateAArray',1.25*ones(n,1));
		hold off;
		end

		function preview(obj)
		close all;
		nFrames=obj.nFrames;
		fig=figure('units','normalized','outerposition',[0 0 1 1]);
		%set(fig,'units','normalized','outerposition',[0 0 1 1]);
		%fig=figure(1);
		for i=1:nFrames
			%3D
			obj.plot3D(i);
			%ylim([-5 75]);
			pause(0.1);
			clf;
		end

		for i=1:nFrames
			%2D
			obj.plot2D(i);
			pause(0.1);
			%ylim([-5 75]);
			clf;
		end
		end

		function makevideo(obj,filename)
		framerate=obj.framerate;
		nFrames=obj.nFrames;

		%save(filename);
		vidObj3D=VideoWriter([filename '_3D.avi']);
		set(vidObj3D,'FrameRate',framerate);
		open(vidObj3D);

		vidObj2D=VideoWriter([filename '_2D.avi']);
		set(vidObj2D,'FrameRate',framerate);
		open(vidObj2D);

		fig=figure('units','normalized','outerposition',[0 0 1 1],'Visible','off');
		%set(fig,'units','normalized','outerposition',[0 0 1 1],'Visible','off');
		%fig=figure(1);
		for i=1:nFrames
			%3D
			%obj.plot3D(i,fig);
			obj.plot3D(i);
			%ylim([-5 75]);
			writeVideo(vidObj3D,getframe(fig));
			clf;
		end
		vidObj3D.close();

		for i=1:nFrames
			%2D
			%obj.plot2D(i,fig);
			obj.plot2D(i);
			%ylim([-5 75]);
			writeVideo(vidObj2D,getframe(fig));
			clf;
		end
		vidObj2D.close();

		close(fig);
		end
			
		%function plotrhoA3D(obj,k,fig)
		function plotrhoA3D(obj,k)
		%figure(fig);
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

		%function plotrhoB3D(obj,k,fig)
		function plotrhoB3D(obj,k)
		%figure(fig);
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

		%function plotAHL3D(obj,k,fig,scaling)
		function plotAHL3D(obj,k,scaling)
		%figure(fig);
		hold on;

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		AHLMax=max(max(max(obj.AHLArray)));

		AHL=obj.AHLArray(:,:,k);
		%multiply for scaling
		%mesh(X,Y,(AHL-AHLMax)*scaling,'facecolor','none');
		mesh(X,Y,AHL*scaling,'facecolor','none');
		%m=meshc(X,Y,AHL*scaling);
		%set(m,'facecolor','none');
		%surf(X,Y,AHL*scaling);
		view(3);
		%view(-37.5,10);
		%view(-40,0);

		%%double plot
		%plot(domain+domain(end),AHL*obj.scaling);

		hold off;
		end

		%function plotleucine3D(obj,k,fig,scaling)
		function plotleucine3D(obj,k,scaling)
		%figure(fig);
		hold on;

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		leucineMax=max(max(max(obj.leucineArray)));

		leucine=obj.leucineArray(:,:,k);
		%multiply for scaling
		%mesh(X,Y,(leucine-leucineMax)*scaling,'facecolor','none');
		mesh(X,Y,leucine*scaling,'facecolor','none');
		%m=meshc(X,Y,leucine*scaling);
		%set(m,'facecolor','none');
		%surf(X,Y,leucine*scaling);
		view(3);
		%view(-37.5,10);
		%view(-40,0);

		%%double plot
		%plot(domain+domain(end),leucine*obj.scaling);

		hold off;
		end

		%function plotbacteriaA3D(obj,k,fig)
		function plotbacteriaA3D(obj,k)
		%figure(fig);
		hold on;

		coordinateAArray=obj.coordinateAMatrix(:,:,k);

		xCoordinateAArray=coordinateAArray(1,:);
		yCoordinateAArray=coordinateAArray(2,:);
		n=length(xCoordinateAArray);

		plot3(xCoordinateAArray,yCoordinateAArray,zeros(1,n),'k.','MarkerSize',20);

		hold off;
		end

		%function plotbacteriaB3D(obj,k,fig)
		function plotbacteriaB3D(obj,k)
		%figure(fig);
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


		function [rhoALimit,rhoBLimit,AHLLimit,leucineLimit]=limitoptimizer(obj,k)
		maxRhoA=max(max(max(obj.rhoAArray)));
		currentMaxRhoA=max(max(obj.rhoAArray(:,:,k)));

		maxRhoB=max(max(max(obj.rhoBArray)));
		currentMaxRhoB=max(max(obj.rhoBArray(:,:,k)));

		maxAHL=max(max(max(obj.AHLArray)));
		currentMaxAHL=max(max(obj.AHLArray(:,:,k)));

		maxleucine=max(max(max(obj.leucineArray)));
		currentMaxleucine=max(max(obj.leucineArray(:,:,k)));

		%rhoALimit=obj.calculatelimit(maxRhoA,currentMaxRhoA);
		%rhoBLimit=obj.calculatelimit(maxRhoB,currentMaxRhoB);
		%AHLLimit=obj.calculatelimit(maxAHL,currentMaxAHL);
		%leucineLimit=obj.calculatelimit(maxleucine,currentMaxleucine);

		rhoALimit=maxRhoA;
		%rhoALimit=1;
		rhoBLimit=maxRhoB;
		%rhoBLimit=1;
		AHLLimit=maxAHL;
		leucineLimit=maxleucine;
		end

		%function plot3D(obj,k,fig)
		function plot3D(obj,k)

		scaling=obj.scaling;
		[rhoALimit,rhoBLimit,AHLLimit,leucineLimit]=obj.limitoptimizer(k);

		%Bacteria A
		subplot(2,2,1);
		%obj.plotrhoA3D(k,fig);
		obj.plotrhoA3D(k);
		%obj.plotbacteriaA3D(k,fig);
		obj.plotbacteriaA3D(k);
		title('Bacteria A');
		%legend('Density','Bacterium');
		zlim([0 rhoALimit]);
		foo = get(gca,'dataaspectratio');
		set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
		xlabel('x');
		ylabel('y');
		zlabel('Density');

		%Bacteria B
		subplot(2,2,2);
		%obj.plotrhoB3D(k,fig);
		obj.plotrhoB3D(k);
		%obj.plotbacteriaB3D(k,fig);
		obj.plotbacteriaB3D(k);
		title('Bacteria B');
		%legend('Density','Bacterium');
		zlim([0 rhoBLimit]);
		foo = get(gca,'dataaspectratio');
		set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
		xlabel('x');
		ylabel('y');
		zlabel('Density');

		%AHL
		subplot(2,2,3);
		%obj.plotAHL3D(k,fig,scaling);
		obj.plotAHL3D(k,scaling);
		title('AHL');
		%legend('Concentration');
		zlim([0 AHLLimit*scaling]);
		foo = get(gca,'dataaspectratio');
		set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
		xlabel('x');
		ylabel('y');
		zlabel('Concentration');

		%leucine
		subplot(2,2,4);
		%obj.plotleucine3D(k,fig,scaling);
		obj.plotleucine3D(k,scaling);
		title('Leucine');
		%legend('Concentration');
		zlim([0 leucineLimit*scaling]);
		foo = get(gca,'dataaspectratio');
		set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
		xlabel('x');
		ylabel('y');
		zlabel('Concentration');
		end

		function plotrhoA2D(obj,k)
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

		function plotrhoB2D(obj,k)
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

		function plotAHL2D(obj,k,scaling)
		hold on;

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		AHL=obj.AHLArray(:,:,k);
		%multiply for scaling
		%mesh(X,Y,AHL*scaling,'facecolor','none');
		surf(X,Y,AHL*scaling);
		shading('flat');
		view(2);

		%%double plot
		%plot(domain+domain(end),AHL*obj.scaling);

		hold off;
		end

		function plotleucine2D(obj,k,scaling)
		hold on;

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		leucine=obj.leucineArray(:,:,k);
		%multiply for scaling
		%mesh(X,Y,leucine*scaling,'facecolor','none');
		surf(X,Y,leucine*scaling);
		shading('flat');
		view(2);

		%%double plot
		%plot(domain+domain(end),leucine*obj.scaling);

		hold off;
		end

		function plotbacteriaA2D(obj,k)
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

		function plotbacteriaB2D(obj,k)
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

		function plot2D(obj,k)

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

		subplot(2,2,1);
		obj.plotrhoA2D(k);
		%obj.plotbacteriaA2D(k,fig);
		title('Bacteria A');
		axis equal;
		xlabel('x');
		ylabel('y');
		zlabel('Density');

		subplot(2,2,2);
		obj.plotrhoB2D(k);
		%obj.plotbacteriaB2D(k,fig);
		title('Bacteria B');
		axis equal;
		xlabel('x');
		ylabel('y');
		zlabel('Density');

		subplot(2,2,3);
		obj.plotAHL2D(k,scaling);
		title('AHL');
		axis equal;
		xlabel('x');
		ylabel('y');
		zlabel('Concentration');

		subplot(2,2,4);
		obj.plotleucine2D(k,scaling);
		title('Leucine');
		axis equal;
		xlabel('x');
		ylabel('y');
		zlabel('Concentration');
		end
	end
end
