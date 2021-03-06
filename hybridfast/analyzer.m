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
		XLength;
		YLength;
	end
	
	methods
		function obj=analyzer(paramAnal,AHLArray,leucineArray,rhoAArray,rhoBArray,coordinateAMatrix,coordinateBMatrix,domain,domainGrid)
		%parameters
		obj.scaling=paramAnal.scaling;
		obj.framerate=paramAnal.framerate;

		%extract history of model
		obj.AHLArray=AHLArray;
		obj.leucineArray=leucineArray;
		obj.rhoAArray=rhoAArray;
		obj.rhoBArray=rhoBArray;
		obj.coordinateAMatrix=coordinateAMatrix;
		obj.coordinateBMatrix=coordinateBMatrix;
		[~,~,obj.nFrames]=size(AHLArray);

		%domain
		obj.domain=domain;
		obj.domainGrid=domainGrid;

		obj.XLength=domain.x(end);
		obj.YLength=domain.y(end);
		end

		function makevideo(obj,filename)
		%make videos in parallel

		obj.make3Dvideo(filename);
		obj.make2Dvideo(filename);

		end

		function preview(obj)
		%make videos in parallel

		obj.preview3Dvideo();
		%obj.preview2Dvideo();

		end

		function [rhoALimit,rhoBLimit,AHLLimit,leucineLimit]=limitoptimizer(obj,k)
		maxRhoA=max(max(max(obj.rhoAArray)));
		%maxRhoA=max([max(max(max(obj.rhoAArray))),1e-2]);
		%currentMaxRhoA=max(max(obj.rhoAArray(:,:,k)));

		maxRhoB=max(max(max(obj.rhoBArray)));
		%maxRhoB=max([max(max(max(obj.rhoBArray))),1e-2]);
		%currentMaxRhoB=max(max(obj.rhoBArray(:,:,k)));

		maxAHL=max(max(max(obj.AHLArray)));
		%maxAHL=max([max(max(max(obj.AHLArray))),1e-2]);
		%currentMaxAHL=max(max(obj.AHLArray(:,:,k)));

		maxleucine=max(max(max(obj.leucineArray)));
		%maxleucine=max([max(max(max(obj.leucineArray))),1e-2]);
		%currentMaxleucine=max(max(obj.leucineArray(:,:,k)));

		rhoALimit=maxRhoA;
		%rhoALimit=1;
		rhoBLimit=maxRhoB;
		%rhoBLimit=1;
		AHLLimit=maxAHL;
		leucineLimit=maxleucine;
		end

		function make3Dvideo(obj,filename)
		framerate=obj.framerate;
		nFrames=obj.nFrames;
		
		XLength=obj.XLength;
		YLength=obj.YLength;

		scaling=obj.scaling;
		[rhoALimit,rhoBLimit,AHLLimit,leucineLimit]=obj.limitoptimizer(1);

		vidObj3D=VideoWriter([filename '_3D.avi']);
		set(vidObj3D,'FrameRate',framerate);
		open(vidObj3D);

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		frameArray=cell(nFrames,1);

		rhoAArray=obj.rhoAArray;
		coordinateAMatrix=obj.coordinateAMatrix;
		rhoBArray=obj.rhoBArray;
		coordinateBMatrix=obj.coordinateBMatrix;

		AHLArray=obj.AHLArray;
		leucineArray=obj.leucineArray;

		%plot coarser grid if given grid is too fine
		Jx=numel(obj.domain.x);
		Jy=numel(obj.domain.y);

		dx=round(Jx/100);
		dy=round(Jy/100);
		if dx<1,dx=1;end;
		if dy<1,dy=1;end;
		dxi=1:dx:Jx;
		dyi=1:dy:Jy;
		X=X(dyi,dxi);
		Y=Y(dyi,dxi);

		%plot only 1000 bacteria
		[nBacteriaA,~,~]=size(coordinateAMatrix);
		[nBacteriaB,~,~]=size(coordinateBMatrix);
		dA=round(nBacteriaA/1000);
		dB=round(nBacteriaB/1000);
		if dA<1,dA=1;end;
		if dB<1,dB=1;end;
		dAi=1:dA:nBacteriaA;
		dBi=1:dB:nBacteriaB;

		parfor i=1:nFrames
			fig=figure('units','normalized','outerposition',[0 0 1 1],'Visible','off');
			%fig=figure('units','normalized','outerposition',[0 0 1 1],'Visible','on');

			%% -- Bacteria A -- %%
			subplot(2,2,1);
			hold on;

			%density
			rhoA=rhoAArray(:,:,i);
			rhoA=rhoA(dyi,dxi);
			mesh(X,Y,rhoA,'facecolor','none');
			view(3);

			%bacteria
			coordinateAArray=coordinateAMatrix(:,:,i);
			coordinateAArray=coordinateAArray(dAi,:);

			xCoordinateAArray=coordinateAArray(:,1);
			yCoordinateAArray=coordinateAArray(:,2);
			n=length(xCoordinateAArray);

			plot3(xCoordinateAArray,yCoordinateAArray,zeros(1,n),'k.','MarkerSize',20);

			%formatting
			title('Bacteria A');
			%xlim([0,XLength]);
			%ylim([0,YLength]);
			zlim([0 rhoALimit])
			foo = get(gca,'dataaspectratio');
			set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
			xlabel('x');
			ylabel('y');
			zlabel('Density');

			%% -- Bacteria B -- %%
			subplot(2,2,2);
			hold on;

			%density
			rhoB=rhoBArray(:,:,i);
			rhoB=rhoB(dyi,dxi);
			mesh(X,Y,rhoB,'facecolor','none');
			view(3);

			%bacteria
			coordinateBArray=coordinateBMatrix(:,:,i);
			coordinateBArray=coordinateBArray(dBi,:);

			xCoordinateBArray=coordinateBArray(:,1);
			yCoordinateBArray=coordinateBArray(:,2);
			n=length(xCoordinateBArray);

			plot3(xCoordinateBArray,yCoordinateBArray,zeros(1,n),'k.','MarkerSize',20);

			%formatting
			title('Bacteria B');
			%xlim([0,XLength]);
			%ylim([0,YLength]);
			zlim([0 rhoBLimit]);
			foo = get(gca,'dataaspectratio');
			set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
			xlabel('x');
			ylabel('y');
			zlabel('Density');

			%% -- AHL -- %%
			subplot(2,2,3);
			hold on;

			%concentration
			AHL=AHLArray(:,:,i);
			AHL=AHL(dyi,dxi);
			mesh(X,Y,AHL*scaling,'facecolor','none');
			view(3);

			%formatting
			title('AHL');
			xlim([0,XLength]);
			ylim([0,YLength]);
			zlim([0 AHLLimit*scaling]);
			foo = get(gca,'dataaspectratio');
			set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
			xlabel('x');
			ylabel('y');
			zlabel('Concentration');

			%% -- leucine -- %%
			subplot(2,2,4);
			hold on;

			%concentration
			leucine=leucineArray(:,:,i);
			leucine=leucine(dyi,dxi);
			mesh(X,Y,leucine*scaling,'facecolor','none');
			view(3);

			%formatting
			title('Leucine');
			xlim([0,XLength]);
			ylim([0,YLength]);
			zlim([0 leucineLimit*scaling]);
			foo = get(gca,'dataaspectratio');
			set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
			xlabel('x');
			ylabel('y');
			zlabel('Concentration');

			frameArray{i}=getframe(fig);
			delete(fig);
		end

		for i=1:nFrames
			writeVideo(vidObj3D,frameArray{i});
		end

		vidObj3D.close();
		end

		function make2Dvideo(obj,filename)
		framerate=obj.framerate;
		nFrames=obj.nFrames;

		scaling=obj.scaling;
		[rhoALimit,rhoBLimit,AHLLimit,leucineLimit]=obj.limitoptimizer(1);

		vidObj2D=VideoWriter([filename '_2D.avi']);
		set(vidObj2D,'FrameRate',framerate);
		open(vidObj2D);

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		frameArray=cell(nFrames,1);

		rhoAArray=obj.rhoAArray;
		coordinateAMatrix=obj.coordinateAMatrix;
		rhoBArray=obj.rhoBArray;
		coordinateBMatrix=obj.coordinateBMatrix;

		AHLArray=obj.AHLArray;
		leucineArray=obj.leucineArray;

		%plot coarser grid if given grid is too fine
		Jx=numel(obj.domain.x);
		Jy=numel(obj.domain.y);

		dx=round(Jx/100);
		dy=round(Jy/100);
		if dx<1,dx=1;end;
		if dy<1,dy=1;end;
		dxi=1:dx:Jx;
		dyi=1:dy:Jy;
		X=X(dyi,dxi);
		Y=Y(dyi,dxi);

		parfor i=1:nFrames
			fig=figure('units','normalized','outerposition',[0 0 1 1],'Visible','off');
			hold on;

			%% -- Bacteria A -- %%
			subplot(2,2,1);

			%density
			rhoA=rhoAArray(:,:,i);
			rhoA=rhoA(dyi,dxi);
			surf(X,Y,rhoA);
			shading('flat');
			view(2);

			%bacteria
			%coordinateAArray=coordinateAMatrix(:,:,i);

			%xCoordinateAArray=coordinateAArray(1,:);
			%yCoordinateAArray=coordinateAArray(2,:);
			%n=length(xCoordinateAArray);

			%plot3(xCoordinateAArray,yCoordinateAArray,zeros(1,n),'k.','MarkerSize',20);

			%formatting
			title('Bacteria A');
			zlim([0 rhoALimit]);
			axis equal;
			xlabel('x');
			ylabel('y');
			zlabel('Density');
			grid off;axis off;

			%% -- Bacteria B -- %%
			subplot(2,2,2);

			%density
			rhoB=rhoBArray(:,:,i);
			rhoB=rhoB(dyi,dxi);
			surf(X,Y,rhoB);
			shading('flat');
			view(2);

			%bacteria
			%coordinateBArray=obj.coordinateBMatrix(:,:,i);

			%xCoordinateBArray=coordinateBArray(1,:);
			%yCoordinateBArray=coordinateBArray(2,:);
			%n=length(xCoordinateBArray);

			%plot3(xCoordinateBArray,yCoordinateBArray,zeros(1,n),'k.','MarkerSize',20);

			%formatting
			title('Bacteria B');
			%legend('Density','Bacterium');
			zlim([0 rhoBLimit]);
			axis equal;
			xlabel('x');
			ylabel('y');
			zlabel('Density');
			grid off;axis off;

			%% -- AHL -- %%
			subplot(2,2,3);

			%concentration
			AHL=AHLArray(:,:,i);
			AHL=AHL(dyi,dxi);
			surf(X,Y,AHL*scaling);
			shading('flat');
			view(2);

			%formatting
			title('AHL');
			%legend('Concentration');
			zlim([0 AHLLimit*scaling]);
			axis equal;
			xlabel('x');
			ylabel('y');
			zlabel('Concentration');
			grid off;axis off;

			%% -- leucine -- %%
			subplot(2,2,4);

			%concentration
			leucine=leucineArray(:,:,i);
			leucine=leucine(dyi,dxi);
			surf(X,Y,leucine*scaling);
			shading('flat');
			view(2);

			%formatting
			title('Leucine');
			%legend('Concentration');
			zlim([0 leucineLimit*scaling]);
			axis equal;
			xlabel('x');
			ylabel('y');
			zlabel('Concentration');
			grid off;axis off;

			frameArray{i}=getframe(fig);
			delete(fig);
		end

		for i=1:nFrames
			writeVideo(vidObj2D,frameArray{i});
		end

		vidObj2D.close();
		end

		function preview3Dvideo(obj)
		framerate=obj.framerate;
		dt=1/framerate;
		nFrames=obj.nFrames;
		
		XLength=obj.XLength;
		YLength=obj.YLength;

		scaling=obj.scaling;
		[rhoALimit,rhoBLimit,AHLLimit,leucineLimit]=obj.limitoptimizer(1);

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		rhoAArray=obj.rhoAArray;
		coordinateAMatrix=obj.coordinateAMatrix;
		rhoBArray=obj.rhoBArray;
		coordinateBMatrix=obj.coordinateBMatrix;

		AHLArray=obj.AHLArray;
		leucineArray=obj.leucineArray;

		%plot coarser grid if given grid is too fine
		Jx=numel(obj.domain.x);
		Jy=numel(obj.domain.y);

		dx=round(Jx/100);
		dy=round(Jy/100);
		if dx<1,dx=1;end;
		if dy<1,dy=1;end;
		dxi=1:dx:Jx;
		dyi=1:dy:Jy;
		X=X(dyi,dxi);
		Y=Y(dyi,dxi);

		%plot only 1000 bacteria
		[nBacteriaA,~,~]=size(coordinateAMatrix);
		[nBacteriaB,~,~]=size(coordinateBMatrix);
		dA=round(nBacteriaA/1000);
		dB=round(nBacteriaB/1000);
		if dA<1,dA=1;end;
		if dB<1,dB=1;end;
		dAi=1:dA:nBacteriaA;
		dBi=1:dB:nBacteriaB;

		fig=figure('units','normalized','outerposition',[0 0 1 1],'Visible','on');
		for i=1:nFrames
			%% -- Bacteria A -- %%
			subplot(2,2,1);
			hold on;

			%density
			rhoA=rhoAArray(:,:,i);
			rhoA=rhoA(dyi,dxi);
			mesh(X,Y,rhoA,'facecolor','none');
			view(3);

			%bacteria
			coordinateAArray=coordinateAMatrix(:,:,i);
			coordinateAArray=coordinateAArray(dAi,:);

			xCoordinateAArray=coordinateAArray(:,1);
			yCoordinateAArray=coordinateAArray(:,2);
			n=length(xCoordinateAArray);

			plot3(xCoordinateAArray,yCoordinateAArray,zeros(1,n),'k.','MarkerSize',20);

			%formatting
			title('Bacteria A');
			%legend('Density','Bacterium');
			xlim([0,XLength]);
			ylim([0,YLength]);
			zlim([0 rhoALimit])
			foo = get(gca,'dataaspectratio');
			set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
			xlabel('x');
			ylabel('y');
			zlabel('Density');

			%% -- Bacteria B -- %%
			subplot(2,2,2);
			hold on;

			%density
			rhoB=rhoBArray(:,:,i);
			rhoB=rhoB(dyi,dxi);
			mesh(X,Y,rhoB,'facecolor','none');
			view(3);

			%bacteria
			coordinateBArray=coordinateBMatrix(:,:,i);
			coordinateBArray=coordinateBArray(dBi,:);

			xCoordinateBArray=coordinateBArray(:,1);
			yCoordinateBArray=coordinateBArray(:,2);
			n=length(xCoordinateBArray);

			plot3(xCoordinateBArray,yCoordinateBArray,zeros(1,n),'k.','MarkerSize',20);

			%formatting
			title('Bacteria B');
			%legend('Density','Bacterium');
			xlim([0,XLength]);
			ylim([0,YLength]);
			zlim([0 rhoBLimit]);
			foo = get(gca,'dataaspectratio');
			set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
			xlabel('x');
			ylabel('y');
			zlabel('Density');

			%% -- AHL -- %%
			subplot(2,2,3);
			hold on;

			%concentration
			AHL=AHLArray(:,:,i);
			AHL=AHL(dyi,dxi);
			mesh(X,Y,AHL*scaling,'facecolor','none');
			view(3);

			%formatting
			title('AHL');
			%legend('Concentration');
			xlim([0,XLength]);
			ylim([0,YLength]);
			zlim([0 AHLLimit*scaling]);
			foo = get(gca,'dataaspectratio');
			set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
			xlabel('x');
			ylabel('y');
			zlabel('Concentration');

			%% -- leucine -- %%
			subplot(2,2,4);
			hold on;

			%concentration
			leucine=leucineArray(:,:,i);
			leucine=leucine(dyi,dxi);
			mesh(X,Y,leucine*scaling,'facecolor','none');
			view(3);

			%formatting
			title('Leucine');
			%legend('Concentration');
			%[0 leucineLimit*scaling]
			xlim([0,XLength]);
			ylim([0,YLength]);
			zlim([0 leucineLimit*scaling]);
			foo = get(gca,'dataaspectratio');
			set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
			xlabel('x');
			ylabel('y');
			zlabel('Concentration');

			pause(dt);
			clf;
		end
		end

		function preview2Dvideo(obj)
		framerate=obj.framerate;
		nFrames=obj.nFrames;

		scaling=obj.scaling;
		[rhoALimit,rhoBLimit,AHLLimit,leucineLimit]=obj.limitoptimizer(1);

		vidObj2D=VideoWriter([filename '_2D.avi']);
		set(vidObj2D,'FrameRate',framerate);
		open(vidObj2D);

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		frameArray=cell(nFrames,1);

		rhoAArray=obj.rhoAArray;
		coordinateAMatrix=obj.coordinateAMatrix;
		rhoBArray=obj.rhoBArray;
		coordinateBMatrix=obj.coordinateBMatrix;

		AHLArray=obj.AHLArray;
		leucineArray=obj.leucineArray;

		%plot coarser grid if given grid is too fine
		Jx=numel(obj.domain.x);
		Jy=numel(obj.domain.y);

		dx=round(Jx/100);
		dy=round(Jy/100);
		if dx<1,dx=1;end;
		if dy<1,dy=1;end;
		dxi=1:dx:Jx;
		dyi=1:dy:Jy;
		X=X(dyi,dxi);
		Y=Y(dyi,dxi);

		fig=figure('units','normalized','outerposition',[0 0 1 1],'Visible','on');
		for i=1:nFrames

			%% -- Bacteria A -- %%
			subplot(2,2,1);
			hold on;

			%density
			rhoA=rhoAArray(:,:,i);
			rhoA=rhoA(dyi,dxi);
			surf(X,Y,rhoA);
			shading('flat');
			view(2);

			%bacteria
			%coordinateAArray=coordinateAMatrix(:,:,i);

			%xCoordinateAArray=coordinateAArray(1,:);
			%yCoordinateAArray=coordinateAArray(2,:);
			%n=length(xCoordinateAArray);

			%plot3(xCoordinateAArray,yCoordinateAArray,zeros(1,n),'k.','MarkerSize',20);

			%formatting
			title('Bacteria A');
			zlim([0 rhoALimit]);
			axis equal;
			xlabel('x');
			ylabel('y');
			zlabel('Density');

			%% -- Bacteria B -- %%
			subplot(2,2,2);
			hold on;

			%density
			rhoB=rhoBArray(:,:,i);
			rhoB=rhoB(dyi,dxi);
			surf(X,Y,rhoB);
			shading('flat');
			view(2);

			%bacteria
			%coordinateBArray=obj.coordinateBMatrix(:,:,i);

			%xCoordinateBArray=coordinateBArray(1,:);
			%yCoordinateBArray=coordinateBArray(2,:);
			%n=length(xCoordinateBArray);

			%plot3(xCoordinateBArray,yCoordinateBArray,zeros(1,n),'k.','MarkerSize',20);

			%formatting
			title('Bacteria B');
			%legend('Density','Bacterium');
			zlim([0 rhoBLimit]);
			axis equal;
			xlabel('x');
			ylabel('y');
			zlabel('Density');

			%% -- AHL -- %%
			subplot(2,2,3);
			hold on;

			%concentration
			AHL=AHLArray(:,:,i);
			AHL=AHL(dyi,dxi);
			surf(X,Y,AHL*scaling);
			shading('flat');
			view(2);

			%formatting
			title('AHL');
			%legend('Concentration');
			zlim([0 AHLLimit*scaling]);
			axis equal;
			xlabel('x');
			ylabel('y');
			zlabel('Concentration');

			%% -- leucine -- %%
			subplot(2,2,4);
			hold on;

			%concentration
			leucine=leucineArray(:,:,i);
			leucine=leucine(dyi,dxi);
			surf(X,Y,leucine*scaling);
			shading('flat');
			view(2);

			%formatting
			title('Leucine');
			%legend('Concentration');
			zlim([0 leucineLimit*scaling]);
			axis equal;
			xlabel('x');
			ylabel('y');
			zlabel('Concentration');

			pause(0.1);
			clf;
		end
		end
	end
end
