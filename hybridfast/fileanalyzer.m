classdef fileanalyzer < handle
	properties
		%array of density arrays & coordinate arrays
		mFile;

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
		function obj=fileanalyzer(filename,paramAnal)
		%parameters
		obj.scaling=paramAnal.scaling;
		obj.framerate=paramAnal.framerate;

		%make matfile
		mFile=matfile([filename '_data.mat']);
		
		%number of frames
		obj.nFrames=round(mFile.tend/mFile.dtPDE)+1;

		%domain
		Jx=mFile.Jx;
		Jy=mFile.Jy;
		XLength=mFile.XLength;
		YLength=mFile.YLength;
		domain.x=linspace(0,XLength,Jx);
		domain.y=linspace(0,YLength,Jy);

		[X,Y]=meshgrid(domain.x,domain.y);
		domainGrid.X=X;
		domainGrid.Y=Y;

		obj.domain=domain;
		obj.domainGrid=domainGrid;

		obj.XLength=domain.x(end);
		obj.YLength=domain.y(end);

		obj.mFile=mFile;
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
		mFile=obj.mFile;

		%maxRhoA=mFile.maxRhoA;
		%minRhoA=mFile.minRhoA;
		minRhoA=min(min(min(mFile.rhoAArray)));
		maxRhoA=max(max(max(mFile.rhoAArray)));
		%currentMaxRhoA=max(max(obj.mFile.rhoAArray(:,:,k)));

		%maxRhoB=mFile.maxRhoB;
		%minRhoB=mFile.minRhoB;
		minRhoB=min(min(min(mFile.rhoBArray)));
		maxRhoB=max(max(max(mFile.rhoBArray)));
		%currentMaxRhoB=max(max(obj.mFile.rhoBArray(:,:,k)));

		%maxAHL=mFile.maxAHL;
		%minAHL=mFile.minAHL;
		minAHL=min(min(min(mFile.AHLArray)));
		maxAHL=max(max(max(mFile.AHLArray)));
		%currentMaxAHL=max(max(obj.mFile.AHLArray(:,:,k)));

		%maxleucine=mFile.maxleucine;
		%minleucine=mFile.minleucine;
		minleucine=min(min(min(mFile.leucineArray)));
		maxleucine=max(max(max(mFile.leucineArray)));
		%currentMaxleucine=max(max(obj.mFile.leucineArray(:,:,k)));

		rhoALimit=[minRhoA,maxRhoA];
		%rhoALimit=1;
		rhoBLimit=[minRhoB,maxRhoB];
		%rhoBLimit=1;
		AHLLimit=[minAHL,maxAHL];
		leucineLimit=[minleucine,maxleucine];
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

		%frameArray=cell(nFrames,1);

		%rhoAArray=obj.rhoAArray;
		%coordinateAMatrix=obj.coordinateAMatrix;
		%rhoBArray=obj.rhoBArray;
		%coordinateBMatrix=obj.coordinateBMatrix;

		%AHLArray=obj.AHLArray;
		%leucineArray=obj.leucineArray;

		%plot coarser grid if given grid is too fine
		mFile=obj.mFile;
		Jx=mFile.Jx;
		Jy=mFile.Jy;

		dx=max([round(Jx/100),1]);
		dy=max([round(Jy/100),1]);
		dxi=1:dx:Jx;
		dyi=1:dy:Jy;
		dx=X(1,2)-X(1,1);
		dy=Y(2,1)-Y(1,1);
		X=X(dyi,dxi);
		Y=Y(dyi,dxi);

		%plot only 1000 bacteria
		nBacteriaA=mFile.nBacteriaA;
		nBacteriaB=mFile.nBacteriaB;
		dA=max([round(nBacteriaA/1000),1]);
		dB=max([round(nBacteriaB/1000),1]);
		dAi=1:dA:nBacteriaA;
		dBi=1:dB:nBacteriaB;

		fig=figure('units','normalized','outerposition',[0 0 1 1],'Visible','off');
		%parfor i=1:nFrames
		for i=1:nFrames
			%fig=figure('units','normalized','outerposition',[0 0 1 1],'Visible','on');

			%% -- Bacteria A -- %%
			subplot(2,2,1);
			hold on;

			%density
			rhoA=mFile.rhoAArray(dyi,dxi,i);
			%size(rhoA)
			%mean(mean(rhoA))
			%DataHash(rhoA)
			surf(X,Y,rhoA);
			mesh(X,Y,rhoA,'facecolor','none');
			view(3);

			%bacteria
			coordinateAArray=mFile.coordinateAMatrix(dAi,:,i);

			xCoordinateAArray=coordinateAArray(:,1);
			yCoordinateAArray=coordinateAArray(:,2);
			n=length(xCoordinateAArray);

			plot3(xCoordinateAArray,yCoordinateAArray,zeros(1,n),'k.','MarkerSize',20);

			%formatting
			title('Bacteria A');
			%periodic
			xlim([0,XLength+dx]);
			ylim([0,YLength+dy]);
			zlim(rhoALimit)
			foo = get(gca,'dataaspectratio');
			set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
			xlabel('x');
			ylabel('y');
			zlabel('Density');

			%% -- Bacteria B -- %%
			subplot(2,2,2);
			hold on;

			%density
			rhoB=mFile.rhoBArray(dyi,dxi,i);
			mesh(X,Y,rhoB,'facecolor','none');
			view(3);

			%bacteria
			coordinateBArray=mFile.coordinateBMatrix(dBi,:,i);

			xCoordinateBArray=coordinateBArray(:,1);
			yCoordinateBArray=coordinateBArray(:,2);
			n=length(xCoordinateBArray);

			plot3(xCoordinateBArray,yCoordinateBArray,zeros(1,n),'k.','MarkerSize',20);

			%formatting
			title('Bacteria B');
			%periodic
			xlim([0,XLength+dx]);
			ylim([0,YLength+dy]);
			zlim(rhoBLimit);
			foo = get(gca,'dataaspectratio');
			set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
			xlabel('x');
			ylabel('y');
			zlabel('Density');

			%% -- AHL -- %%
			subplot(2,2,3);
			hold on;

			%concentration
			AHL=mFile.AHLArray(dyi,dxi,i);
			mesh(X,Y,AHL*scaling,'facecolor','none');
			view(3);

			%formatting
			title('AHL');
			xlim([0,XLength]);
			ylim([0,YLength]);
			zlim(AHLLimit*scaling);
			foo = get(gca,'dataaspectratio');
			set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
			xlabel('x');
			ylabel('y');
			zlabel('Concentration');

			%% -- leucine -- %%
			subplot(2,2,4);
			hold on;

			%concentration
			leucine=mFile.leucineArray(dyi,dxi,i);
			mesh(X,Y,leucine*scaling,'facecolor','none');
			view(3);

			%formatting
			title('Leucine');
			xlim([0,XLength]);
			ylim([0,YLength]);
			zlim(leucineLimit*scaling);
			foo = get(gca,'dataaspectratio');
			set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
			xlabel('x');
			ylabel('y');
			zlabel('Concentration');

			writeVideo(vidObj3D,getframe(fig));
			%frameArray{i}=getframe(fig);
			%delete(fig);
			clf;
		end

		%for i=1:nFrames
		%end

		vidObj3D.close();
		end

		function make2Dvideo(obj,filename)
		framerate=obj.framerate;
		nFrames=obj.nFrames;

		XLength=obj.XLength;
		YLength=obj.YLength;

		scaling=obj.scaling;
		[rhoALimit,rhoBLimit,AHLLimit,leucineLimit]=obj.limitoptimizer(1);

		vidObj2D=VideoWriter([filename '_2D.avi']);
		set(vidObj2D,'FrameRate',framerate);
		open(vidObj2D);

		X=obj.domainGrid.X;
		Y=obj.domainGrid.Y;

		%frameArray=cell(nFrames,1);

		%plot coarser grid if given grid is too fine
		mFile=obj.mFile;
		Jx=mFile.Jx;
		Jy=mFile.Jy;

		dx=max([round(Jx/100),1]);
		dy=max([round(Jy/100),1]);
		dxi=1:dx:Jx;
		dyi=1:dy:Jy;
		dx=X(1,2)-X(1,1);
		dy=Y(2,1)-Y(1,1);
		X=X(dyi,dxi);
		Y=Y(dyi,dxi);

		fig=figure('units','normalized','outerposition',[0 0 1 1],'Visible','off');
		%fig=figure('units','normalized','outerposition',[0 0 1 1],'Visible','on');
		for i=1:nFrames
		%parfor i=1:nFrames
			%fig=figure('units','normalized','outerposition',[0 0 1 1],'Visible','off');
			hold on;

			%% -- Bacteria A -- %%
			subplot(2,2,1);

			%density
			rhoA=mFile.rhoAArray(dyi,dxi,i);
			%size(rhoA)
			%mean(mean(rhoA))
			%DataHash(rhoA)
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
			xlim([0,XLength+dx]);
			ylim([0,YLength+dy]);
			zlim(rhoALimit);
			axis image;
			xlabel('x');
			ylabel('y');
			zlabel('Density');
			grid off;
			axis off;

			%% -- Bacteria B -- %%
			subplot(2,2,2);

			%density
			rhoB=mFile.rhoBArray(dyi,dxi,i);
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
			xlim([0,XLength+dx]);
			ylim([0,YLength+dy]);
			zlim([rhoBLimit]);
			axis image;
			xlabel('x');
			ylabel('y');
			zlabel('Density');
			grid off;
			axis off;

			%% -- AHL -- %%
			subplot(2,2,3);

			%concentration
			AHL=mFile.AHLArray(dyi,dxi,i);
			surf(X,Y,AHL*scaling);
			shading('flat');
			view(2);

			%formatting
			title('AHL');
			%legend('Concentration');
			zlim([AHLLimit*scaling]);
			axis image;
			xlabel('x');
			ylabel('y');
			zlabel('Concentration');
			grid off;
			axis off;

			%% -- leucine -- %%
			subplot(2,2,4);

			%concentration
			leucine=mFile.leucineArray(dyi,dxi,i);
			surf(X,Y,leucine*scaling);
			shading('flat');
			view(2);

			%formatting
			title('Leucine');
			%legend('Concentration');
			zlim([leucineLimit*scaling]);
			axis image;
			xlabel('x');
			ylabel('y');
			zlabel('Concentration');
			grid off;
			axis off;

			writeVideo(vidObj2D,getframe(fig));
			%frameArray{i}=getframe(fig);
			%delete(fig);
			clf;
		end

%		for i=1:nFrames
%			writeVideo(vidObj2D,frameArray{i});
%		end

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

		%rhoAArray=obj.rhoAArray;
		%coordinateAMatrix=obj.coordinateAMatrix;
		%rhoBArray=obj.rhoBArray;
		%coordinateBMatrix=obj.coordinateBMatrix;

		%AHLArray=obj.AHLArray;
		%leucineArray=obj.leucineArray;

		%plot coarser grid if given grid is too fine
		mFile=obj.mFile;
		Jx=mFile.Jx;
		Jy=mFile.Jy;

		dx=max([round(Jx/100),1]);
		dy=max([round(Jy/100),1]);
		dxi=1:dx:Jx;
		dyi=1:dy:Jy;
		X=X(dyi,dxi);
		Y=Y(dyi,dxi);

		%plot only 1000 bacteria
		nBacteriaA=mFile.nBacteriaA;
		nBacteriaB=mFile.nBacteriaB;
		dA=max([round(nBacteriaA/1000),1]);
		dB=max([round(nBacteriaB/1000),1]);
		dAi=1:dA:nBacteriaA;
		dBi=1:dB:nBacteriaB;

		fig=figure('units','normalized','outerposition',[0 0 1 1],'Visible','on');
		for i=1:nFrames
			%% -- Bacteria A -- %%
			subplot(2,2,1);
			hold on;

			%density
			rhoA=mFile.rhoAArray(dyi,dxi,i);
			mesh(X,Y,rhoA,'facecolor','none');
			view(3);

			%bacteria
			coordinateAArray=mFile.coordinateAMatrix(dAi,:,i);

			xCoordinateAArray=coordinateAArray(:,1);
			yCoordinateAArray=coordinateAArray(:,2);
			n=length(xCoordinateAArray);

			plot3(xCoordinateAArray,yCoordinateAArray,zeros(1,n),'k.','MarkerSize',20);

			%formatting
			title('Bacteria A');
			%legend('Density','Bacterium');
			xlim([0,XLength]);
			ylim([0,YLength]);
			zlim([rhoALimit])
			foo = get(gca,'dataaspectratio');
			set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
			xlabel('x');
			ylabel('y');
			zlabel('Density');

			%% -- Bacteria B -- %%
			subplot(2,2,2);
			hold on;

			%density
			rhoB=mFile.rhoBArray(dyi,dxi,i);
			mesh(X,Y,rhoB,'facecolor','none');
			view(3);

			%bacteria
			coordinateBArray=mFile.coordinateBMatrix(dBi,:,i);

			xCoordinateBArray=coordinateBArray(:,1);
			yCoordinateBArray=coordinateBArray(:,2);
			n=length(xCoordinateBArray);

			plot3(xCoordinateBArray,yCoordinateBArray,zeros(1,n),'k.','MarkerSize',20);

			%formatting
			title('Bacteria B');
			%legend('Density','Bacterium');
			xlim([0,XLength]);
			ylim([0,YLength]);
			zlim([rhoBLimit]);
			foo = get(gca,'dataaspectratio');
			set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
			xlabel('x');
			ylabel('y');
			zlabel('Density');

			%% -- AHL -- %%
			subplot(2,2,3);
			hold on;

			%concentration
			AHL=mFile.AHLArray(dyi,dxi,i);
			mesh(X,Y,AHL*scaling,'facecolor','none');
			view(3);

			%formatting
			title('AHL');
			%legend('Concentration');
			xlim([0,XLength]);
			ylim([0,YLength]);
			zlim([AHLLimit*scaling]);
			foo = get(gca,'dataaspectratio');
			set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)]);
			xlabel('x');
			ylabel('y');
			zlabel('Concentration');

			%% -- leucine -- %%
			subplot(2,2,4);
			hold on;

			%concentration
			leucine=mFile.leucineArray(dyi,dxi,i);
			mesh(X,Y,leucine*scaling,'facecolor','none');
			view(3);

			%formatting
			title('Leucine');
			%legend('Concentration');
			%[0 leucineLimit*scaling]
			xlim([0,XLength]);
			ylim([0,YLength]);
			zlim([leucineLimit*scaling]);
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

		%frameArray=cell(nFrames,1);

		rhoAArray=obj.rhoAArray;
		coordinateAMatrix=obj.coordinateAMatrix;
		rhoBArray=obj.rhoBArray;
		coordinateBMatrix=obj.coordinateBMatrix;

		AHLArray=obj.AHLArray;
		leucineArray=obj.leucineArray;

		%plot coarser grid if given grid is too fine
		mFile=obj.mFile;
		Jx=mFile.Jx;
		Jy=mFile.Jy;

		dx=max([round(Jx/100),1]);
		dy=max([round(Jy/100),1]);
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
			rhoA=mFile.rhoAArray(dyi,dxi,i);
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
			zlim([rhoALimit]);
			axis image;
			xlabel('x');
			ylabel('y');
			zlabel('Density');

			%% -- Bacteria B -- %%
			subplot(2,2,2);
			hold on;

			%density
			rhoB=mFile.rhoBArray(dyi,dxi,i);
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
			zlim([rhoBLimit]);
			axis image;
			xlabel('x');
			ylabel('y');
			zlabel('Density');

			%% -- AHL -- %%
			subplot(2,2,3);
			hold on;

			%concentration
			AHL=mFile.AHLArray(dyi,dxi,i);
			surf(X,Y,AHL*scaling);
			shading('flat');
			view(2);

			%formatting
			title('AHL');
			%legend('Concentration');
			zlim([AHLLimit*scaling]);
			axis image;
			xlabel('x');
			ylabel('y');
			zlabel('Concentration');

			%% -- leucine -- %%
			subplot(2,2,4);
			hold on;

			%concentration
			leucine=mFile.leucineArray(dyi,dxi,i);
			surf(X,Y,leucine*scaling);
			shading('flat');
			view(2);

			%formatting
			title('Leucine');
			%legend('Concentration');
			zlim([leucineLimit*scaling]);
			axis image;
			xlabel('x');
			ylabel('y');
			zlabel('Concentration');

			pause(0.1);
			clf;
		end
		end

		function makeinteractionvideo(obj,filename);
		%plot bacteria and their radii

		framerate=obj.framerate;
		nFrames=obj.nFrames;

		vidObjinteraction=VideoWriter([filename '_interaction.avi']);
		set(vidObjinteraction,'FrameRate',framerate);
		open(vidObjinteraction);
		
		domain=obj.domain;

		XLength=domain.x(end);
		dx=domain.x(2)-domain.x(1);
		YLength=domain.y(end);
		dy=domain.y(2)-domain.y(1);

		mFile=obj.mFile;
		r0=mFile.r0;

		fig=figure('units','normalized','outerposition',[0 0 1 1],'Visible','off');
		for i=1:nFrames
			hold on;

			coordinateAArray=mFile.coordinateAMatrix(:,:,i);
			coordinateBArray=mFile.coordinateBMatrix(:,:,i);

			xCoordinateAArray=coordinateAArray(:,1);
			yCoordinateAArray=coordinateAArray(:,2);

			xCoordinateBArray=coordinateBArray(:,1);
			yCoordinateBArray=coordinateBArray(:,2);

			xCoordinateArray=[xCoordinateAArray;xCoordinateBArray];
			yCoordinateArray=[yCoordinateAArray;yCoordinateBArray];

			[N,~]=size(xCoordinateArray);

			%plot bacteria
			scatter(xCoordinateAArray,yCoordinateAArray);

			%plot radii
			viscircles([xCoordinateArray,yCoordinateArray],ones(N,1)*r0);

			%formatting
			title('Bacteria A');
			axis image;
			xlim([0 XLength+dx]);
			ylim([0 YLength+dy]);
			xlabel('x');
			ylabel('y');

			writeVideo(vidObjinteraction,getframe(fig));
			clf;
		end

		vidObjinteraction.close();
		end
	end
end
