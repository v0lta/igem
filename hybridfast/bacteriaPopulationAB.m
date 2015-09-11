classdef bacteriaPopulationAB < handle
	properties
		%density
		kernelfun;
		bandwidth;

		%bacteria A
		muA;
		VthA;
		kappaA;

		%bacteria B
		muB;
		VthB;
		kappaB;

		%common bacteria parameters
		r0;			%cell radius
		rcut;		%cut off radius for attraction
		rsearch;	%search radius
		rsearch2;	%square of search radius
		k1;	%spring constant
		k2;	%spring constant
		k3;	%spring constant
		gamma;	%friction parameter

		%bacteria arrays
		bactX;
		bactY;
		bactNb;
		numA;
		numB;
		%periodic
		virtBactNb;

		%domain
		domain;
		domainGrid;

		%extended domain and index arrays for periodic boundary conditions
		domainExtended;
		domainGridExtended;

		dxi;
		dsx1;
		dsxend;

		dyi;
		dsy1;
		dsyend;

		Sx;
		Sy;

		Jx;
		Jy;

		%neighbor search helper variables
		counter;
		modulo;
	end

	methods
		function obj=bacteriaPopulationAB(paramAB,bacteriaA,bacteriaB,domain,domainGrid)
		%parameters
		obj.kernelfun=paramAB.kernelfun;
		obj.bandwidth=paramAB.bandwidth;
		obj.muA=paramAB.muA;
		obj.VthA=paramAB.VthA;
		obj.kappaA=paramAB.kappaA;
		obj.muB=paramAB.muB;
		obj.VthB=paramAB.VthB;
		obj.kappaB=paramAB.kappaB;
		obj.r0=paramAB.r0;
		obj.rcut=paramAB.rcut;
		obj.rsearch=paramAB.rsearch;
		obj.rsearch2=paramAB.rsearch^2;
		obj.k1=paramAB.k1;
		obj.k2=paramAB.k2;
		obj.k3=paramAB.k3;
		obj.gamma=paramAB.gamma;
		obj.modulo=paramAB.modulo;

		%domain
		obj.domain=domain;
		obj.domainGrid=domainGrid;

		%extended domain for periodic boundary conditions
		dx=domain.x(2)-domain.x(1);
		Jx=numel(domain.x);
		Sx=floor(paramAB.bandwidth/dx);
		Lx=domain.x(end);
		domainExtended.x=linspace(-Sx*dx,Lx+(Sx+1)*dx,Jx+2*Sx+1);

		dy=domain.y(2)-domain.y(1);
		Jy=numel(domain.y);
		Sy=floor(paramAB.bandwidth/dy);
		Ly=domain.y(end);
		domainExtended.y=linspace(-Sy*dy,Ly+(Sy+1)*dy,Jy+2*Sy+1);

		[X,Y]=meshgrid(domainExtended.x,domainExtended.y);
		domainGridExtended.X=X;
		domainGridExtended.Y=Y;

		obj.domainExtended=domainExtended;
		obj.domainGridExtended=domainGridExtended;

		obj.dxi=Sx+1:Sx+Jx;
		obj.dsx1=1:Sx;
		obj.dsxend=Sx+Jx+1:Jx+2*Sx+1;

		obj.dyi=Sy+1:Sy+Jy;
		obj.dsy1=1:Sy;
		obj.dsyend=Sy+Jy+1:Jy+2*Sy+1;

		obj.Sx=Sx;
		obj.Sy=Sy;

		obj.Jx=Jx;
		obj.Jy=Jy;

		%neighbor variables
		obj.counter=0;

		%initialize bacteria arrays & counters
		obj.bactX=[];
		obj.bactY=[];
		obj.numA=0;
		obj.numB=0;
		obj.bactNb={};

		%add given bacteria
		obj.addbacteria(bacteriaA,bacteriaB);
		end

		function addbacteria(obj,bacteriaA,bacteriaB)
		%add list of bacteria

		numA=obj.numA;
		numB=obj.numB;
		N=numA+numB;

		[nA,~]=size(bacteriaA);
		[nB,~]=size(bacteriaB);

		if N==0
			obj.bactX=[bacteriaA(:,1);bacteriaB(:,1)];
			obj.bactY=[bacteriaA(:,2);bacteriaB(:,2)];
		else
			obj.bactX=[obj.bactX(1:numA);bacteriaA(:,1);obj.bactX(numA+1:end);bacteriaB(:,1)];
			obj.bactY=[obj.bactY(1:numA);bacteriaA(:,2);obj.bactY(numA+1:end);bacteriaB(:,2)];
		end

		obj.numA=numA+nA;
		obj.numB=numB+nB;
		end


		function addbacteriumA(obj,x,y)
		%add a bacterium of type A

		numA=obj.numA;
		bactX=obj.bactX;
		bactY=obj.bactY;
		if numA~=0
			precedingAX=bactX(1:numA);
			precedingAY=bactY(1:numA);
		else
			precedingAX=[];
			precedingAY=[];
		end

		if numB~=0
			succeedingBX=bactX(numA+1:end);
			succeedingBY=bactY(numA+1:end);
		else
			succeedingBX=[];
			succeedingBY=[];
		end

		obj.bactX=[precedingAX;x;succeedingBX];
		obj.bactY=[precedingAY;y;succeedingBY];

		obj.numA=numA+1;
		end

		function addbacteriumB(obj,x,y)
		%add a bacterium of type B

		obj.bactX=[obj.bactX;x];
		obj.bactY=[obj.bactY;y];
		obj.numB=numB+1;
		end

		function coordinateArrayA=coordinatesA(obj);
		%return coordinate array of bacteria A
		numA=obj.numA;
		coordinateArrayA=[obj.bactX(1:numA),obj.bactY(1:numA)];
		end
	
		function coordinateArrayB=coordinatesB(obj);
		%return coordinate array of bacteria B
		numA=obj.numA;
		coordinateArrayB=[obj.bactX(numA+1:end),obj.bactY(numA+1:end)];
		end

		function coordinateArray=coordinates(obj);
		%return coordinate array of all bacteria
		coordinateArray=[obj.bactX,obj.bactY];
		end

		function rhoA=bacteriadensityA(obj)
		%return bacteria A density array evaluated on grid points of domain
		%coordinateArray=obj.coordinatesA();
		%rhoA=KDE2D(coordinateArray,obj.kernelfun,obj.domainGrid.X,obj.domainGrid.Y,obj.bandwidth);
		rhoA=obj.bacteriadensityAperiodic();
		end

		function rhoB=bacteriadensityB(obj)
		%return bacteria B density array evaluated on grid points of domain
		%coordinateArray=obj.coordinatesB();
		%rhoB=KDE2D(coordinateArray,obj.kernelfun,obj.domainGrid.X,obj.domainGrid.Y,obj.bandwidth);
		rhoB=obj.bacteriadensityBperiodic();
		end

		function rho=bacteriadensity(obj)
		%return total bacteria density array evaluated on grid points of domain
		%coordinateArray=obj.coordinates();
		%rho=KDE2D(coordinateArray,obj.kernelfun,obj.domainGrid.X,obj.domainGrid.Y,obj.bandwidth);
		rho=obj.bacteriadensityperiodic();
		end

		function rhoA=bacteriadensityAperiodic(obj)
		%return bacteria A density array evaluated on grid points of periodic domain
		coordinateArray=obj.coordinatesA();
		rhoAExtended=KDE2D(coordinateArray,obj.kernelfun,...
			obj.domainGridExtended.X,obj.domainGridExtended.Y,obj.bandwidth);

		%load index arrays
		dxi=obj.dxi;
		dsx1=obj.dsx1;
		dsxend=obj.dsxend;

		dyi=obj.dyi;
		dsy1=obj.dsy1;
		dsyend=obj.dsyend;

		Sx=obj.Sx;
		Sy=obj.Sy;

		Jx=obj.Jx;
		Jy=obj.Jy;

		%center part
		rhoA=rhoAExtended(dyi,dxi);

		%left
		%Sx
		%dsxend
		rhoA(:,1:Sx+1)=rhoA(:,1:Sx+1)+rhoAExtended(dyi,dsxend);
		%right
		rhoA(:,Jx-Sx+1:Jx)=rhoA(:,Jx-Sx+1:Jx)+rhoAExtended(dyi,dsx1);

		%top
		rhoA(1:Sy+1,:)=rhoA(1:Sy+1,:)+rhoAExtended(dsyend,dxi);
		%bottom
		rhoA(Jy-Sy+1:Jy,:)=rhoA(Jy-Sy+1:Jy,:)+rhoAExtended(dsy1,dxi);

		%top left corner
		rhoA(1:Sy+1,1:Sx+1)=rhoA(1:Sy+1,1:Sx+1)+rhoAExtended(dsyend,dsxend);
		%top right corner
		rhoA(1:Sy+1,Jx-Sx+1:Jx)=rhoA(1:Sy+1,Jx-Sx+1:Jx)+rhoAExtended(dsyend,dsx1);
		%bottom left corner
		rhoA(Jy-Sy+1:Jy,1:Sx+1)=rhoA(Jy-Sy+1:Jy,1:Sx+1)+rhoAExtended(dsy1,dsxend);
		%bottom right corner
		rhoA(Jy-Sy+1:Jy,Jx-Sx+1:Jx)=rhoA(Jy-Sy+1:Jy,Jx-Sx+1:Jx)+rhoAExtended(dsy1,dsx1);

		end

		function rhoB=bacteriadensityBperiodic(obj)
		%return bacteria B density array evaluated on grid points of domain
		coordinateArray=obj.coordinatesB();
		rhoBExtended=KDE2D(coordinateArray,obj.kernelfun,...
			obj.domainGridExtended.X,obj.domainGridExtended.Y,obj.bandwidth);

		%load index arrays
		dxi=obj.dxi;
		dsx1=obj.dsx1;
		dsxend=obj.dsxend;

		dyi=obj.dyi;
		dsy1=obj.dsy1;
		dsyend=obj.dsyend;

		Sx=obj.Sx;
		Sy=obj.Sy;

		Jx=obj.Jx;
		Jy=obj.Jy;

		%center part
		rhoB=rhoBExtended(dyi,dxi);

		%left
		rhoB(:,1:Sx+1)=rhoB(:,1:Sx+1)+rhoBExtended(dyi,dsxend);
		%right
		rhoB(:,Jx-Sx+1:Jx)=rhoB(:,Jx-Sx+1:Jx)+rhoBExtended(dyi,dsx1);

		%top
		rhoB(1:Sy+1,:)=rhoB(1:Sy+1,:)+rhoBExtended(dsyend,dxi);
		%bottom
		rhoB(Jy-Sy+1:Jy,:)=rhoB(Jy-Sy+1:Jy,:)+rhoBExtended(dsy1,dxi);

		%top left corner
		rhoB(1:Sy+1,1:Sx+1)=rhoB(1:Sy+1,1:Sx+1)+rhoBExtended(dsyend,dsxend);
		%top right corner
		rhoB(1:Sy+1,Jx-Sx+1:Jx)=rhoB(1:Sy+1,Jx-Sx+1:Jx)+rhoBExtended(dsyend,dsx1);
		%bottom left corner
		rhoB(Jy-Sy+1:Jy,1:Sx+1)=rhoB(Jy-Sy+1:Jy,1:Sx+1)+rhoBExtended(dsy1,dsxend);
		%bottom right corner
		rhoB(Jy-Sy+1:Jy,Jx-Sx+1:Jx)=rhoB(Jy-Sy+1:Jy,Jx-Sx+1:Jx)+rhoBExtended(dsy1,dsx1);
		end

		function rho=bacteriadensityperiodic(obj)
		%return total bacteria density array evaluated on grid points of domain
		coordinateArray=obj.coordinates();
		rhoExtended=KDE2D(coordinateArray,obj.kernelfun,...
			obj.domainGridExtended.X,obj.domainGridExtended.Y,obj.bandwidth);

		%load index arrays
		dxi=obj.dxi;
		dsx1=obj.dsx1;
		dsxend=obj.dsxend;

		dyi=obj.dyi;
		dsy1=obj.dsy1;
		dsyend=obj.dsyend;

		Sx=obj.Sx;
		Sy=obj.Sy;

		Jx=obj.Jx;
		Jy=obj.Jy;

		%center part
		rho=rhoExtended(dyi,dxi);

		%left
		rho(:,1:Sx+1)=rho(:,1:Sx+1)+rhoExtended(dyi,dsxend);
		%right
		rho(:,Jx-Sx+1:Jx)=rho(:,Jx-Sx+1:Jx)+rhoExtended(dyi,dsx1);

		%top
		rho(1:Sy+1,:)=rho(1:Sy+1,:)+rhoExtended(dsyend,dxi);
		%bottom
		rho(Jy-Sy+1:Jy,:)=rho(Jy-Sy+1:Jy,:)+rhoExtended(dsy1,dxi);

		%top left corner
		rho(1:Sy+1,1:Sx+1)=rho(1:Sy+1,1:Sx+1)+rhoExtended(dsyend,dsxend);
		%top right corner
		rho(1:Sy+1,Jx-Sx+1:Jx)=rho(1:Sy+1,Jx-Sx+1:Jx)+rhoExtended(dsyend,dsx1);
		%bottom left corner
		rho(Jy-Sy+1:Jy,1:Sx+1)=rho(Jy-Sy+1:Jy,1:Sx+1)+rhoExtended(dsy1,dsxend);
		%bottom right corner
		rho(Jy-Sy+1:Jy,Jx-Sx+1:Jx)=rho(Jy-Sy+1:Jy,Jx-Sx+1:Jx)+rhoExtended(dsy1,dsx1);
		end

		function refreshneighborscells(obj)
		%searches for neighbors using cells algorithm, serial implementation

		rsearch=obj.rsearch;
		rsearch2=obj.rsearch2;
		domain=obj.domain;
		XLength=domain.x(end);
		YLength=domain.y(end);
		%zero flux
		%Kx=floor(XLength/rsearch)+1;
		%periodic
		dx=domain.x(2)-domain.x(1);
		dy=domain.y(2)-domain.y(1);
		Kx=floor((XLength+dx)/rsearch)+1;
		Ky=floor((YLength+dy)/rsearch)+1;

			function cellid=determinecellid(x,y)
				%determine cell id based on coordinates

				%grid x coordinates
				xi=x-mod(x,rsearch);
				yi=y-mod(y,rsearch);

				%indices
				i=xi/rsearch+1;
				j=yi/rsearch+1;

				%cellid
				cellid=i+(j-1)*Kx;
			end

		f=@determinecellid;

		numA=obj.numA;
		numB=obj.numB;
		N=numA+numB;
		bactX=obj.bactX;
		bactY=obj.bactY;
		cellidArray=zeros(N,1);

		%determine cellid for all bacteria
		parfor i=1:N
		%for i=1:N
			x=bactX(i);
			y=bactY(i);

			cellidArray(i)=f(x,y);
		end

		%sort cellids
		[cellidArraySorted,sortedIndices]=sort(cellidArray);
		sortedBactX=bactX(sortedIndices);
		sortedBactY=bactY(sortedIndices);

		%get unique cellids
		uniqueCellidArray=unique(cellidArraySorted);
		%Kx
		n=numel(uniqueCellidArray);

		%determine starting and ending index of cell ids
		limitIndices=zeros(n,2);

		limitIndices(1,1)=1;
		lastLimitNumber=1;
		lastCellid=uniqueCellidArray(1);

		for i=2:N
			currentCellid=cellidArraySorted(i);

			if lastCellid~=currentCellid
				limitIndices(lastLimitNumber,2)=i-1;
				limitIndices(lastLimitNumber+1,1)=i;
				
				lastLimitNumber=lastLimitNumber+1;
				lastCellid=currentCellid;
			end
		end

		limitIndices(n,2)=N;

		%initialize neighbor list
		bactNb=cell(N,1);
		%periodic
		virtBactNb=cell(N,1);
		lastCellidNumber=1;
		currentCellidNumber=1;
		lastCellid=uniqueCellidArray(1);
		finalCellid=uniqueCellidArray(end);

		%iterate over bacteria, create list of potential neighbors and check for neighbors
		for i=1:N
			%i
			%initialise current cell and cell id
			indexArray=[];
			currentCellid=cellidArraySorted(i);

			%keep track of current cellid number
			if lastCellid~=currentCellid
				lastCellidNumber=lastCellidNumber+1;
				lastCellid=currentCellid;
				currentCellidNumber=lastCellidNumber;

				%supercomputer: comment this error check for improved performance!
				%if lastCellid~=uniqueCellidArray(lastCellidNumber)
				%	warning('Error in cell algorithm');
				%end
			end

			%add current cell potential neighbors
			endCurrentCellIndex=limitIndices(currentCellidNumber,2);
			indexArray=[indexArray i+1:endCurrentCellIndex];

			%check if k+1, k+Kx-1, k+Kx and k+Kx+1 cellids exists and if it exists, add to potential neighbors
			%custom method
			if currentCellid<finalCellid
				otherCellidArray=[currentCellid+1 currentCellid+Kx-1 currentCellid+Kx currentCellid+Kx+1];
			else
				otherCellidArray=[];
			end

			lastOtherCellidNumber=currentCellidNumber+1;

			for otherCellid=otherCellidArray
				while lastOtherCellidNumber<n && otherCellid>uniqueCellidArray(lastOtherCellidNumber)
					lastOtherCellidNumber=lastOtherCellidNumber+1;
				end

				if otherCellid==uniqueCellidArray(lastOtherCellidNumber)
					cellidNumber=lastOtherCellidNumber;
					startCellIndex=limitIndices(cellidNumber,1);
					endCellIndex=limitIndices(cellidNumber,2);

					indexArray=[indexArray startCellIndex:endCellIndex];
				end
			end

			%inner loop over potential neighbors
			for otherIndex=indexArray
				currentIndex=i;
				otherIndex=otherIndex;

				currentX=sortedBactX(currentIndex);
				currentY=sortedBactY(currentIndex);

				otherX=sortedBactX(otherIndex);
				otherY=sortedBactY(otherIndex);

				rsquare=(currentX-otherX)^2+(currentY-otherY)^2;

				if rsquare<rsearch2
					r=sqrt(rsquare);

					%convert sorted index to index of actual list
					currentIndex=sortedIndices(currentIndex);
					otherIndex=sortedIndices(otherIndex);

					bactNb{currentIndex}=[bactNb{currentIndex},[otherIndex;otherX;otherY;r]];
					bactNb{otherIndex}=[bactNb{otherIndex},[currentIndex;currentX;currentY;r]];
				end
			end

			%%% -- PERIODIC BOUNDARY CONDITIONS -- %%
			%%left edge neighbors
			if mod(currentCellid-1,Kx)==0		%left edge
				virtType=1;

				otherCellidArray=[currentCellid+Kx-2,...
					currentCellid+Kx-1,...
					currentCellid+2*Kx-2,...
					currentCellid+2*Kx-1];

				lastOtherCellidNumber=currentCellidNumber+1;

				for otherCellid=otherCellidArray
					while lastOtherCellidNumber<n && otherCellid>uniqueCellidArray(lastOtherCellidNumber)
						lastOtherCellidNumber=lastOtherCellidNumber+1;
					end

					if otherCellid==uniqueCellidArray(lastOtherCellidNumber)
						cellidNumber=lastOtherCellidNumber;
						startCellIndex=limitIndices(cellidNumber,1);
						endCellIndex=limitIndices(cellidNumber,2);

						indexArray=[indexArray startCellIndex:endCellIndex];
					end
				end

				%inner loop over potential neighbors
				for otherIndex=indexArray
					currentIndex=i;
					otherIndex=otherIndex;

					currentX=sortedBactX(currentIndex);
					currentY=sortedBactY(currentIndex);

					otherX=sortedBactX(otherIndex);
					otherY=sortedBactY(otherIndex);

					otherX=otherX-XLength-dx;
					otherY=otherY;

					rsquare=(currentX-otherX)^2+(currentY-otherY)^2;

					currentX=currentX+XLength+dx;
					currentY=currentY;

					if rsquare<rsearch2
						r=sqrt(rsquare);

						%convert sorted index to index of actual list
						currentIndex=sortedIndices(currentIndex);
						otherIndex=sortedIndices(otherIndex);

						virtBactNb{currentIndex}=[virtBactNb{currentIndex},...
							[otherIndex;otherX;otherY;r;virtType]];
						virtBactNb{otherIndex}=[virtBactNb{otherIndex},...
							[currentIndex;currentX;currentY;r;virtType+1]];
					end
				end
			end
			
			%%top edge neighbors
			if currentCellid<=Kx
				virtType=3;

				otherCellidArray=[currentCellid+(Ky-2)*Kx,...
					currentCellid+(Ky-2)*Kx+1,...
					currentCellid+(Ky-1)*Kx,...
					currentCellid+(Ky-1)*Kx+1];

				lastOtherCellidNumber=currentCellidNumber+1;

				for otherCellid=otherCellidArray
					while lastOtherCellidNumber<n && otherCellid>uniqueCellidArray(lastOtherCellidNumber)
						lastOtherCellidNumber=lastOtherCellidNumber+1;
					end

					if otherCellid==uniqueCellidArray(lastOtherCellidNumber)
						cellidNumber=lastOtherCellidNumber;
						startCellIndex=limitIndices(cellidNumber,1);
						endCellIndex=limitIndices(cellidNumber,2);

						indexArray=[indexArray startCellIndex:endCellIndex];
					end
				end

				%inner loop over potential neighbors
				for otherIndex=indexArray
					currentIndex=i;
					otherIndex=otherIndex;

					currentX=sortedBactX(currentIndex);
					currentY=sortedBactY(currentIndex);

					otherX=sortedBactX(otherIndex);
					otherY=sortedBactY(otherIndex);

					otherX=otherX;
					otherY=otherY-YLength-dy;

					rsquare=(currentX-otherX)^2+(currentY-otherY)^2;

					currentX=currentX;
					currentY=currentY+YLength+dy;

					if rsquare<rsearch2
						r=sqrt(rsquare);

						%convert sorted index to index of actual list
						currentIndex=sortedIndices(currentIndex);
						otherIndex=sortedIndices(otherIndex);

						virtBactNb{currentIndex}=[virtBactNb{currentIndex},...
							[otherIndex;otherX;otherY;r;virtType]];
						virtBactNb{otherIndex}=[virtBactNb{otherIndex},...
							[currentIndex;currentX;currentY;r;virtType+1]];
					end
				end
			end

			%%topleft neighbors
			if currentCellid==1
				virtType=5;

				otherCellidArray=[(Ky-1)*Kx-1,...
					(Ky-1)*Kx,...
					Ky*Kx-1,...
					Ky*Kx];

				lastOtherCellidNumber=currentCellidNumber+1;

				for otherCellid=otherCellidArray
					while lastOtherCellidNumber<n && otherCellid>uniqueCellidArray(lastOtherCellidNumber)
						lastOtherCellidNumber=lastOtherCellidNumber+1;
					end

					if otherCellid==uniqueCellidArray(lastOtherCellidNumber)
						cellidNumber=lastOtherCellidNumber;
						startCellIndex=limitIndices(cellidNumber,1);
						endCellIndex=limitIndices(cellidNumber,2);

						indexArray=[indexArray startCellIndex:endCellIndex];
					end
				end

				%inner loop over potential neighbors
				for otherIndex=indexArray
					currentIndex=i;
					otherIndex=otherIndex;

					currentX=sortedBactX(currentIndex);
					currentY=sortedBactY(currentIndex);

					otherX=sortedBactX(otherIndex);
					otherY=sortedBactY(otherIndex);

					otherX=otherX-XLength-dx;
					otherY=otherY-YLength-dy;

					rsquare=(currentX-otherX)^2+(currentY-otherY)^2;

					currentX=currentX+XLength+dx;
					currentY=currentY+YLength+dy;

					if rsquare<rsearch2
						r=sqrt(rsquare);

						%convert sorted index to index of actual list
						currentIndex=sortedIndices(currentIndex);
						otherIndex=sortedIndices(otherIndex);

						virtBactNb{currentIndex}=[virtBactNb{currentIndex},...
							[otherIndex;otherX;otherY;r;virtType]];
						virtBactNb{otherIndex}=[virtBactNb{otherIndex},...
							[currentIndex;currentX;currentY;r;virtType+1]];
					end
				end
			end

			%%topright neighbors
			if currentCellid==Kx||currentCellid==Kx-1
				virtType=7;

				otherCellidArray=[(Ky-2)*Kx+1,...
					(Ky-1)*Kx+1];

				lastOtherCellidNumber=currentCellidNumber+1;

				for otherCellid=otherCellidArray
					while lastOtherCellidNumber<n && otherCellid>uniqueCellidArray(lastOtherCellidNumber)
						lastOtherCellidNumber=lastOtherCellidNumber+1;
					end

					if otherCellid==uniqueCellidArray(lastOtherCellidNumber)
						cellidNumber=lastOtherCellidNumber;
						startCellIndex=limitIndices(cellidNumber,1);
						endCellIndex=limitIndices(cellidNumber,2);

						indexArray=[indexArray startCellIndex:endCellIndex];
					end
				end

				%inner loop over potential neighbors
				for otherIndex=indexArray
					currentIndex=i;
					otherIndex=otherIndex;

					currentX=sortedBactX(currentIndex);
					currentY=sortedBactY(currentIndex);

					otherX=sortedBactX(otherIndex);
					otherY=sortedBactY(otherIndex);

					otherX=otherX+XLength+dx;
					otherY=otherY-YLength-dy;

					rsquare=(currentX-otherX)^2+(currentY-otherY)^2;

					currentX=currentX-XLength-dx;
					currentY=currentY+YLength+dy;

					if rsquare<rsearch2
						r=sqrt(rsquare);

						%convert sorted index to index of actual list
						currentIndex=sortedIndices(currentIndex);
						otherIndex=sortedIndices(otherIndex);

						virtBactNb{currentIndex}=[virtBactNb{currentIndex},...
							[otherIndex;otherX;otherY;r;virtType]];
						virtBactNb{otherIndex}=[virtBactNb{otherIndex},...
							[currentIndex;currentX;currentY;r;virtType+1]];
					end
				end
			end
		end
		obj.bactNb=bactNb;
		obj.virtBactNb=virtBactNb;
		end

		function refreshneighborsprojection(obj)
		%searches for neighbors using projection algorithm, serial implementation
		bactX=obj.bactX;
		bactY=obj.bactY;
		varX=var(bactX);
		varY=var(bactY);

		%select dimension with largest variance
		if varX>varY
			[coordArray,indices]=sort(bactX);
		else
			[coordArray,indices]=sort(bactY);
		end

		%initialize neighbor list
		N=obj.numA+obj.numB;
		bactNb=cell(N,1);

		%initialize pointers
		i=1;
		j=1;

		rsearch=obj.rsearch;
		rsearch2=obj.rsearch2;

		%outer loop
		while j<=N
			currentCoord=coordArray(i);

			%set second pointer
			while coordArray(j)-currentCoord<rsearch && j<N
				j=j+1;
			end

			%inner loop over subset of potential neighbors
			for k=i+1:j
				currentIndex=indices(i);
				otherIndex=indices(k);

				currentX=bactX(currentIndex);
				currentY=bactY(currentIndex);

				otherX=bactX(otherIndex);
				otherY=bactY(otherIndex);

				rsquare=(currentX-otherX)^2+(currentY-otherY)^2;

				if rsquare<rsearch2
					r=sqrt(rsquare);

					bactNb{currentIndex}=[bactNb{currentIndex},[otherIndex;otherX;otherY;r]];
					bactNb{otherIndex}=[bactNb{otherIndex},[currentIndex;currentX;currentY;r]];
				end
			end

			i=i+1;
			if j==N
				j=j+1;
			end
		end
		obj.bactNb=bactNb;
		end

		function updateneighbors(obj)
		%update coordinates of neighbors and distance

		numA=obj.numA;
		numB=obj.numB;
		N=numA+numB;

		bactX=obj.bactX;
		bactY=obj.bactY;
		bactNb=obj.bactNb;
		virtBactNb=obj.virtBactNb;

		domain=obj.domain;
		XLength=domain.x(end);
		YLength=domain.y(end);
		dx=domain.x(2)-domain.x(1);
		dy=domain.y(2)-domain.y(1);
		Xdx=XLength+dx;
		Ydy=YLength+dy;

		%parallel or serial? faster for large bateria population
		parfor i=1:N
		%for i=1:N
			currentX=bactX(i);
			currentY=bactY(i);

			A=bactNb{i};
			[~,m]=size(A);

			for j=1:m
				otherIndex=A(1,j);
				otherX=bactX(otherIndex);
				otherY=bactY(otherIndex);

				r=sqrt((currentX-otherX)^2+(currentY-otherY)^2);

				A(2:4,j)=[otherX;otherY;r];
			end
			bactNb{i}=A;

			A=virtBactNb{i};
			[~,m]=size(A);

			for j=1:m
				otherIndex=A(1,j);
				otherX=bactX(otherIndex);
				otherY=bactY(otherIndex);

				virtType=A(5,j);

				switch virtType
				case 1	%left edge
					otherX=otherX-Xdx;
					otherY=otherY;
				case 2	%right edge
					otherX=otherX+Xdx;
					otherY=otherY;
				case 3	%top edge
					otherX=otherX;
					otherY=otherY-Ydy;
				case 4	%bottom edge
					otherX=otherX;
					otherY=otherY+Ydy;
				case 5	%top left corner
					otherX=otherX-Xdx;
					otherY=otherY-Ydy;
				case 6	%bottom right corner
					otherX=otherX+Xdx;
					otherY=otherY+Ydy;
				case 7	%top right corner
					otherX=otherX+Xdx;
					otherY=otherY-Ydy;
				case 8	%bottom left corner
					otherX=otherX-Xdx;
					otherY=otherY+Ydy;
				otherwise
					error('Unkown edge case');
				end

				r=sqrt((currentX-otherX)^2+(currentY-otherY)^2);

				A(2:4,j)=[otherX;otherY;r];
			end
			virtBactNb{i}=A;
		end
		obj.bactNb=bactNb;
		obj.virtBactNb=virtBactNb;
		end

		function refreshorupdateneighbors(obj)
		%either refresh or update neighbors

		%update neighbors
		if mod(obj.counter,obj.modulo)==0
			%obj.refreshneighborsprojection();
			obj.refreshneighborscells();
		else
			obj.updateneighbors();
		end

		obj.counter=obj.counter+1;
		end

		function [currentMuArray]=determinemu(obj,AHLField)
		%calculate mobilities for every bacterium

		domain=obj.domain;
			function interpolatedValue=interpol(field,coordinateVector)
				%disp(['# of zeros in field: ' num2str(sum(sum(field==0)))]);
				xCoordinate=coordinateVector(1);
				yCoordinate=coordinateVector(2);
				m=length(domain.y);
				n=length(domain.x);

				dx=domain.x(2)-domain.x(1);
				dy=domain.y(2)-domain.y(1);

				%search for i
				yDiff=yCoordinate-domain.y(1);
				iDiff=yDiff/dy;
				i=floor(iDiff)+1;
				iplus=i+1;
				if i==m
					%zero flux boundary condition
					%i=m-1;
					%periodic boundary condition
					iplus=1;
				end

				%search for j
				xDiff=xCoordinate-domain.x(1);
				jDiff=xDiff/dx;
				j=floor(jDiff)+1;
				jplus=j+1;
				if j==n
					%zero flux boundary condition
					%j=n-1;
					%periodic boundary condition
					jplus=1;
				end

				%delta x and delta y
				deltaY=yCoordinate-domain.y(i);
				deltaX=xCoordinate-domain.x(j);
				sigmaY=1-deltaY;
				sigmaX=1-deltaX;

				%disp(['i: ' num2str(i)]);
				%disp(['j: ' num2str(j)]);
				%disp(['i+1: ' num2str(i+1)]);
				%disp(['j+1): ' num2str(j+1)]);
				%disp(['field(i,j): ' num2str(field(i,j))]);
				%disp(['field(i+1,j): ' num2str(field(i+1,j))]);
				%disp(['field(i,j+1): ' num2str(field(i,j+1))]);
				%disp(['field(i+1,j+1): ' num2str(field(i+1,j+1))]);
				interpolatedValue=1/(dx*dy)*(sigmaX*sigmaY*field(i,j)+...
					sigmaY*deltaX*field(i,jplus)+sigmaX*deltaY*field(iplus,j)+...
					deltaX*deltaY*field(iplus,jplus));
			end
		f=@interpol;

		conc=AHLField.getconcentration();

		numA=obj.numA;
		numB=obj.numB;
		N=numA+numB;
		AHL=zeros(N,1);
		bactX=obj.bactX;
		bactY=obj.bactY;

		%parallel or serial?
		%parfor k=1:N
		for k=1:N
			currentX=bactX(k);
			currentY=bactY(k);

			AHL(k)=f(conc,[currentX currentY]);
		end

		muA=obj.muA;
		VthA=obj.VthA;
		muB=obj.muB;
		VthB=obj.VthB;
		currentMuArray=zeros(N,1);

		%parallel or serial?
		%parfor i=1:N
		%for i=1:N
		%	%bacteria A
		%	if i<=numA
		%		if AHL(i)>VthA
		%			currentMuArray(i)=muA.low;
		%		else
		%			currentMuArray(i)=muA.high;
		%		end
		%	%bacteria B
		%	else
		%		if AHL(i)<VthB
		%			currentMuArray(i)=muB.low;
		%		else
		%			currentMuArray(i)=muB.high;
		%		end
		%	end
		%end

		%alternative
		currentMuAArray=(AHL(1:numA)>VthA)*muA.low+(AHL(1:numA)<=VthA)*muA.high;
		currentMuBArray=(AHL(numA+1:end)<VthB)*muB.low+(AHL(numA+1:end)>=VthB)*muB.high;
		currentMuArray=[currentMuAArray;currentMuBArray];
		end

		function [chemodx,chemody]=calculatechemo(obj,leucineField,currentMuArray,dt)
		%calculate displacement due to chemotaxis

		domain=obj.domain;
			function interpolatedValue=interpol(field,coordinateVector)
				%disp(['# of zeros in field: ' num2str(sum(sum(field==0)))]);
				xCoordinate=coordinateVector(1);
				yCoordinate=coordinateVector(2);
				m=length(domain.y);
				n=length(domain.x);

				dx=domain.x(2)-domain.x(1);
				dy=domain.y(2)-domain.y(1);

				%search for i
				yDiff=yCoordinate-domain.y(1);
				iDiff=yDiff/dy;
				i=floor(iDiff)+1;
				iplus=i+1;
				if i==m
					%zero flux boundary condition
					%i=m-1;
					%periodic boundary condition
					iplus=1;
				end

				%search for j
				xDiff=xCoordinate-domain.x(1);
				jDiff=xDiff/dx;
				j=floor(jDiff)+1;
				jplus=j+1;
				if j==n
					%zero flux boundary condition
					%j=n-1;
					%periodic boundary condition
					jplus=1;
				end

				%delta x and delta y
				deltaY=yCoordinate-domain.y(i);
				deltaX=xCoordinate-domain.x(j);
				sigmaY=1-deltaY;
				sigmaX=1-deltaX;

				%disp(['i: ' num2str(i)]);
				%disp(['j: ' num2str(j)]);
				%disp(['i+1: ' num2str(i+1)]);
				%disp(['j+1): ' num2str(j+1)]);
				%disp(['field(i,j): ' num2str(field(i,j))]);
				%disp(['field(i+1,j): ' num2str(field(i+1,j))]);
				%disp(['field(i,j+1): ' num2str(field(i,j+1))]);
				%disp(['field(i+1,j+1): ' num2str(field(i+1,j+1))]);
				interpolatedValue=1/(dx*dy)*(sigmaX*sigmaY*field(i,j)+...
					sigmaY*deltaX*field(i,jplus)+sigmaX*deltaY*field(iplus,j)+...
					deltaX*deltaY*field(iplus,jplus));
			end
		f=@interpol;

		numA=obj.numA;
		numB=obj.numB;
		N=numA+numB;

		xgradient=leucineField.gradient.x;
		ygradient=leucineField.gradient.y;
		concentration=leucineField.concentration;
		%disp(['# of zeros in concentration: ' num2str(sum(sum(concentration==0)))]);

		dleucinex=zeros(N,1);
		dleuciney=zeros(N,1);
		leucine=zeros(N,1);

		bactX=obj.bactX;
		bactY=obj.bactY;

		%parallel or serial?
		%parfor i=numA+1:N
		for k=numA+1:N
			currentX=bactX(k);
			currentY=bactY(k);

			dleucinex(k)=f(xgradient,[currentX currentY]);
			dleuciney(k)=f(ygradient,[currentX currentY]);
			leucine(k)=f(concentration,[currentX currentY]);
			%pause;
			%disp(['assigned leucine value in for: ' num2str(temp)]);
		end

		%leucine

		kappaB=obj.kappaB;
		chemodx=zeros(N,1);
		chemody=zeros(N,1);

		%for i=numA+1:N
		%	chemodx(i)=-currentMuArray(i)*kappaB/leucine(i)*dleucinex(i)*dt;
		%	chemody(i)=-currentMuArray(i)*kappaB/leucine(i)*dleuciney(i)*dt;
		%end

		%disp(['# of NaN numbers in leucine: ' num2str(sum(isnan(leucine)))]);
		%disp(['# of zeros in leucine: ' num2str(sum(leucine==0))]);
		%disp(['# of zeros in leucine(numA+1:N): ' num2str(sum(leucine(numA+1:N)==0))]);
		%disp(['# of NaN numbers in currentMuArray: ' num2str(sum(isnan(currentMuArray)))]);
		%disp(['# of NaN numbers in dleucinex: ' num2str(sum(isnan(dleucinex)))]);
		%disp(['# of NaN numbers in dleuciney: ' num2str(sum(isnan(dleuciney)))]);
		%disp(['# of NaN numbers in kappaB: ' num2str(sum(isnan(kappaB)))]);
		%disp(['# of NaN numbers in dt: ' num2str(sum(isnan(dt)))]);
		%disp(['# of NaN numbers in dt: ' num2str(sum(isnan(dt)))]);

		chemodx(numA+1:N)=(-1).*currentMuArray(numA+1:N).*kappaB./leucine(numA+1:N).*dleucinex(numA+1:N).*dt;
		chemody(numA+1:N)=(-1).*currentMuArray(numA+1:N).*kappaB./leucine(numA+1:N).*dleuciney(numA+1:N).*dt;
		end

		function [randdx,randdy]=calculaterand(obj,currentMuArray,dt)
		%calculate displacement due to brownian movement

		numA=obj.numA;
		numB=obj.numB;
		N=numA+numB;
		randdx=zeros(N,1);
		randdy=zeros(N,1);

		%parallel or serial?
		%parfor i=1:N
		%for i=1:N
		%	randdx(i)=sqrt(2*currentMuArray(i)*dt)*normrnd(0,1);
		%	randdy(i)=sqrt(2*currentMuArray(i)*dt)*normrnd(0,1);
		%end

		%faster?
		randdx=sqrt(2*currentMuArray*dt).*normrnd(0,1,N,1);
		randdy=sqrt(2*currentMuArray*dt).*normrnd(0,1,N,1);
		end
		
		function [celldx,celldy]=calculatecell(obj,dt)
		%calculate displacement due to interaction with neighboring cells

		numA=obj.numA;
		numB=obj.numB;
		N=numA+numB;

		k1=obj.k1;
		k2=obj.k2;
		r0=obj.r0;
		gamma=obj.gamma;

			function [dx,dy]=cellforce(x1,y1,x2,y2,r,rcut,k3)
			%Calculates displacement of bacterium 1 due to bacterium 2

				%out of influence
				if r>=rcut*2
					dx=0;
					dy=0;
					return
				end

				%calculate unit vector pointing away from bacterium 2
				if r==0
					theta=rand*2*pi;
					ex=cos(theta);
					ey=sin(theta);
				else
					ex=(x1-x2)/r;
					ey=(y1-y2)/r;
				end

				%calculate force
				if r>=r0*2
					Fr=k3*(2*r0-r);
				elseif r>=r0
					Fr=k2*(2*r0-r);
				else
					Fr=k1*((k1+k2)/k1*r0-r);
				end

				%calculate displacement
				dx=dt/gamma*Fr*ex;
				dy=dt/gamma*Fr*ey;
			end
		f=@cellforce;

		rcutA=obj.rcut;
		k3a=obj.k3;
		bactNb=obj.bactNb;
		bactX=obj.bactX;
		bactY=obj.bactY;

		celldx=zeros(N,1);
		celldy=zeros(N,1);

		bactNb=obj.bactNb;
		virtBactNb=obj.virtBactNb;

		%parallel or serial? -> parallel
		parfor i=1:N
		%for i=1:N
			%current coordinates
			currentX=bactX(i);
			currentY=bactY(i);

			%load neighbor matrix
			nbMatrix=bactNb{i};
			[~,m1]=size(nbMatrix);

			%load virtual neighbor matrix
			virtNbMatrix=virtBactNb{i};
			[~,m2]=size(virtNbMatrix);

			%total neighbor matrix
			nbMatrix=[nbMatrix,virtNbMatrix(1:end-1,:)];
			m=m1+m2;

			%initialise total displacement
			totaldx=0;
			totaldy=0;

			%iterate over neighbors
			for j=1:m
				otherIndex=nbMatrix(1,j);
				otherX=nbMatrix(2,j);
				otherY=nbMatrix(3,j);
				r=nbMatrix(4,j);

				%distinguish between A <-> A interactions and other interactions
				if i<=numA && otherIndex<=numA
					rcut=rcutA;
					k3=k3a;
				else
					rcut=r0;
					k3=0;
				end

				[dx,dy]=f(currentX,currentY,otherX,otherY,r,rcut,k3);

				totaldx=totaldx+dx;
				totaldy=totaldy+dy;
			end
			
			celldx(i)=totaldx;
			celldy(i)=totaldy;
		end
		end

		function updatefast(obj,AHLField,leucineField,dt)
		%update based on AHL field, leuceine field and timestep

		%refresh or update neighbors 
		obj.refreshorupdateneighbors();

		%determine mu for every bacterium
		%currentMuArray=obj.determinemu(AHLField);

		%calculate displacement due to chemotaxis, neighboring cell interactions and brownian motion
		%[chemodx,chemody]=calculatechemo(obj,leucineField,currentMuArray,dt);
		chemodx=0;
		chemody=0;
		%[randdx,randdy]=calculaterand(obj,currentMuArray,dt);
		randdx=0;
		randdy=0;
		[celldx,celldy]=calculatecell(obj,dt);
		%celldx=0;
		%celldy=0;

		%test periodic boundary conditions
		%vx=1;
		%vy=0;
		%vx=0;
		%vy=0;
		%testx=vx*dt;
		%testy=vy*dt;

		%disp(['# of NaN numbers in chemodx: ' num2str(sum(isnan(chemodx)))]);
		%disp(['# of NaN numbers in chemody: ' num2str(sum(isnan(chemody)))]);
		%disp(['# of NaN numbers in randdx: ' num2str(sum(isnan(randdx)))]);
		%disp(['# of NaN numbers in randdy: ' num2str(sum(isnan(randdy)))]);
		%disp(['# of NaN numbers in celldx: ' num2str(sum(isnan(celldx)))]);
		%disp(['# of NaN numbers in celldy: ' num2str(sum(isnan(celldy)))]);
		%calculate new position
		newBactX=obj.bactX+chemodx+randdx+celldx;
		newBactY=obj.bactY+chemody+randdy+celldy;

		%newBactX=obj.bactX+testx;
		%newBactY=obj.bactY+testy;

		%correct for going out of boundary
		N=obj.numA+obj.numB;

		domain=obj.domain;
		dx=domain.x(2)-domain.x(1);
		domainX1=domain.x(1);
		domainXEnd=domain.x(end);
		domainXEndPeriodic=domainXEnd+dx;
		dy=domain.y(2)-domain.y(1);
		domainY1=domain.y(1);
		domainYEnd=domain.y(end);
		domainYEndPeriodic=domainYEnd+dy;

		%%zero flux boundary condition
		%%parallel or serial?
		%%parfor i=1:N
		%for i=1:N
		%	currentX=newBactX(i);
		%	currentY=newBactY(i);

		%	if currentX<domainX1
		%		newBactX(i)=domainX1;
		%	elseif currentX>domainXEnd
		%		newBactX(i)=domainXEnd;
		%	else
		%		newBactX(i)=currentX;
		%	end

		%	if currentY<domainY1
		%		newBactY(i)=domainY1;
		%	elseif currentY>domainYEnd
		%		newBactY(i)=domainYEnd;
		%	else
		%		newBactY(i)=currentY;
		%	end
		%end

		%Periodic boundary condition
		%parallel or serial?
		%parfor i=1:N
		for i=1:N
			currentX=newBactX(i);
			currentY=newBactY(i);

			newBactX(i)=mod(currentX,domainXEndPeriodic);
			newBactY(i)=mod(currentY,domainYEndPeriodic);
		end

		obj.bactX=newBactX;
		obj.bactY=newBactY;

		%disp(['# of NaN numbers in x: ' num2str(sum(isnan(obj.bactX)))]);
		%disp(['# of NaN numbers in y: ' num2str(sum(isnan(obj.bactY)))]);
		end
	end
end
