classdef bacterium < handle
	properties
		xCoordinate;
		yCoordinate;

		nbArray;
		nbXArray;
		nbYArray;
		nbIsBlackArray;
		nbRArray;
	end

	methods
		function obj=bacterium(xCoordinate,yCoordinate)
		obj.xCoordinate=xCoordinate;
		obj.yCoordinate=yCoordinate;
		end

		function xCoordinate=getxcoordinate(obj)
		xCoordinate=obj.xCoordinate;
		end

		function setxcoordinate(obj,xCoordinate)
		obj.xCoordinate=xCoordinate;
		end

		function yCoordinate=getycoordinate(obj)
		yCoordinate=obj.yCoordinate;
		end

		function setycoordinate(obj,yCoordinate)
		obj.yCoordinate=yCoordinate;
		end

		function addneighbor(obj,neighbor,xCoordinate,yCoordinate,isBlack,r)
		obj.nbArray=[obj.nbArray neighbor];
		obj.nbXArray=[obj.nbXArray xCoordinate];
		obj.nbYArray=[obj.nbYArray yCoordinate];
		obj.nbIsBlackArray=[obj.nbIsBlackArray isBlack];
		obj.nbRArray=[obj.nbRArray r];
		end

		function addneighborarray(obj,nbArray,nbXArray,nbYArray,nbIsBlackArray,nbRArray)
		n=length(nbArray);
		for i=1:n
			obj.addneighbor(nbArray(i),nbXArray(i),nbYArray(i),nbIsBlackArray(i),nbRArray(i));
		end
		end

		function deleteithneighbor(obj,i)
		obj.nbArray(i)=[];
		obj.nbXArray(i)=[];
		obj.nbYArray(i)=[];
		obj.nbIsBlackArray(i)=[];
		obj.nbRArray(i)=[];
		end

		function deleteneighbor(obj,neighborToDelete)
		n=length(obj.nbArray);
		i=1;
		currentNeighbor=obj.nbArray(i);

		while currentNeighbor~=neighborToDelete && i<=n
			i=i+1;
			currentNeighbor=obj.nbArray(i);
		end

		if i<= n
			obj.deleteithneighbor(i);
		else
			disp('Error: neighbor to delete not found');
		end
		end

		function [nbArray,nbXArray,nbYArray,nbIsBlackArray,nbRArray]=getneighborarray(obj)
			nbArray=obj.nbArray;
			nbXArray=obj.nbXArray;
			nbYArray=obj.nbYArray;
			nbIsBlackArray=obj.nbIsBlackArray;
			nbRArray=obj.nbRArray;
		end

		function purge(obj,currentIsBlack)
		n=length(obj.nbArray);
		%disp('length nbArray');
		%length(obj.nbArray)
		%disp('length nbIsBlackArray');
		%length(obj.nbIsBlackArray)

		for i=n:-1:1
			if obj.nbIsBlackArray(i)~=currentIsBlack
				obj.deleteithneighbor(i);
			end
		end

		end

		function updateneighborscoordinates(obj,nbXArray,nbYArray,nbRArray)
		obj.nbXArray=nbXArray;
		obj.nbYArray=nbYArray;
		obj.nbRArray=nbRArray;
		end

		function [nbXArray,nbYArray,nbRArray]=getneighborscoordinates(obj)
		nbXArray=obj.nbXArray;
		nbYArray=obj.nbYArray;
		nbRArray=obj.nbRArray;
		end
				
	end
end
