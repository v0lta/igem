classdef bacteriumB < handle
	properties
		xCoordinate;
		yCoordinate;

		nbAArray;
		nbAXArray;
		nbAYArray;
		nbAIsBlackArray;
		nbARArray;

		nbBArray;
		nbBXArray;
		nbBYArray;
		nbBIsBlackArray;
		nbBRArray;
	end

	methods
		function obj=bacteriumB(xCoordinate,yCoordinate)
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

		function addneighbora(obj,neighbor,xCoordinate,yCoordinate,isBlack,r)
		obj.nbAArray=[obj.nbAArray neighbor];
		obj.nbAXArray=[obj.nbAXArray xCoordinate];
		obj.nbAYArray=[obj.nbAYArray yCoordinate];
		obj.nbAIsBlackArray=[obj.nbAIsBlackArray isBlack];
		obj.nbARArray=[obj.nbARArray r];
		end

		function addneighboraarray(obj,nbArray,nbXArray,nbYArray,nbIsBlackArray,nbRArray)
		n=length(nbArray);
		for i=1:n
			obj.addneighbora(nbArray(i),nbXArray(i),nbYArray(i),nbIsBlackArray(i),nbRArray(i));
		end
		end

		function deleteithneighbora(obj,i)
		obj.nbAArray(i)=[];
		obj.nbAXArray(i)=[];
		obj.nbAYArray(i)=[];
		obj.nbAIsBlackArray(i)=[];
		obj.nbARArray(i)=[];
		end

		function deleteneighbora(obj,neighborToDelete)
		n=length(obj.nbAArray);
		i=1;
		currentNeighbor=obj.nbAArray(i);

		while currentNeighbor~=neighborToDelete && i<=n
			i=i+1;
			currentNeighbor=obj.nbAArray(i);
		end

		if i<= n
			obj.deleteithneighbora(i);
		else
			disp('Error: neighbor to delete not found');
		end
		end

		function [nbAArray,nbAXArray,nbAYArray,nbAIsBlackArray,nbARArray]=getneighboraarray(obj)
			nbAArray=obj.nbAArray;
			nbAXArray=obj.nbAXArray;
			nbAYArray=obj.nbAYArray;
			nbAIsBlackArray=obj.nbAIsBlackArray;
			nbARArray=obj.nbARArray;
		end

		function updateneighborsacoordinates(obj,nbXArray,nbYArray,nbRArray)
		obj.nbAXArray=nbXArray;
		obj.nbAYArray=nbYArray;
		obj.nbARArray=nbRArray;
		end

		function [nbXArray,nbYArray,nbRArray]=getneighborsacoordinates(obj)
		nbXArray=obj.nbAXArray;
		nbYArray=obj.nbAYArray;
		nbRArray=obj.nbARArray;
		end

		function addneighborb(obj,neighbor,xCoordinate,yCoordinate,isBlack,r)
		obj.nbBArray=[obj.nbBArray neighbor];
		obj.nbBXArray=[obj.nbBXArray xCoordinate];
		obj.nbBYArray=[obj.nbBYArray yCoordinate];
		obj.nbBIsBlackArray=[obj.nbBIsBlackArray isBlack];
		obj.nbBRArray=[obj.nbBRArray r];
		end

		function addneighborbarray(obj,nbArray,nbXArray,nbYArray,nbIsBlackArray,nbRArray)
		n=length(nbArray);
		for i=1:n
			obj.addneighborb(nbArray(i),nbXArray(i),nbYArray(i),nbIsBlackArray(i),nbRArray(i));
		end
		end

		function deleteithneighborb(obj,i)
		obj.nbBArray(i)=[];
		obj.nbBXArray(i)=[];
		obj.nbBYArray(i)=[];
		obj.nbBIsBlackArray(i)=[];
		obj.nbBRArray(i)=[];
		end

		function deleteneighborb(obj,neighborToDelete)
		n=length(obj.nbBArray);
		i=1;
		currentNeighbor=obj.nbBArray(i);

		while currentNeighbor~=neighborToDelete && i<=n
			i=i+1;
			currentNeighbor=obj.nbBArray(i);
		end

		if i<= n
			obj.deleteithneighborb(i);
		else
			disp('Error: neighbor to delete not found');
		end
		end

		function [nbBArray,nbBXArray,nbBYArray,nbBIsBlackArray,nbBRArray]=getneighborbarray(obj)
			nbBArray=obj.nbBArray;
			nbBXArray=obj.nbBXArray;
			nbBYArray=obj.nbBYArray;
			nbBIsBlackArray=obj.nbBIsBlackArray;
			nbBRArray=obj.nbBRArray;
		end

		function updateneighborsbcoordinates(obj,nbXArray,nbYArray,nbRArray)
		obj.nbBXArray=nbXArray;
		obj.nbBYArray=nbYArray;
		obj.nbBRArray=nbRArray;
		end

		function [nbXArray,nbYArray,nbRArray]=getneighborsbcoordinates(obj)
		nbXArray=obj.nbBXArray;
		nbYArray=obj.nbBYArray;
		nbRArray=obj.nbBRArray;
		end

		function purge(obj,currentIsBlack)
		n=length(obj.nbAArray);
		m=length(obj.nbBArray);
		%disp('length nbArray');
		%length(obj.nbArray)
		%disp('length nbIsBlackArray');
		%length(obj.nbIsBlackArray)

		%purge cells A
		for i=n:-1:1
			if obj.nbAIsBlackArray(i)~=currentIsBlack
				obj.deleteithneighbora(i);
			end
		end

		%purge cells B
		for i=m:-1:1
			if obj.nbBIsBlackArray(i)~=currentIsBlack
				obj.deleteithneighborb(i);
			end
		end

		end
	end
end
