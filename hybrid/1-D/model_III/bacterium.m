%Model II
classdef bacterium < handle
	properties
		xCoordinate;
		%path;
		direction;	%0 for moving to left, 1 for moving to right
		memory;
	end

	methods
		function obj=bacterium(xCoordinate)
		obj.xCoordinate=xCoordinate;
		%initialize with random direction
		f=rand;
		if f<0.5
			obj.direction=0;
		else
			obj.direction=1;
		end
		%initialize without memory
		%obj.memory=0;
		obj.memory=0.1;
		end

		function xCoordinate=getxcoordinate(obj)
		xCoordinate=obj.xCoordinate;
		end

		function setxcoordinate(obj,xCoordinate)
		obj.xCoordinate=xCoordinate;
		end

		function direction=getdirection(obj)
		direction=obj.direction;
		end

		function setdirection(obj,direction)
		obj.direction=direction;
		end

		function memory=getmemory(obj)
		memory=obj.memory;
		end

		function setmemory(obj,memory)
		obj.memory=memory;
		end

	end
end
