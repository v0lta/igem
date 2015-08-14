classdef bacterium < handle
	properties
		xCoordinate;
		yCoordinate;
		direction;
	end

	methods
		function obj=bacterium(xCoordinate,yCoordinate,direction)
		obj.xCoordinate=xCoordinate;
		obj.yCoordinate=yCoordinate;
		obj.direction=direction;
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

		function direction=getdirection(obj)
		direction=obj.direction;
		end

		function setdirection(obj,direction)
		obj.direction=direction;
		end
	end
end
