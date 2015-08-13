classdef bacterium < handle
	properties
		xCoordinate;
		yCoordinate;
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
	end
end



