classdef bacterium < handle
	properties
		xCoordinate;
		path;
	end

	methods
		function obj=bacterium(xCoordinate)
		obj.xCoordinate=xCoordinate;
		end

		function xCoordinate=getxcoordinate(obj)
		xCoordinate=obj.xCoordinate;
		end

		function setxcoordinate(obj,xCoordinate)
		obj.xCoordinate=xCoordinate;
		end

	end
end
