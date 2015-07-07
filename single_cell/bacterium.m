classdef bacterium < handle
	properties
		xCoordinate;
		yCoordinate;
		%1 for run, 0 for tumble
		runOrTumble;
		angle;
		runCounter;
		path;
		highCounter=0;
	end
	
	properties (Constant)
		MAX_RUN_COUNTER=10;
		%10 µm per step size
		STEP_SIZE=1000;
		%STEP_SIZE=10;
	end

	methods
		function obj=bacterium(xCoordinate,yCoordinate)
			obj.xCoordinate=xCoordinate;
			obj.yCoordinate=yCoordinate;
			obj.angle=rand*2*pi;
			obj.runCounter=obj.MAX_RUN_COUNTER;
			path=[xCoordinate,yCoordinate];
		end

		function resetToTumble(obj)
			obj.runOrTumble=0;
			obj.runCounter=obj.MAX_RUN_COUNTER;
		end

		function runOnce(obj)
			%Randomize at start of run?
			%if obj.runCounter==obj.MAX_RUN_COUNTER
			%	obj.randomizeAngle();
			%end

			obj.xCoordinate=obj.xCoordinate+obj.STEP_SIZE*cos(obj.angle);
			obj.yCoordinate=obj.yCoordinate+obj.STEP_SIZE*sin(obj.angle);
			obj.runCounter=obj.runCounter-1;
			if obj.runCounter==0
				%disp('test')
				obj.resetToTumble();
			end
		end

		function randomizeAngle(obj)
			obj.angle=rand*2*pi;
		end

		function tumbleOnce(obj)
			obj.randomizeAngle();
			%disp('pretumble');
			%obj.xCoordinate
			%obj.yCoordinate
			obj.xCoordinate=obj.xCoordinate+obj.STEP_SIZE*cos(obj.angle);
			obj.yCoordinate=obj.yCoordinate+obj.STEP_SIZE*sin(obj.angle);
			%disp('posttumble');
			%obj.xCoordinate
			%obj.yCoordinate
		end

		function runOrTumbleOnce(obj)
			%disp('test');
			if obj.runOrTumble==1
				%if running
				obj.runOnce();
			else
				%disp('test');
				%if tumbling
				obj.tumbleOnce();
			end
		end

		function [x,y]=updateCoordinates(obj,pattern)
			%disp('test');
			%Check if density is high
			densityHigh=pattern.checkDensity(obj.xCoordinate,obj.yCoordinate);

			if densityHigh
				%disp('density is high');
				%then reset to tumble
				obj.resetToTumble();
				obj.highCounter=obj.highCounter+1;
			else
				%set to run
				obj.runOrTumble=1;
			end

			%run or tumble
			obj.runOrTumbleOnce();

			%correct for border
			pattern.border.correctBorder(obj);

			x=obj.xCoordinate;
			y=obj.yCoordinate;
			[m,n]=size(obj.path);
			obj.path(m+1,:)=[x y];
		end

		function plotPath(obj,fig)
			hold on;

			[m,n]=size(obj.path);

			for i=1:m-1
				plot([obj.path(i,1) obj.path(i+1,1)],[obj.path(i,2) obj.path(i+1,2)],'b');
			end
		end
	end
end
