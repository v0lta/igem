classdef circleBorder
	properties
		radius;
	end

	methods
		function obj=circleBorder(radius)
			obj.radius=radius;
		end

		function insideBorder=checkBorder(obj,x,y)
			r=sqrt(x.^2+y.^2);
			insideBorder=r<obj.radius
		end

		function correctBorder(obj,bacterium)
			x=bacterium.xCoordinate;
			y=bacterium.yCoordinate;

			r=sqrt(x^2+y^2);

			%if outside border
			if r>=obj.radius
				%rescale to border
				bacterium.xCoordinate=x*obj.radius/r;
				bacterium.yCoordinate=y*obj.radius/r;
			
				%randomize angle away from border
				theta=acos(1-bacterium.STEP_SIZE^2/(2*obj.radius^2));
				bacterium.angle=bacterium.angle+pi/2+theta/2+rand*(pi-theta);

				%reset counter
				bacterium.runCounter=bacterium.MAX_RUN_COUNTER;
			end
		end


		function plot(obj,fig)
			figure(fig);
			hold on;

			x=-obj.radius:obj.radius/100:obj.radius;
			y1=sqrt(obj.radius^2-x.^2);
			y2=-sqrt(obj.radius^2-x.^2);

			plot(x,y1,'k');
			plot(x,y2,'k');
		end
	end
end
