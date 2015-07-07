classdef stripe
	properties
		%Inner and outer radius of high cell density
		innerRadius=0;
		outerRadius
	end
	methods
		function obj=stripe(innerRadius,outerRadius)
			obj.innerRadius=innerRadius;
			obj.outerRadius=outerRadius;
		end
		
		function densityHigh=checkDensity(obj,x,y)
			r=sqrt(x.^2+y.^2);
			if r>=obj.innerRadius && r<=obj.outerRadius
				densityHigh=1;
			else
				densityHigh=0;
			end
		end

		function plot(obj,fig)
			figure(fig);

			hold on;
			%upper part
			x=-obj.outerRadius:obj.outerRadius/100:obj.outerRadius;
			y1=sqrt(obj.outerRadius^2-x.^2);
			y2=real(sqrt(obj.innerRadius^2-x.^2));
			X=[x,fliplr(x)];
			Y=[y1,fliplr(y2)];
			fill(X,Y,'r');
			
			%bottom part
			y1=-sqrt(obj.outerRadius^2-x.^2);
			y2=real(-sqrt(obj.innerRadius^2-x.^2));
			X=[x,fliplr(x)];
			Y=[y1,fliplr(y2)]; 
			fill(X,Y,'r');

			hold off;
		end
	end
end
