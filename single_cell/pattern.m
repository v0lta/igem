classdef pattern<handle
	properties
		components;
		border;
	end

	methods
		function addComp(obj,comp)
			obj.components=[obj.components comp];
		end
		
		function addBorder(obj,border)
			obj.border=border;
		end

		function densityHigh=checkDensity(obj,x,y)
			densityHigh=0;
			for comp=obj.components
				if comp.checkDensity(x,y)
					densityHigh=1;
				end
			end
		end

		function insideBorder=checkBorder(obj,x,y)
			insideBorder=obj.border.checkBorder(x,y);
		end

		function plot(obj,fig)
			figure(fig);
			hold on;
			for comp=obj.components
				comp.plot(fig);
			end
			obj.border.plot(fig);
			hold off;
		end
	end
end
