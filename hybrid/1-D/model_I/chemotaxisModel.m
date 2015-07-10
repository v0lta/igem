classdef chemotaxisModel<handle
	properties
		%bacterial population and nutrient field
		bacteriaPop;
		nutrientField;

		mu;			%diffusion constant
		kappa;		%chemotactic sensitivity constant
		d;			%consumption rate

		%kernel function and bandwidth
		kernelfun;
		bandwidth;

		%array of density array & nutrient fields
		rhoArray;
		nutrientArray;
		coordinateMatrix;
	end

	methods
		function obj=chemotaxisModel(bacteriaPop,nutrientField,mu,kappa,d,kernelfun,bandwidth)
			obj.bacteriaPop=bacteriaPop;
			obj.nutrientField=nutrientField;
			obj.mu=mu;
			obj.kappa=kappa;
			obj.d=d;
			obj.kernelfun=kernelfun;
			obj.bandwidth=bandwidth;

			%record initial density function, coordinates & nutrient fields
			domain=nutrientField.getdomain();
			rho=obj.bacteriaPop.bacteriadensity(obj.kernelfun,obj.bandwidth);
			obj.rhoArray=rho(domain);	%density
			obj.coordinateMatrix=bacteriaPop.coordinates();	%coordinates
			obj.nutrientArray=nutrientField.getconcentration();	%nutrient field
		end

		function update(obj)
			%calculate bacteria density
			rho=obj.bacteriaPop.bacteriadensity(obj.kernelfun,obj.bandwidth);
			%update nutrient field
			obj.nutrientField.update(rho,obj.d);
			%record nutrient field
			obj.nutrientArray(end+1,:)=obj.nutrientField.getconcentration();

			%update bacteria positions
			obj.bacteriaPop.update(obj.nutrientField,obj.mu,obj.kappa);
			%record rho
			domain=obj.nutrientField.getdomain();
			rho=obj.bacteriaPop.bacteriadensity(obj.kernelfun,obj.bandwidth);
			obj.rhoArray(end+1,:)=rho(domain);
			%record coordinates
			obj.coordinateMatrix(end+1,:)=obj.bacteriaPop.coordinates();
		end

		function plotrho(obj,k,fig)
			figure(fig);
			hold on;

			domain=obj.nutrientField.getdomain();
			rho=obj.rhoArray(k,:);
			plot(domain,rho);

			%pause(5);
			hold off;
		end

		function plotnutrient(obj,k,fig)
			figure(fig);
			hold on;

			domain=obj.nutrientField.getdomain();
			nutrient=obj.nutrientArray(k,:);
			%plot(domain,nutrient);
			plot(domain,nutrient*50);

			hold off;
		end

		function plotbacteria(obj,k,fig)
			figure(fig);
			hold on;

			coordinateArray=obj.coordinateMatrix(k,:);
			y=coordinateArray*0;

			plot(coordinateArray,y,'*');

			hold off;
		end
		

		function plot(obj,k,fig)
			obj.plotrho(k,fig);
			obj.plotnutrient(k,fig);
			obj.plotbacteria(k,fig);
		end
	end
end
