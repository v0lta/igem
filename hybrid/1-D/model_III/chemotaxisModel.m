%Model II
classdef chemotaxisModel<handle
	properties
		%bacterial population and nutrient field
		bacteriaPop;
		nutrientField;

		%mu;		%equal to speed^2/(2*lambda)
		kappa;		%chemotactic sensitivity constant
		d;			%consumption rate
		lambda0;	%base turning frequency
		speed;		%constant speed
		Ds;			%Diffusion constant of nutrients

		%kernel function and bandwidth
		kernelfun;
		bandwidth;

		%time stepsize
		timestep;

		%array of density array & nutrient fields
		rhoArray;
		nutrientArray;
		coordinateMatrix;
	end

	methods
		function obj=chemotaxisModel(bacteriaPop,nutrientField,kappa,d,lambda0,speed,Ds,kernelfun,bandwidth,timestep)
			obj.bacteriaPop=bacteriaPop;
			obj.nutrientField=nutrientField;
			obj.kappa=kappa;
			obj.d=d;
			obj.lambda0=lambda0;
			obj.speed=speed;
			obj.Ds=Ds;
			obj.kernelfun=kernelfun;
			obj.bandwidth=bandwidth;
			obj.timestep=timestep;

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
			obj.nutrientField.update(rho,obj.d,obj.Ds,obj.timestep);
			%record nutrient field
			obj.nutrientArray(end+1,:)=obj.nutrientField.getconcentration();

			%update bacteria positions
			obj.bacteriaPop.update(obj.nutrientField,obj.lambda0,obj.kappa,obj.speed,obj.timestep);
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
