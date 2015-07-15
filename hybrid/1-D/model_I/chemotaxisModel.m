classdef chemotaxisModel<handle
	properties
		%bacterial population and attractant field
		bacteriaPop;
		attractantField;
		nutrientField;

		mu;			%diffusion constant
		Vth;		%threshold concentration of attractant
		kappa;		%chemotactic sensitivity constant
		eta;		%bacteria production rate
		Da;			%Diffusion constant of attractant
		alpha;		%degradation rate of attractant
		Ds;			%Diffusion constant of attractant
		beta;		%degradation rate of attractant

		%kernel function and bandwidth
		kernelfun;
		bandwidth;

		%timestep
		timestep;

		%array of density array & attractant fields
		rhoArray;
		attractantArray;
		nutrientArray;
		coordinateMatrix;

		%scaling for plotting concentrations
		scaling;
	end

	methods
		function obj=chemotaxisModel(bacteriaPop,attractantField,nutrientField,mu,Vth,kappa,eta,Da,alpha,...
		Ds,beta,kernelfun,bandwidth,timestep)
			obj.bacteriaPop=bacteriaPop;
			obj.attractantField=attractantField;
			obj.nutrientField=nutrientField;
			obj.mu=mu;
			obj.Vth=Vth;
			obj.kappa=kappa;
			obj.eta=eta;
			obj.Da=Da;
			obj.alpha=alpha;
			obj.Ds=Ds;
			obj.beta=beta;
			obj.kernelfun=kernelfun;
			obj.bandwidth=bandwidth;
			obj.timestep=timestep;
			obj.scaling=50;

			%record initial density function, coordinates, attractant & nutrient fields
			domain=attractantField.getdomain();
			rho=obj.bacteriaPop.bacteriadensity(obj.kernelfun,obj.bandwidth);
			obj.rhoArray=rho(domain);	%density
			obj.coordinateMatrix=bacteriaPop.coordinates();	%coordinates
			obj.attractantArray=attractantField.getconcentration();	%attractant field
			obj.nutrientArray=nutrientField.getconcentration();	%nutrient field
		end

		function update(obj)
			%update bacteria positions
			obj.bacteriaPop.update(obj.attractantField,obj.mu,obj.Vth,obj.kappa);
			%calculate bacteria density
			rho=obj.bacteriaPop.bacteriadensity(obj.kernelfun,obj.bandwidth);

			%update attractant field
			obj.attractantField.update(rho,obj.Da,obj.eta,obj.alpha,obj.timestep);
			%record attractant field
			obj.attractantArray(end+1,:)=obj.attractantField.getconcentration();

			%update nutrient field
			obj.nutrientField.update(rho,obj.Ds,obj.beta,obj.timestep);
			%record nutrient field
			obj.nutrientArray(end+1,:)=obj.nutrientField.getconcentration();

			%record rho
			domain=obj.attractantField.getdomain();
			rho=obj.bacteriaPop.bacteriadensity(obj.kernelfun,obj.bandwidth);
			obj.rhoArray(end+1,:)=rho(domain);
			%record coordinates
			obj.coordinateMatrix(end+1,:)=obj.bacteriaPop.coordinates();
		end

		function n=getlength(obj)
		[n,m]=size(obj.rhoArray);
		end

		function plotrho(obj,k,fig)
			figure(fig);
			hold on;

			domain=obj.attractantField.getdomain();
			rho=obj.rhoArray(k,:);
			plot(domain,rho);

			%pause(5);
			hold off;
		end

		function plotattractant(obj,k,fig)
			figure(fig);
			hold on;

			domain=obj.attractantField.getdomain();
			attractant=obj.attractantArray(k,:);
			%plot(domain,attractant);
			%multiply for scaling
			plot(domain,attractant*obj.scaling);

			hold off;
		end

		function plotnutrient(obj,k,fig)
			figure(fig);
			hold on;

			domain=obj.nutrientField.getdomain();
			nutrient=obj.nutrientArray(k,:);
			%plot(domain,nutrient);
			%multiply for scaling
			plot(domain,nutrient*obj.scaling);

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

		function plotlines(obj,fig)
			figure(fig);
			hold on;

			domain=obj.nutrientField.getdomain();
			n=length(domain);
			%disp('test');
			plot(domain,0*ones(n,1));
			plot(domain,obj.Vth*obj.scaling*ones(n,1));
			
			hold off;
		end
		
		function plot(obj,k,fig)
			obj.plotrho(k,fig);
			obj.plotattractant(k,fig);
			obj.plotnutrient(k,fig);
			obj.plotbacteria(k,fig);
			obj.plotlines(fig);
		end
	end
end
