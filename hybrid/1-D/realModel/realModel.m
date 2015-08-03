classdef realModel<handle
	properties
		%bacterial population A and B, AHL and leucine fields
		bacteriaPopA;
		bacteriaPopB;
		AHLField;
		leucineField;

		alpha;		%production rate of AHL
		beta;		%production rate of leucine
		k1;			%degradation rate of AHL
		k2;			%degradation rate of leucine
		DAHL;		%Diffusion constant of AHL
		Dleucine;	%Diffusion constant of leucine
		muA;		%diffusion constant of bacteria A
		muB;		%diffusion constant of bacteria B
		kappaA;		%chemotactic sensitivity constant of bacteria A
		kappaB;		%chemotactic sensitivity constant of bacteria B
		VthA;		%threshold concentration of AHL of bacteria A
		VthB;		%threshold concentration of AHL of bacteria B

		%kernel function and bandwidth
		kernelfun;
		bandwidth;

		%timestep
		timestep;

		%array of density array & AHL fields
		rhoAArray;
		rhoBArray;
		AHLArray;
		leucineArray;
		coordinateAMatrix;
		coordinateBMatrix;

		%scaling for plotting concentrations
		scaling;
	end

	methods
		function obj=realModel(bacteriaPopA,bacteriaPopB,AHLField,leucineField,...
		alpha,beta,k1,k2,DAHL,Dleucine,muA,muB,kappaA,kappaB,VthA,VthB,...
		kernelfun,bandwidth,timestep,scaling)
			obj.bacteriaPopA=bacteriaPopA;
			obj.bacteriaPopB=bacteriaPopB;
			obj.AHLField=AHLField;
			obj.leucineField=leucineField;
			obj.alpha=alpha;
			obj.beta=beta;
			obj.k1=k1;
			obj.k2=k2;
			obj.DAHL=DAHL;
			obj.Dleucine=Dleucine;
			obj.muA=muA;
			obj.muB=muB;
			obj.kappaA=kappaA;
			obj.kappaB=kappaB;
			obj.VthA=VthA;
			obj.VthB=VthB;
			obj.kernelfun=kernelfun;
			obj.bandwidth=bandwidth;
			obj.timestep=timestep;
			obj.scaling=scaling;

			%record initial density functions, coordinates, AHL & leucine fields
			domain=AHLField.getdomain();

			%bacteria A
			rhoA=obj.bacteriaPopA.bacteriadensity(obj.kernelfun,obj.bandwidth);
			obj.rhoAArray=rhoA(domain);							%density
			obj.coordinateAMatrix=bacteriaPopA.coordinates();	%coordinates

			%bacteria B
			rhoB=obj.bacteriaPopB.bacteriadensity(obj.kernelfun,obj.bandwidth);
			obj.rhoBArray=rhoB(domain);							%density
			obj.coordinateBMatrix=bacteriaPopB.coordinates();	%coordinates

			%AHL and leucine fields
			obj.AHLArray=AHLField.getconcentration();			%AHL field
			obj.leucineArray=leucineField.getconcentration();	%leucine field
		end

		function update(obj)
			domain=obj.AHLField.getdomain();

			%bacteria A
			%update bacteria positions
			obj.bacteriaPopA.update(obj.AHLField,obj.muA,obj.VthA,obj.kappaA);
			%calculate bacteria density
			rhoA=obj.bacteriaPopA.bacteriadensity(obj.kernelfun,obj.bandwidth);
			%record rho
			obj.rhoAArray(end+1,:)=rhoA(domain);
			%record coordinates
			obj.coordinateAMatrix(end+1,:)=obj.bacteriaPopA.coordinates();

			%bacteria B
			%update bacteria positions
			obj.bacteriaPopB.update(obj.AHLField,obj.leucineField,obj.muB,obj.VthB,obj.kappaB);
			%calculate bacteria density
			rhoB=obj.bacteriaPopB.bacteriadensity(obj.kernelfun,obj.bandwidth);
			%record rho
			obj.rhoBArray(end+1,:)=rhoB(domain);
			%record coordinates
			obj.coordinateBMatrix(end+1,:)=obj.bacteriaPopB.coordinates();

			%update AHL field
			obj.AHLField.update(rhoA,obj.DAHL,obj.alpha,obj.k1,obj.timestep);
			%record AHL field
			obj.AHLArray(end+1,:)=obj.AHLField.getconcentration();

			%update leucine field
			obj.leucineField.update(rhoA,obj.Dleucine,obj.beta,obj.k2,obj.timestep);
			%record leucine field
			obj.leucineArray(end+1,:)=obj.leucineField.getconcentration();
		end

		function n=getlength(obj)
		[n,m]=size(obj.rhoAArray);
		end

		function plotrhos(obj,k,fig)
			figure(fig);
			hold on;

			domain=obj.AHLField.getdomain();
			rhoA=obj.rhoAArray(k,:);
			rhoB=obj.rhoBArray(k,:);
			plot(domain,rhoA);
			plot(domain,rhoB);
			%double plot
			%plot(domain+domain(end),rho);

			%pause(5);
			hold off;
		end

		function plotAHL(obj,k,fig)
			figure(fig);
			hold on;

			domain=obj.AHLField.getdomain();
			AHL=obj.AHLArray(k,:);
			%plot(domain,AHL);
			%multiply for scaling
			plot(domain,AHL*obj.scaling);
			%%double plot
			%plot(domain+domain(end),AHL*obj.scaling);

			hold off;
		end

		function plotleucine(obj,k,fig)
			figure(fig);
			hold on;

			domain=obj.leucineField.getdomain();
			leucine=obj.leucineArray(k,:);
			%plot(domain,leucine);
			%multiply for scaling
			plot(domain,leucine*obj.scaling);
			%double plot
			%plot(domain+domain(end),leucine*obj.scaling);

			hold off;
		end

		function plotbacterias(obj,k,fig)
			figure(fig);
			hold on;

			coordinateAArray=obj.coordinateAMatrix(k,:);
			coordinateBArray=obj.coordinateBMatrix(k,:);
			yA=coordinateAArray*0+1;
			yB=coordinateAArray*0+2;

			plot(coordinateAArray,yA,'*');
			plot(coordinateBArray,yB,'*');
			%double plot
			%plot(coordinateArray+1,y,'*');

			hold off;
		end

		function plotlines(obj,fig)
			figure(fig);
			hold on;

			domain=obj.AHLField.getdomain();
			n=length(domain);
			%disp('test');
			plot(domain,0*ones(n,1));
			plot(domain,obj.VthA*obj.scaling*ones(n,1));
			plot(domain,obj.VthB*obj.scaling*ones(n,1));
			%double plot
			%plot(domain+domain(end),0*ones(n,1));
			%plot(domain+domain(end),obj.Vth*obj.scaling*ones(n,1));
			
			hold off;
		end
		
		function plot(obj,k,fig)
			obj.plotrhos(k,fig);
			obj.plotAHL(k,fig);
			obj.plotleucine(k,fig);
			obj.plotlines(fig);
			obj.plotbacterias(k,fig);
			legend('Density bacteria A','Density bacteria B','AHL','Leucine');%...
			%'Threshold bacteria A','Threshold bacteria B','Bacteria A','Bacteria B');
		end
	end
end
