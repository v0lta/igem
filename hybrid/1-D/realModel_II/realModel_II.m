classdef realModel_II<handle
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
		lambda0A;	%Base turning frequency of bacteria A
		lambda0B;	%Base turning frequency of bacteria B
		speedA;		%Constant speed of bacteria A
		speedB;		%Constant speed of bacteria B
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
		function obj=realModel_II(bacteriaPopA,bacteriaPopB,AHLField,leucineField,...
		alpha,beta,k1,k2,DAHL,Dleucine,lambda0A,lambda0B,speedA,speedB,kappaA,kappaB,VthA,VthB,...
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
			obj.lambda0A=lambda0A;
			obj.lambda0B=lambda0B;
			obj.speedA=speedA;
			obj.speedB=speedB;
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
			%obj.rhoAArray=rhoA(domain);						%density
			obj.rhoAArray=obj.periodic(rhoA,domain);			%density periodic
			obj.coordinateAMatrix=bacteriaPopA.coordinates();	%coordinates

			%bacteria B
			rhoB=obj.bacteriaPopB.bacteriadensity(obj.kernelfun,obj.bandwidth);
			%obj.rhoBArray=rhoB(domain);						%density
			obj.rhoBArray=obj.periodic(rhoB,domain);			%density periodic
			obj.coordinateBMatrix=bacteriaPopB.coordinates();	%coordinates

			%AHL and leucine fields
			obj.AHLArray=AHLField.getconcentration();			%AHL field
			obj.AHLArray(end+1)=obj.AHLArray(1);				%AHL field periodic
			obj.leucineArray=leucineField.getconcentration();	%leucine field
			obj.leucineArray(end+1)=obj.leucineArray(1);		%leucine field periodic
		end

		function periodicRho=periodic(obj,rho,domain)
		periodicRho=rho(domain);
		n=length(domain);
		dx=domain(2)-domain(1);

		%Add right tail to left side
		rightx=domain(end);
		i=1;
		foo=rho(rightx+i*dx);

		while foo~=0
			periodicRho(i)=periodicRho(i)+foo;
			i=i+1;
			foo=rho(rightx+i*dx);
		end

		%Add left tail to right side
		leftx=domain(1);
		i=1;
		foo=rho(leftx-i*dx);

		while foo~=0
			periodicRho(n+1-i)=periodicRho(n+1-i)+foo;
			i=i+1;
			foo=rho(leftx-i*dx);
		end
		periodicRho(end+1)=periodicRho(1);
		end

		function update(obj)
			domain=obj.AHLField.getdomain();

			%bacteria A
			%update bacteria positions
			obj.bacteriaPopA.update(obj.AHLField,obj.lambda0A,obj.speedA,obj.kappaA,obj.VthA,obj.timestep);
			%calculate bacteria density
			rhoA=obj.bacteriaPopA.bacteriadensity(obj.kernelfun,obj.bandwidth);
			%record rho
			%obj.rhoAArray(end+1,:)=rhoA(domain);
			%record rho with periodic boundary conditions
			obj.rhoAArray(end+1,:)=obj.periodic(rhoA,domain);
			%record coordinates
			obj.coordinateAMatrix(end+1,:)=obj.bacteriaPopA.coordinates();

			%bacteria B
			%update bacteria positions
			obj.bacteriaPopB.update(obj.AHLField,obj.leucineField,...
				obj.lambda0B,obj.speedB,obj.kappaB,obj.VthB,obj.timestep);
			%calculate bacteria density
			rhoB=obj.bacteriaPopB.bacteriadensity(obj.kernelfun,obj.bandwidth);
			%record rho
			%obj.rhoBArray(end+1,:)=rhoB(domain);
			%record rho with periodic boundary conditions
			obj.rhoBArray(end+1,:)=obj.periodic(rhoB,domain);

			%record coordinates
			obj.coordinateBMatrix(end+1,:)=obj.bacteriaPopB.coordinates();

			%update AHL field
			obj.AHLField.update(rhoA,obj.DAHL,obj.alpha,obj.k1,obj.timestep);
			%record AHL field
			%obj.AHLArray(end+1,:)=obj.AHLField.getconcentration();
			%record AHL field periodic
			foo=obj.AHLField.getconcentration();
			foo(end+1)=foo(1);
			obj.AHLArray(end+1,:)=foo;

			%update leucine field
			obj.leucineField.update(rhoA,obj.Dleucine,obj.beta,obj.k2,obj.timestep);
			%record leucine field
			%obj.leucineArray(end+1,:)=obj.leucineField.getconcentration();
			%record leucine field periodic
			foo=obj.leucineField.getconcentration();
			foo(end+1)=foo(1);
			obj.leucineArray(end+1,:)=foo;
		end

		function n=getlength(obj)
		[n,m]=size(obj.rhoAArray);
		end

		function plotrhos(obj,k,fig)
			figure(fig);
			hold on;

			domain=obj.AHLField.getdomain();
			%periodic
			dx=domain(2)-domain(1);
			domain(end+1)=domain(end)+dx;
			%periodic end
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
			%periodic
			dx=domain(2)-domain(1);
			domain(end+1)=domain(end)+dx;
			%periodic end
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
			%periodic
			dx=domain(2)-domain(1);
			domain(end+1)=domain(end)+dx;
			%periodic end
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
			%periodic
			dx=domain(2)-domain(1);
			domain(end+1)=domain(end)+dx;
			%periodic end
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
			%periodic
			domain=obj.leucineField.getdomain();
			dx=domain(2)-domain(1);
			domain(end+1)=domain(end)+dx;
			xlim([domain(1) domain(end)]);
			%periodic end
			%'Threshold bacteria A','Threshold bacteria B','Bacteria A','Bacteria B');
		end
	end
end
