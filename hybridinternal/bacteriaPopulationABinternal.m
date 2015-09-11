classdef bacteriaPopulationABInternal < bacteriaPopulationAB
	properties	
		%Internal concentrations of bacteria A & B, n x N matrix, every column represent one bacterium
		concA;
		concB;

		%Governing equations
		ODEA;
		ODEB;

		%Solver
		solver;
	end

	methods
		function obj=bacteriaPopulationABInternal(paramAB,bacteriaA,bacteriaB,domain,domainGrid,...
			ODEA,ODEB,initA,initB,solver)
		%superclass definition
		obj@bacteriaPopulationAB(paramAB,bacteriaA,bacteriaB,domain,domainGrid);

		%record initial concentrations
		[~,nA]=size(initA);
		[~,nB]=size(initB);
		if nA==1
			obj.concA=repmat(initA,1,obj.numA);
		else
			obj.concA=initA;
		end
		if nB==1
			obj.concB=repmat(initB,1,obj.numB);
		else
			obj.concB=initB;
		end

		%ode system for A and B
		obj.ODEA=ODEA;
		obj.ODEB=ODEB;

		%ODE solver
		obj.solver=solver;
		end

		function addbacteria(obj,bacteriaA,bacteriaB,concA,concB)
		%add list of bacteria and their concentrations

		obj@bacteriaPopulationAB.addbacteria(bacteriaA,bacteriaB)

		obj.concA=[obj.concA,concA];
		obj.concB=[obj.concB,concB];
		end

		function addbacteriumA(obj,x,y,conc)
		%add a bacterium of type A

		obj@bacteriaPopulationAB.addbacteriumA(x,y);

		obj.concA=[obj.concA,conc];
		end

		function addbacteriumB(obj,x,y,conc)
		%add a bacterium of type B

		obj@bacteriaPopulationAB.addbacteriumB(x,y);

		obj.concB=[obj.concB,conc];
		end

		function concA=concentrationA(obj)
		%return internal concentrations of bacteria A

		concA=obj.concA;
		end

		function concB=concentrationB(obj)
		%return internal concentrations of bacteria B

		concB=obj.concB;
		end

		function updatefast(obj,AHLField,leucineField,dt)
		%update based on AHL field, leuceine field and timestep

		obj@bacteriaPopulationAB.updatefast(AHLField,leucineField,dt);

		numA=obj.numA;
		solver=obj.solver;
		ODEA=obj.ODEA;
		ODEB=obj.ODEB;
		concA=obj.concA;
		concB=obj.concB;

		parfor i=1:numA
			currentConc=concA(:,i);

			[T,newConc]=solver(ODEA,[0 dt],currentConc)
			newConc=(newConc(end,:))';

			concA(:,i)=newConc;
		end

		obj.concA=concA;

		parfor i=1:numB
			currentConc=concB(:,i);

			[T,newConc]=solver(ODEB,[0 dt],currentConc)
			newConc=(newConc(end,:))';

			concB(:,i)=newConc;
		end

		obj.concB=concB;
		end
	end
end
