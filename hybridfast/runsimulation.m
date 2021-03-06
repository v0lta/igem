function runsimulation(filename,simulationCounter,poolsize,...
	dtPDE,dtBact,tend,nBacteriaA,nBacteriaB,initialpattern,...
	XLength,YLength,Jx,Jy,...
	bandwidth,...
	kappa,r0,k1,k2,k3,gamma,modulo,...
	muHighA,muLowA,muHighB,muLowB,VthA,VthB,...
	alpha,beta,kAHL,kleucine,DAHL,Dleucine,minConcentration)
	%run simulation with given parameters and output data in filename.m
	
	%convert str to numeric values
	if (ischar(simulationCounter)), simulationCounter=str2num(simulationCounter);end;
	if (ischar(poolsize)), poolsize=str2num(poolsize);end;
	if (ischar(dtPDE)), dtPDE=str2num(dtPDE);end;
    if (ischar(dtBact)), dtBact=str2num(dtBact);end;
	if (ischar(tend)), tend=str2num(tend);end;
	if (ischar(nBacteriaA)), nBacteriaA=str2num(nBacteriaA);end;
	if (ischar(nBacteriaB)), nBacteriaB=str2num(nBacteriaB);end;
	if (ischar(XLength)), XLength=str2num(XLength);end;
	if (ischar(YLength)), YLength=str2num(YLength);end;
	if (ischar(Jx)), Jx=str2num(Jx);end;
	if (ischar(Jy)), Jy=str2num(Jy);end;
	if (ischar(bandwidth)), bandwidth=str2num(bandwidth);end;
	if (ischar(kappa)), kappa=str2num(kappa);end;
	if (ischar(r0)), r0=str2num(r0);end;
	if (ischar(k1)), k1=str2num(k1);end;
	if (ischar(k2)), k2=str2num(k2);end;
	if (ischar(k3)), k3=str2num(k3);end;
	if (ischar(gamma)), gamma=str2num(gamma);end;
	if (ischar(modulo)), modulo=str2num(modulo);end;
	if (ischar(muHighA)), muHighA=str2num(muHighA);end;
	if (ischar(muLowA)), muLowA=str2num(muLowA);end;
	if (ischar(muHighB)), muHighB=str2num(muHighB);end;
	if (ischar(muLowB)), muLowB=str2num(muLowB);end;
	if (ischar(VthA)), VthA=str2num(VthA);end;
	if (ischar(VthB)), VthB=str2num(VthB);end;
	if (ischar(alpha)), alpha=str2num(alpha);end;
	if (ischar(beta)), beta=str2num(beta);end;
	if (ischar(kAHL)), kAHL=str2num(kAHL);end;
	if (ischar(kleucine)), kleucine=str2num(kleucine);end;
	if (ischar(DAHL)), DAHL=str2num(DAHL);end;
	if (ischar(Dleucine)), Dleucine=str2num(Dleucine);end;
	if (ischar(minConcentration)), minConcentration=str2num(minConcentration);end;

	disp('Saving parameters in matlab file');
	t2=tic;
	save([filename '_data.mat'],'-v7.3',...
		'dtPDE',...
		'dtBact',...
		'tend',...
		'nBacteriaA',...
		'nBacteriaB',...
		'initialpattern',...
		'XLength',...
		'YLength',...
		'Jx',...
		'Jy',...
		'bandwidth',...
		'kappa',...
		'r0',...
		'k1',...
		'k2',...
		'k3',...
		'gamma',...
		'modulo',...
		'muHighA',...
		'muLowA',...
		'muHighB',...
		'muLowB',...
		'VthA',...
		'VthB',...
		'alpha',...
		'beta',...
		'kAHL',...
		'kleucine',...
		'DAHL',...
		'Dleucine',...
		'minConcentration');
	disp('Parameters saved');
	toc(t2);

	%supercomputer: uncomment the lines below
	%myCluster=parcluster('local');
	%delete(myCluster.Jobs);

	p=gcp('nocreate');
	if isempty(p)
		parpool(poolsize);
	elseif p.NumWorkers~=poolsize
		delete(p);
		parpool(poolsize);
	end

	%% process simulation parameters
	N=tend/dtPDE;			%time steps

	%% process domain
	domain.x=linspace(0,XLength,Jx);
	domain.y=linspace(0,YLength,Jy);

	%% process kernel functions and bandwidth
	kernelfun=@epanechnikov2DNorm;

	paramAB.kernelfun=kernelfun;
	paramAB.bandwidth=bandwidth;

	%% define constants
	%Bacteria A and B
	paramAB.r0=r0;			%cell radius
	paramAB.rcut=r0*1.25;		%cut off radius for attraction 
	paramAB.rsearch=r0*1.25*2;	%search radius
	paramAB.k1=k1;			%spring constant
	paramAB.k2=k2;			%spring constant
	paramAB.k3=k3;			%spring constant
	paramAB.gamma=gamma;			%cell sensitivity
	paramAB.modulo=modulo;		%number of iterations between refreshes

	paramAB.muA.low=muLowA;
	paramAB.muA.high=muHighA;
	paramAB.VthA=VthA;		%threshold concentration of AHL for bacteria A
	paramAB.kappaA=kappa;		%chemotactic sensitivity constant of bacteria A

	paramAB.muB.low=muLowB;
	paramAB.muB.high=muHighB;
	paramAB.VthB=VthB;		%threshold concentration of AHL for bacteria B
	paramAB.kappaB=kappa;		%chemotactic sensitivity constant of bacteria B

	%AHL and leucine
	paramAHL.alpha=alpha;		%production rate of AHL
	paramleucine.beta=beta;		%production rate of leucine
	paramAHL.k1=kAHL;		%degradation rate of AHL
	paramleucine.k2=kleucine;		%degradation rate of leucine
	paramAHL.DAHL=DAHL;		%Diffusion constant of AHL
	paramleucine.Dleucine=Dleucine;	%Diffusion constant of leucine

	%% Initialize bacteria A
	%periodic
	dx=domain.x(2)-domain.x(1);
	dy=domain.y(2)-domain.y(1);
	disp('Initializing bacteria');
	t1=tic;
	bacteriaA=zeros(nBacteriaA,2);

	switch initialpattern
	case 'center_spot'
		%bacteriaA=[ones(nBacteriaA,1)*XLength*1/2,ones(nBacteriaA,1)*YLength/2];
		%periodic
		bacteriaA=[ones(nBacteriaA,1)*(XLength+dx)*1/2,ones(nBacteriaA,1)*(YLength+dy)/2];
	case 'corner_spot'
		%bacteriaA=[ones(nBacteriaA,1)*XLength,ones(nBacteriaA,1)*YLength];
		bacteriaA=[ones(nBacteriaA,1)*XLength-3,ones(nBacteriaA,1)*YLength-3];
	case 'gaussian'
		bacteriaA=[normrnd(1/2*XLength,XLength/10,nBacteriaA,1),normrnd(1/2*YLength,YLength/10,nBacteriaA,1)];
		%periodic
		bacteriaA(:,1)=mod(bacteriaA(:,1),XLength+dx);
		bacteriaA(:,2)=mod(bacteriaA(:,2),YLength+dy);
		%parfor i=1:nBacteriaA
		%	XY=bacteriaA(i,:);
		%	x=XY(1);
		%	y=XY(2);

		%	if x>XLength
		%		x=XLength;
		%	elseif x<0
		%		x=0;
		%	else
		%		x=x;
		%	end

		%	if y>YLength
		%		y=YLength;
		%	elseif y<0
		%		y=0;
		%	else
		%		y=y
		%	end

		%	bacteriaA(i,:)=[x,y];
		%end
	case 'uniform_random'
		%bacteriaA=[rand(nBacteriaA,1)*XLength,rand(nBacteriaA,1)*YLength];
		%periodic
		bacteriaA=[rand(nBacteriaA,1)*(XLength+dx),rand(nBacteriaA,1)*(YLength+dy)];
	otherwise
		warning('Unknown distribution, defaulting to uniform random');
		%bacteriaA=[rand(nBacteriaA,1)*XLength,rand(nBacteriaA,1)*YLength];
		%periodic
		bacteriaA=[rand(nBacteriaA,1)*(XLength+dx),rand(nBacteriaA,1)*(YLength+dy)];
	end

%	bacteriaA=[8.5,8.5;...
%		1.5,1.5];

	%% Initialize bacteria B
	bacteriaB=zeros(nBacteriaB,2);

	switch initialpattern
	case 'center_spot'
		%bacteriaB=[ones(nBacteriaB,1)*XLength*1/2,ones(nBacteriaB,1)*YLength/2];
		%periodic
		bacteriaB=[ones(nBacteriaB,1)*(XLength+dx)*1/2,ones(nBacteriaB,1)*(YLength+dy)/2];
	case 'corner_spot'
		%bacteriaB=[ones(nBacteriaB,1)*XLength,ones(nBacteriaB,1)*YLength];
		bacteriaB=[ones(nBacteriaB,1)*XLength-3,ones(nBacteriaB,1)*YLength-3];
	case 'gaussian'
		bacteriaB=[normrnd(1/2*XLength,XLength/10,nBacteriaB,1),normrnd(1/2*YLength,YLength/10,nBacteriaB,1)];
		%periodic
		bacteriaB(:,1)=mod(bacteriaB(:,1),XLength+dx);
		bacteriaB(:,2)=mod(bacteriaB(:,2),YLength+dy);
		%parfor i=1:nBacteriaB
		%	XY=bacteriaB(i,:);
		%	x=XY(1);
		%	y=XY(2);

		%	if x>XLength
		%		x=XLength;
		%	elseif x<0
		%		x=0;
		%	else
		%		x=x;
		%	end

		%	if y>YLength
		%		y=YLength;
		%	elseif y<0
		%		y=0;
		%	else
		%		y=y
		%	end

		%	bacteriaB(i,:)=[x,y];
		%end
	case 'uniform_random'
		%bacteriaB=[rand(nBacteriaB,1)*XLength,rand(nBacteriaB,1)*YLength];
		%periodic
		bacteriaB=[rand(nBacteriaB,1)*(XLength+dx),rand(nBacteriaB,1)*(YLength+dy)];
	otherwise
		warning('Unknown distribution, defaulting to uniform random');
		%bacteriaB=[rand(nBacteriaB,1)*XLength,rand(nBacteriaB,1)*YLength];
		%periodic
		bacteriaB=[rand(nBacteriaB,1)*(XLength+dx),rand(nBacteriaB,1)*(YLength+dy)];
	end
	%bacteriaB=[XLength,YLength];
	%bacteriaB=[XLength/2,YLength/2];

	disp('Bacteria initialization finished');
	toc(t1);

	%bacteriaA
	%bacteriaB

	%% initialize AHL field
	disp('Initializing fields');
	t1=tic;

	m=length(domain.y);
	n=length(domain.x);
	%[X,Y]=meshgrid(domain.x,domain.y);

	%uniform
	%AHLconcentration=zeros(m,n)+1e-5;
	AHLconcentration=ones(m,n)*minConcentration;
    %AHLconcentration(end-10:end,end-10:end) = 10;

	%boundary conditions
	west=AHLconcentration(:,1);
	east=AHLconcentration(:,end);
	south=AHLconcentration(1,:);
	north=AHLconcentration(end,:);
	if Jx<Jy
		foo=0;
		bar=Jy-Jx;
	else
		foo=Jx-Jy;
		bar=0;
	end
	AHLboundaries=[[west;zeros(foo,1)],[east;zeros(foo,1)],[south';zeros(bar,1)],[north';zeros(bar,1)]];

	%% initialize leucine field
	m=length(domain.y);
	n=length(domain.x);
	%[X,Y]=meshgrid(domain.x,domain.y);

	%uniform
	leucineconcentration=ones(m,n)*minConcentration;
    %leucineconcentration(end-10:end,end-10:end) = 10;

	%boundary conditions
	west=leucineconcentration(:,1);
	east=leucineconcentration(:,end);
	south=leucineconcentration(1,:);
	north=leucineconcentration(end,:);
	if Jx<Jy
		foo=0;
		bar=Jy-Jx;
	else
		foo=Jx-Jy;
		bar=0;
	end
	leucineboundaries=[[west;zeros(foo,1)],[east;zeros(foo,1)],[south';zeros(bar,1)],[north';zeros(bar,1)]];

	disp('Fields initialization finished');
	toc(t1);

	%% Set up model
	disp('Setting up model');
	t1=tic;
	%model=model1(paramA,paramB,paramAHL,paramleucine,...
	model=modelfast1(filename,paramAB,paramAHL,paramleucine,...
		bacteriaA,bacteriaB,...
		AHLconcentration,AHLboundaries,leucineconcentration,leucineboundaries,...
		domain,N);
	disp('Model setup finished');
	toc(t1);

	%% run simulation
	disp(['Running simulation ' num2str(simulationCounter)]);
	t1=tic;
	for i=1:N
		%if mod(i,1)==0
		%if mod(i,N/10)==0 || i==10
		if mod(i,N/10)==0
			disp('----- %%%%% -----');
			%now
			disp(datestr(now));
			disp(['Current iteration: ' num2str(i) '/' num2str(N)]);
		end
		model.update(dtPDE,dtBact);
		%if mod(i,1)==0
		%if mod(i,N/10)==0 || i==10
		if mod(i,N/10)==0
			%elapsed
			tElapsed=toc(t1);
			timeString1=displaytime(tElapsed);
			disp([timeString1 ' have elapsed']);

			%ETA
			tPerIter=tElapsed/i;
			tETA=tPerIter*(N-i);
			timeString2=displaytime(tETA);
			disp([timeString2 ' left']);

			%total
			tTotal=tElapsed/i*N;
			timeString3=displaytime(tTotal);
			disp([timeString3 ' total time']);
		end
	end

	disp('Simulation finished');
	toc(t1);

	%% export data
	disp('Saving data in matlab file');
	t2=tic;

	mFile=matfile([filename '_data.mat'],'Writable',true);
	mFile.AHLArray=model.AHLArray;
	mFile.leucineArray=model.leucineArray;
	mFile.rhoAArray=model.rhoAArray;
	mFile.rhoBArray=model.rhoBArray;
	mFile.coordinateAMatrix=model.coordinateAMatrix;
	mFile.coordinateBMatrix=model.coordinateBMatrix;

	disp('Data saved');
	toc(t2);
end
