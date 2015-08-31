function runsimulation(filename,simulationCounter,poolsize,...
	dt,tend,nBacteriaA,nBacteriaB,initialpattern,...
	XLength,YLength,Jx,Jy,...
	bandwidth,...
	kappa,r0,k1,k2,k3,gamma,modulo,...
	muHighA,muLowA,muHighB,muLowB,VthA,VthB,...
	alpha,beta,kAHL,kleucine,DAHL,Dleucine)
	%run simulation with given parameters and output data in filename.csv
	
	if (ischar(simulationCounter)), simulationCounter=str2num(simulationCounter),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;
	if (ischar(ab)), ab=str2num(ab),end;

	p=gcp('nocreate');
	if isempty(p)
		parpool(poolsize);
	elseif p.NumWorkers~=poolsize
		delete(p);
		parpool(poolsize);
	end

	%% process simulation parameters
	N=tend/dt;			%time steps

	%% process domain
	domain.x=linspace(0,XLength,Jx);
	domain.y=linspace(0,YLength,Jy);

	%% process kernel functions and bandwidth
	kernelfun=@epanechnikov2DNorm;

	paramAB.kernelfun=kernelfun;
	paramAB.bandwidth=bandwidth;

	%% define constants
	%Bacteria A and B

	%% -- fast -- %%
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
	bacteriaA=[];

	for i=1:nBacteriaA
		switch initialpattern
		case 'spot'
			x=XLength/2;
			y=YLength/2;
		case 'gaussian'
			%gaussian
			x=normrnd(1/2*XLength,1);
			y=normrnd(1/2*YLength,1);
			%x=normrnd(9/10*XLength,1);
			%y=normrnd(9/10*YLength,1);
		case 'uniform_random'
			%uniform random
			x=rand*XLength;
			y=rand*YLength;
		otherwise
			warning('Unknown distribution, defaulting to uniform random');
			%uniform random
			x=rand*XLength;
			y=rand*YLength;
		end

		if x>XLength
			x=XLength;
		elseif x<0
			x=0;
		end

		if y>YLength
			y=YLength;
		elseif y<0
			y=0;
		end

		bacteriaA=[bacteriaA;x y];
	end

	%% Initialize bacteria B
	bacteriaB=[];

	for i=1:nBacteriaB
		switch initialpattern
		case 'spot'
			x=XLength/2;
			y=YLength/2;
		case 'gaussian'
			%gaussian
			x=normrnd(1/2*XLength,1);
			y=normrnd(1/2*YLength,1);
			%x=normrnd(9/10*XLength,1);
			%y=normrnd(9/10*YLength,1);
		case 'uniform_random'
			%uniform random
			x=rand*XLength;
			y=rand*YLength;
		otherwise
			warning('Unknown distribution, defaulting to uniform random');
			%uniform random
			x=rand*XLength;
			y=rand*YLength;
		end

		if x>XLength
			x=XLength;
		elseif x<0
			x=0;
		end

		if y>YLength
			y=YLength;
		elseif y<0
			y=0;
		end

		bacteriaB=[bacteriaB;x y];
	end

	%% initialize AHL field
	m=length(domain.y);
	n=length(domain.x);
	[X,Y]=meshgrid(domain.x,domain.y);

	%uniform
	AHLconcentration=zeros(m,n)+1e-5;

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
	[X,Y]=meshgrid(domain.x,domain.y);

	%uniform
	leucineconcentration=zeros(m,n)+1e-5;

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

	%model=model1(paramA,paramB,paramAHL,paramleucine,...
	model=modelfast1(paramAB,paramAHL,paramleucine,...
		bacteriaA,bacteriaB,...
		AHLconcentration,AHLboundaries,leucineconcentration,leucineboundaries,...
		domain);

	%% run simulation
	disp(['Running simulation ' num2str(simulationCounter)]);
	t1=tic;
	for i=1:N
		%if mod(i,1)==0
		if mod(i,N/10)==0
			disp('----- %%%%% -----');
			%now
			disp(datestr(now));
			disp(['Current iteration: ' num2str(i) '/' num2str(N)]);
		end
		model.update(dt);
		%if mod(i,1)==0
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

	beep on;beep;beep off;

	disp('Simulation finished');
	toc(t1);

	%export data
	AHLArray=model.AHLArray;
	leucineArray=model.leucineArray;
	rhoAArray=model.rhoAArray;
	rhoBArray=model.rhoBArray;
	coordinateAMatrix=model.coordinateAMatrix;
	coordinateBMatrix=model.coordinateBMatrix;

	%write arguments
	fid=fopen([filename '_arguments.csv'],'w');
	header=['dt,tend,nBacteriaA,nBacteriaB,initialpattern,',...
		'XLength,YLength,Jx,Jy,',...
		'bandwidth,',...
		'kappa,r0,k1,k2,k3,gamma,modulo,',...
		'muHighA,muLowA,muHighB,muLowB,VthA,VthB,',...
		'alpha,beta,kAHL,kleucine,DAHL,Dleucine\r\n'];
	fprintf(fid,header);

	floatingSpec='%.5e,';
	formatSpec=[repmat(floatingSpec,1,2) '%d,%d,%s,',...
		repmat(floatingSpec,1,2),'%d,%d,',...
		floatingSpec,...
		repmat(floatingSpec,1,6),'%d,',...
		repmat(floatingSpec,1,6),...
		repmat(floatingSpec,1,6)];
	formatSpec(end)='';
	formatSpec=[formatSpec '\r\n'];
	fprintf(fid,formatSpec,dt,tend,nBacteriaA,nBacteriaB,initialpattern,...
		XLength,YLength,Jx,Jy,...
		bandwidth,...
		kappa,r0,k1,k2,k3,gamma,modulo,...
		muHighA,muLowA,muHighB,muLowB,VthA,VthB,...
		alpha,beta,kAHL,kleucine,DAHL,Dleucine);
	fclose(fid);
		
	%write AHLArray
	fid=fopen([filename '_AHLArray.csv'],'w');

	outputstr='';
	formatSpec=repmat(floatingSpec,1,Jx);
	formatSpec(end)='';
	formatSpec=[formatSpec ';'];
	parfor i=1:N+1
		temp=sprintf(['t=' floatingSpec(1:end-1) ';\r\n'],(i-1)*dt);
		temp=[temp 'AHLArray=['];
		temp=[temp sprintf(formatSpec,AHLArray(:,:,i)') '];\r\n'];
		outputstr=[outputstr temp];
	end
	fprintf(fid,outputstr);
	fclose(fid);

	%write leucineArray
	fid=fopen([filename '_leucineArray.csv'],'w');

	outputstr='';
	formatSpec=repmat(floatingSpec,1,Jx);
	formatSpec(end)='';
	formatSpec=[formatSpec ';'];
	parfor i=1:N+1
		temp=sprintf(['t=' floatingSpec(1:end-1) ';\r\n'],(i-1)*dt);
		temp=[temp 'leucineArray=['];
		temp=[temp sprintf(formatSpec,leucineArray(:,:,i)') '];\r\n'];
		outputstr=[outputstr temp];
	end
	fprintf(fid,outputstr);
	fclose(fid);

	%write rhoAArray
	fid=fopen([filename '_rhoAArray.csv'],'w');

	outputstr='';
	formatSpec=repmat(floatingSpec,1,Jx);
	formatSpec(end)='';
	formatSpec=[formatSpec ';'];
	parfor i=1:N+1
		temp=sprintf(['t=' floatingSpec(1:end-1) ';\r\n'],(i-1)*dt);
		temp=[temp 'rhoAArray=['];
		temp=[temp sprintf(formatSpec,rhoAArray(:,:,i)') '];\r\n'];
		outputstr=[outputstr temp];
	end
	fprintf(fid,outputstr);
	fclose(fid);

	%write rhoBArray
	fid=fopen([filename '_rhoBArray.csv'],'w');

	outputstr='';
	formatSpec=repmat(floatingSpec,1,Jx);
	formatSpec(end)='';
	formatSpec=[formatSpec ';'];
	parfor i=1:N+1
		temp=sprintf(['t=' floatingSpec(1:end-1) ';\r\n'],(i-1)*dt);
		temp=[temp 'rhoBArray=['];
		temp=[temp sprintf(formatSpec,rhoBArray(:,:,i)') '];\r\n'];
		outputstr=[outputstr temp];
	end
	fprintf(fid,outputstr);
	fclose(fid);

	%write coordinateAMatrix
	fid=fopen([filename '_coordinateAMatrix.csv'],'w');

	outputstr='';
	formatSpec=repmat(floatingSpec,1,nBacteriaA);
	formatSpec(end)='';
	formatSpec=[formatSpec '];\r\n'];
	parfor i=1:N+1
		temp=sprintf(['t=' floatingSpec(1:end-1) ';\r\n'],(i-1)*dt);
		temp=[temp 'bactX=[' sprintf(formatSpec,coordinateAMatrix(1,:,i)')];
		temp=[temp 'bactY=[' sprintf(formatSpec,coordinateAMatrix(2,:,i)')];
		outputstr=[outputstr temp];
	end
	fprintf(fid,outputstr);
	fclose(fid);

	%write coordinateBMatrix
	fid=fopen([filename '_coordinateBMatrix.csv'],'w');

	outputstr='';
	formatSpec=repmat(floatingSpec,1,nBacteriaB);
	formatSpec(end)='';
	formatSpec=[formatSpec '];\r\n'];
	parfor i=1:N+1
		temp=sprintf(['t=' floatingSpec(1:end-1) ';\r\n'],(i-1)*dt);
		temp=[temp 'bactX=[' sprintf(formatSpec,coordinateBMatrix(1,:,i)')];
		temp=[temp 'bactY=[' sprintf(formatSpec,coordinateBMatrix(2,:,i)')];
		outputstr=[outputstr temp];
	end
	fprintf(fid,outputstr);
	fclose(fid);
end
