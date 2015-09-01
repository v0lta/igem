function runsimulation(filename,simulationCounter,poolsize,...
	dt,tend,nBacteriaA,nBacteriaB,initialpattern,...
	XLength,YLength,Jx,Jy,...
	bandwidth,...
	kappa,r0,k1,k2,k3,gamma,modulo,...
	muHighA,muLowA,muHighB,muLowB,VthA,VthB,...
	alpha,beta,kAHL,kleucine,DAHL,Dleucine)
	%run simulation with given parameters and output data in filename.m
	
	%convert str to numeric values
	if (ischar(simulationCounter)), simulationCounter=str2num(simulationCounter);,end;
	if (ischar(poolsize)), poolsize=str2num(poolsize);,end;
	if (ischar(dt)), dt=str2num(dt);,end;
	if (ischar(tend)), tend=str2num(tend);,end;
	if (ischar(nBacteriaA)), nBacteriaA=str2num(nBacteriaA);,end;
	if (ischar(nBacteriaB)), nBacteriaB=str2num(nBacteriaB);,end;
	if (ischar(XLength)), XLength=str2num(XLength);,end;
	if (ischar(YLength)), YLength=str2num(YLength);,end;
	if (ischar(Jx)), Jx=str2num(Jx);,end;
	if (ischar(Jy)), Jy=str2num(Jy);,end;
	if (ischar(bandwidth)), bandwidth=str2num(bandwidth);,end;
	if (ischar(kappa)), kappa=str2num(kappa);,end;
	if (ischar(r0)), r0=str2num(r0);,end;
	if (ischar(k1)), k1=str2num(k1);,end;
	if (ischar(k2)), k2=str2num(k2);,end;
	if (ischar(k3)), k3=str2num(k3);,end;
	if (ischar(gamma)), gamma=str2num(gamma);,end;
	if (ischar(modulo)), modulo=str2num(modulo);,end;
	if (ischar(muHighA)), muHighA=str2num(muHighA);,end;
	if (ischar(muLowA)), muLowA=str2num(muLowA);,end;
	if (ischar(muHighB)), muHighB=str2num(muHighB);,end;
	if (ischar(muLowB)), muLowB=str2num(muLowB);,end;
	if (ischar(VthA)), VthA=str2num(VthA);,end;
	if (ischar(VthB)), VthB=str2num(VthB);,end;
	if (ischar(alpha)), alpha=str2num(alpha);,end;
	if (ischar(beta)), beta=str2num(beta);,end;
	if (ischar(kAHL)), kAHL=str2num(kAHL);,end;
	if (ischar(kleucine)), kleucine=str2num(kleucine);,end;
	if (ischar(DAHL)), DAHL=str2num(DAHL);,end;
	if (ischar(Dleucine)), Dleucine=str2num(Dleucine);,end;

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
	bacteriaA=zeros(nBacteriaA,2);

	switch initialpattern
	case 'spot'
		bacteriaA=[ones(nBacteriaA,1)*XLength/2,ones(nBacteriaA,1)*YLength/2];
	case 'gaussian'
		bacteriaA=[normrnd(1/2*XLength,1,nBacteriaA,1),normrnd(1/2*YLength,1,nBacteriaA,1)];
		parfor i=1:nBacteriaA
			XY=bacteriaA(i,:);
			x=XY(1);
			y=XY(2);

			if x>XLength
				x=XLength;
			elseif x<0
				x=0;
			else
				x=x;
			end

			if y>YLength
				y=YLength;
			elseif y<0
				y=0;
			else
				y=y
			end

			bacteriaA(i,:)=[x,y];
		end
	case 'uniform_random'
		bacteriaA=[rand(nBacteriaA,1)*XLength,rand(nBacteriaA,1)*YLength];
	otherwise
		warning('Unknown distribution, defaulting to uniform random');
		bacteriaA=[rand(nBacteriaA,1)*XLength,rand(nBacteriaA,1)*YLength];
	end

	%% Initialize bacteria B
	bacteriaB=zeros(nBacteriaB,2);

	switch initialpattern
	case 'spot'
		bacteriaB=[ones(nBacteriaB,1)*XLength/2,ones(nBacteriaB,1)*YLength/2];
	case 'gaussian'
		bacteriaB=[normrnd(1/2*XLength,1,nBacteriaB,1),normrnd(1/2*YLength,1,nBacteriaB,1)];
		parfor i=1:nBacteriaB
			XY=bacteriaB(i,:);
			x=XY(1);
			y=XY(2);

			if x>XLength
				x=XLength;
			elseif x<0
				x=0;
			else
				x=x;
			end

			if y>YLength
				y=YLength;
			elseif y<0
				y=0;
			else
				y=y
			end

			bacteriaB(i,:)=[x,y];
		end
	case 'uniform_random'
		bacteriaB=[rand(nBacteriaB,1)*XLength,rand(nBacteriaB,1)*YLength];
	otherwise
		warning('Unknown distribution, defaulting to uniform random');
		bacteriaB=[rand(nBacteriaB,1)*XLength,rand(nBacteriaB,1)*YLength];
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
		if mod(i,N/10)==0 || i==10
			disp('----- %%%%% -----');
			%now
			disp(datestr(now));
			disp(['Current iteration: ' num2str(i) '/' num2str(N)]);
		end
		model.update(dt);
		%if mod(i,1)==0
		if mod(i,N/10)==0 || i==10
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

	%beep on;beep;beep off;

	disp('Simulation finished');
	toc(t1);

	%export data
	AHLArray=model.AHLArray;
	leucineArray=model.leucineArray;
	%disp(max(max(max(leucineArray))));
	rhoAArray=model.rhoAArray;
	rhoBArray=model.rhoBArray;
	coordinateAMatrix=model.coordinateAMatrix;
	coordinateBMatrix=model.coordinateBMatrix;

	%write arguments
	%fid=fopen([filename '_arguments.m'],'w');
	%header=['dt,tend,nBacteriaA,nBacteriaB,initialpattern,',...
	%	'XLength,YLength,Jx,Jy,',...
	%	'bandwidth,',...
	%	'kappa,r0,k1,k2,k3,gamma,modulo,',...
	%	'muHighA,muLowA,muHighB,muLowB,VthA,VthB,',...
	%	'alpha,beta,kAHL,kleucine,DAHL,Dleucine\r\n'];
	%fprintf(fid,header);

	%formatSpec=[repmat(floatingSpec,1,2) '%d,%d,%s,',...
	%	repmat(floatingSpec,1,2),'%d,%d,',...
	%	floatingSpec,...
	%	repmat(floatingSpec,1,6),'%d,',...
	%	repmat(floatingSpec,1,6),...
	%	repmat(floatingSpec,1,6)];
	%formatSpec(end)='';
	%formatSpec=[formatSpec '\r\n'];

	fid=fopen([filename '_arguments.m'],'w');
	floatingSpec='%.5e,';
	fSpec='%.5e;';
	formatSpec=['dt=',fSpec,'\r\n',...
		'tend=',fSpec,'\r\n',...
		'nBacteriaA=%d;\r\n',...
		'nBacteriaB=%d;\r\n',...
		'initialpattern=''%s'';\r\n',...
		'XLength=',fSpec,'\r\n',...
		'YLength=',fSpec,'\r\n',...
		'Jx=',fSpec,'\r\n',...
		'Jy=',fSpec,'\r\n',...
		'bandwidth=',fSpec,'\r\n',...
		'kappa=',fSpec,'\r\n',...
		'r0=',fSpec,'\r\n',...
		'k1=',fSpec,'\r\n',...
		'k2=',fSpec,'\r\n',...
		'k3=',fSpec,'\r\n',...
		'gamma=',fSpec,'\r\n',...
		'modulo=%d;\r\n',...
		'muHighA=',fSpec,'\r\n',...
		'muLowA=',fSpec,'\r\n',...
		'muHighB=',fSpec,'\r\n',...
		'muLowB=',fSpec,'\r\n',...
		'VthA=',fSpec,'\r\n',...
		'VthB=',fSpec,'\r\n',...
		'alpha=',fSpec,'\r\n',...
		'beta=',fSpec,'\r\n',...
		'kAHL=',fSpec,'\r\n',...
		'kleucine=',fSpec,'\r\n',...
		'DAHL=',fSpec,'\r\n',...
		'Dleucine=',fSpec,'\r\n'];

	fprintf(fid,formatSpec,dt,tend,nBacteriaA,nBacteriaB,initialpattern,...
		XLength,YLength,Jx,Jy,...
		bandwidth,...
		kappa,r0,k1,k2,k3,gamma,modulo,...
		muHighA,muLowA,muHighB,muLowB,VthA,VthB,...
		alpha,beta,kAHL,kleucine,DAHL,Dleucine);
	fclose(fid);
		
	%write AHLArray
	fid=fopen([filename '_AHLArray.m'],'w');

	outputstr=sprintf('AHLArray=zeros(%d,%d,%d);',Jy,Jx,N+1);
	formatSpec=repmat(floatingSpec,1,Jx);
	formatSpec(end)='';
	formatSpec=[formatSpec ';'];
	parfor i=1:N+1
		temp=sprintf(['t=' floatingSpec(1:end-1) ';\r\n'],(i-1)*dt);
		temp=[temp 'AHLArray(:,:,' num2str(i) ')=['];
		temp=[temp sprintf(formatSpec,AHLArray(:,:,i)') '];\r\n'];
		outputstr=[outputstr temp];
	end
	fprintf(fid,outputstr);
	fclose(fid);

	%write leucineArray
	fid=fopen([filename '_leucineArray.m'],'w');

	outputstr=sprintf('leucineArray=zeros(%d,%d,%d);',Jy,Jx,N+1);
	formatSpec=repmat(floatingSpec,1,Jx);
	formatSpec(end)='';
	formatSpec=[formatSpec ';'];
	%parfor i=1:N+1
	for i=1:N+1
		temp=sprintf(['t=' floatingSpec(1:end-1) ';\r\n'],(i-1)*dt);
		temp=[temp 'leucineArray(:,:,' num2str(i) ')=['];
		temp=[temp sprintf(formatSpec,leucineArray(:,:,i)') '];\r\n'];
		outputstr=[outputstr temp];
		clf;
	end
	fprintf(fid,outputstr);
	fclose(fid);

	%write rhoAArray
	fid=fopen([filename '_rhoAArray.m'],'w');

	outputstr=sprintf('rhoAArray=zeros(%d,%d,%d);',Jy,Jx,N+1);
	formatSpec=repmat(floatingSpec,1,Jx);
	formatSpec(end)='';
	formatSpec=[formatSpec ';'];
	parfor i=1:N+1
		temp=sprintf(['t=' floatingSpec(1:end-1) ';\r\n'],(i-1)*dt);
		temp=[temp 'rhoAArray(:,:,' num2str(i) ')=['];
		temp=[temp sprintf(formatSpec,rhoAArray(:,:,i)') '];\r\n'];
		outputstr=[outputstr temp];
	end
	fprintf(fid,outputstr);
	fclose(fid);

	%write rhoBArray
	fid=fopen([filename '_rhoBArray.m'],'w');

	outputstr=sprintf('rhoBArray=zeros(%d,%d,%d);',Jy,Jx,N+1);
	formatSpec=repmat(floatingSpec,1,Jx);
	formatSpec(end)='';
	formatSpec=[formatSpec ';'];
	parfor i=1:N+1
		temp=sprintf(['t=' floatingSpec(1:end-1) ';\r\n'],(i-1)*dt);
		temp=[temp 'rhoBArray(:,:,' num2str(i) ')=['];
		temp=[temp sprintf(formatSpec,rhoBArray(:,:,i)') '];\r\n'];
		outputstr=[outputstr temp];
	end
	fprintf(fid,outputstr);
	fclose(fid);

	%write coordinateAMatrix
	fid=fopen([filename '_coordinateAMatrix.m'],'w');

	outputstr=sprintf('coordinateAMatrix=zeros(%d,%d,%d);\r\n',nBacteriaA,2,N+1);
	formatSpec=repmat(floatingSpec,1,nBacteriaA);
	formatSpec(end)='';
	formatSpec=[formatSpec ']'';\r\n'];
	parfor i=1:N+1
		temp=sprintf(['t=' floatingSpec(1:end-1) ';\r\n'],(i-1)*dt);
		temp=[temp 'bactX=[' sprintf(formatSpec,coordinateAMatrix(1,:,i)')];
		temp=[temp 'bactY=[' sprintf(formatSpec,coordinateAMatrix(2,:,i)')];
		temp=[temp 'coordinateAMatrix(:,:,' num2str(i) ')=[bactX,bactY];\r\n'];
		outputstr=[outputstr temp];
	end
	fprintf(fid,outputstr);
	fclose(fid);

	%write coordinateBMatrix
	fid=fopen([filename '_coordinateBMatrix.m'],'w');

	outputstr=sprintf('coordinateBMatrix=zeros(%d,%d,%d);\r\n',nBacteriaB,2,N+1);
	formatSpec=repmat(floatingSpec,1,nBacteriaB);
	formatSpec(end)='';
	formatSpec=[formatSpec ']'';\r\n'];
	parfor i=1:N+1
		temp=sprintf(['t=' floatingSpec(1:end-1) ';\r\n'],(i-1)*dt);
		temp=[temp 'bactX=[' sprintf(formatSpec,coordinateBMatrix(1,:,i)')];
		temp=[temp 'bactY=[' sprintf(formatSpec,coordinateBMatrix(2,:,i)')];
		temp=[temp 'coordinateBMatrix(:,:,' num2str(i) ')=[bactX,bactY];\r\n'];
		outputstr=[outputstr temp];
	end
	fprintf(fid,outputstr);
	fclose(fid);
end
