%% Written bij KU Leuven iGEM team
%% Model I
close all;clear all;clc;

poolsize=4;
p=gcp('nocreate');
if isempty(p)
	parpool(poolsize);
elseif p.NumWorkers~=poolsize
	delete(p);
	parpool(poolsize);
end

numSimulation=1;
for simulationCounter=1:numSimulation
%numSimulation=0;
%while 1
%numSimulation=numSimulation+1;
%filename=['zero_flux_',...
%	'zero_initial_concentration_',...
%	'high_degradation_',...
%	'uniform_A_',...
%	'uniform_B_',...
%	'rectangular_domain_',...
%	'small_simulation'];
%filename='cell_interaction_large_simulation_test2';
filename=['fasttest' num2str(simulationCounter)];
framerate=25;
scaling=20;

%% Simulation parameters
dt=0.1;				%time step
%dt=1;				%time step
tend=50;			%end time of simulation
N=tend/dt;			%time steps
%N=200;				%time steps
nBacteriaA=50;	%number of bacteria A
%nBacteriaA=300;		%number of bacteria A
%nBacteriaA=1;		%number of bacteria A
nBacteriaB=50;	%number of bacteria B
%nBacteriaB=300;		%number of bacteria B
%nBacteriaB=1;		%number of bacteria B
%initialpattern='gaussian';
initialpattern='uniform_random';
%initialpattern='spot';

%% define domain
XLength=10;			%Length of domain
YLength=10;			%Length of domain
Jx=101;				%# of subdivisions
Jy=101;				%# of subdivisions
domain.x=linspace(0,XLength,Jx);
domain.y=linspace(0,YLength,Jy);
%domain=[domainx',domainy'];

%% define kernel functions and bandwidth
kernelfun=@epanechnikov2DNorm;
bandwidth=.5;
%bandwidth=10;

paramAB.kernelfun=kernelfun;
paramAB.bandwidth=bandwidth;
paramA.kernelfun=kernelfun;
paramB.kernelfun=kernelfun;
paramA.bandwidth=bandwidth;
paramB.bandwidth=bandwidth;

%% define constants
%Bacteria A and B
kappa=2;
r0=0.5;
k1=10;
k2=5;
k3=20;
gamma=60;
modulo=10;

%% -- fast -- %%
paramAB.r0=r0;			%cell radius
paramAB.rcut=r0*1.25;		%cut off radius for attraction 
paramAB.rsearch=r0*1.25*2;	%search radius
paramAB.k1=k1;			%spring constant
paramAB.k2=k2;			%spring constant
paramAB.k3=k3;			%spring constant
paramAB.gamma=gamma;			%cell sensitivity
paramAB.modulo=modulo;		%number of iterations between refreshes

muHighA=1/30;	%high diffusion constant of bacteria A
muLowA=1/300;	%low diffusion constant of bacteria A
%muLowA=1/30;	%low diffusion constant of bacteria A
paramAB.muA.low=muLowA;
paramAB.muA.high=muHighA;
%paramAB.VthA=1.2;		%threshold concentration of AHL for bacteria A
paramAB.VthA=0.2;		%threshold concentration of AHL for bacteria A
paramAB.kappaA=kappa;		%chemotactic sensitivity constant of bacteria A

muHighB=1/30;	%high diffusion constant of bacteria B
muLowB=1/300;	%low diffusion constant of bacteria B
%muLowB=1/30;	%low diffusion constant of bacteria B
paramAB.muB.low=muLowB;
paramAB.muB.high=muHighB;
%paramAB.VthB=1.0;		%threshold concentration of AHL for bacteria B
paramAB.VthB=0.1;		%threshold concentration of AHL for bacteria B
paramAB.kappaB=kappa;		%chemotactic sensitivity constant of bacteria B

%% -- normal -- %%
muHighA=1/30;	%high diffusion constant of bacteria A
muLowA=1/300;	%low diffusion constant of bacteria A
%muLowA=1/30;	%low diffusion constant of bacteria A
paramA.muA.low=muLowA;
paramA.muA.high=muHighA;
paramA.kappaA=kappa;		%chemotactic sensitivity constant of bacteria A
%paramA.VthA=1.2;		%threshold concentration of AHL for bacteria A
paramA.VthA=0.2;		%threshold concentration of AHL for bacteria A
paramA.r0=r0;			%cell radius
paramA.rcut=r0*1.25;		%cut off radius for attraction 
paramA.k1=k1;			%spring constant
paramA.k2=k2;			%spring constant
paramA.k3=k3;			%spring constant
paramA.gamma=gamma;			%cell sensitivity
paramA.modulo=modulo;		%number of iterations between refreshes

muHighB=1/30;	%high diffusion constant of bacteria B
muLowB=1/300;	%low diffusion constant of bacteria B
%muLowB=1/30;	%low diffusion constant of bacteria B
paramB.muB.low=muLowB;
paramB.muB.high=muHighB;
paramB.kappaB=kappa;		%chemotactic sensitivity constant of bacteria B
%paramB.VthB=1.0;		%threshold concentration of AHL for bacteria B
paramB.VthB=0.1;		%threshold concentration of AHL for bacteria B
paramB.r0=r0;			%cell radius
paramB.k1=k1;			%spring constant
paramB.k2=k2;			%spring constant
paramB.gamma=gamma;			%cell sensitivity
paramB.modulo=modulo;		%number of iterations between refreshes


%AHL and leucine
paramAHL.alpha=1e-3;		%production rate of AHL
%alpha=0;		%production rate of AHL
%alpha=-3e-3;	%production rate of AHL
paramleucine.beta=2e-3;		%production rate of leucine
%k1=5e-3;		%degradation rate of AHL
%k1=0;			%degradation rate of AHL
paramAHL.k1=5e-1;		%degradation rate of AHL
paramleucine.k2=3e-1;		%degradation rate of leucine
%DAHL=1/300;	%Diffusion constant of AHL
paramAHL.DAHL=1/30;		%Diffusion constant of AHL
paramleucine.Dleucine=1/20;	%Diffusion constant of leucine

%% Initialize bacteria A
bacteriaA=[];

%b=bacteriumA(XLength/2-r0*1.1,YLength/2);
%b=bacteriumA(XLength,YLength);
%bacteriaA=[bacteriaA b];
%b=bacteriumA(XLength/2+r0*1.1,YLength/2);
%bacteriaA=[bacteriaA b];

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

	%b=bacteriumA(x,y);
	%bacteriaA=[bacteriaA b];
	bacteriaA=[bacteriaA;x y];
end

%% Initialize bacteria B
bacteriaB=[];

%b=bacteriumB(XLength/2+r0*0.3,YLength/2);
%bacteriaB=[bacteriaB b];
%b=bacteriumB(XLength/2-r0*0.3,YLength/2);
%bacteriaB=[bacteriaB b];

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

	%b=bacteriumB(x,y);
	%bacteriaB=[bacteriaB b];
	bacteriaB=[bacteriaB;x y];
end

%% initialize AHL field
m=length(domain.y);
n=length(domain.x);
[X,Y]=meshgrid(domain.x,domain.y);

%uniform
%concentration=zeros(m,n)+1;
AHLconcentration=zeros(m,n)+1e-5;

%block
%a=floor(J/3);
%idy=a:J-a;
%idx=a:J-a;
%AHLconcentration(idy,idx)=ones(length(idy),length(idx));

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
%leucineconcentration=zeros(m,n)+1;
leucineconcentration=zeros(m,n)+1e-5;

%block
%a=floor(J/3);
%idy=a:J-a;
%idx=a:J-a;
%leucineconcentration(idy,idx)=ones(length(idy),length(idx));

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

%% Analytic parameters
paramAnal.scaling=scaling;			%scaling for plotting concentrations
paramAnal.framerate=framerate;	%framerate

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
		disp(['Current iteration: ' num2str(i)]);

		%elapsed
		tElapsed=toc(t1);
		timeString1=displaytime(tElapsed);
		disp([timeString1 ' have elapsed']);

		%ETA
		tPerIter=tElapsed/i;
		tETA=tPerIter*(N-i);
		timeString2=displaytime(tETA);
		disp([timeString2 ' left']);
	end
	model.update(dt);
end

beep on;beep;beep off;

disp('Simulation finished');
toc(t1);

analObject=analyzer(paramAnal,model);

%% preview cell cell
%analObject.interactionpreview();

%% preview
%disp('Preview of simulation');
%analObject.preview();

%% save workspace and make videos
disp('Saving workspace and videos');
t2=tic;
save(filename);
analObject.makevideo(filename);
disp('Workspace and videos saved');
toc(t2);

%% end!
disp('End!');
end
disp('End of all simulations!');
