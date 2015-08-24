close all;clear all;clc;
filename=['gaussian_with_interaction_large_simulation'];
paramAnal.framerate=10;
paramAnal.tPause=0.1;

%% Simulation parameters
dt=1;				%time step
tend=20;			%end time of simulation
N=tend/dt;			%time steps
%N=200;				%time steps
nBacteriaA=300;	%number of bacteria A

%% define domain
XLength=10;			%Length of domain
YLength=10;			%Length of domain
Jx=101;				%# of subdivisions
Jy=101;				%# of subdivisions
domain.x=linspace(0,XLength,Jx);
domain.y=linspace(0,YLength,Jy);
domainGrid=meshgrid(domain.x,domain.y);

%% define kernel functions and bandwidth
kernelfun=@epanechnikov2DNorm;
bandwidth=.5;

paramA.kernelfun=kernelfun;
paramA.bandwidth=bandwidth;

%% define constants
%Bacteria A and B
%muA=1/30;
muHighA=1/300;	%high diffusion constant of bacteria A
paramA.muA=muHighA;
paramA.kappaA=2;		%chemotactic sensitivity constant of bacteria A
paramA.VthA=1.2;		%threshold concentration of AHL for bacteria A
paramA.r0=2;			%cell radius
paramA.k=2;			%spring constant
paramA.gamma=2;			%cell sensitivity
paramA.modulo=5;		%number of iterations between refreshes

%% Initialize bacteria A
bacteriaA=[];

for i=1:nBacteriaA
	%x=XLength/2;
	%y=YLength/2;

	%gaussian
	x=normrnd(1/2*XLength,1);
	y=normrnd(1/2*YLength,1);
	%x=normrnd(9/10*XLength,1);
	%y=normrnd(9/10*YLength,1);

	%uniform random
	%x=rand*XLength;
	%y=rand*YLength;

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

	b=bacterium(x,y);
	bacteriaA=[bacteriaA b];
end

model=cellmodel(paramA,bacteriaA,domain);

%% run simulation
disp('Running simulation');
tic;
for i=1:N
	%if mod(i,N/10)==0
	if mod(i,1)==0
		disp(['iteration ' num2str(i)]);
	end
	model.update(dt);
end

beep on;beep;beep off;

disp('Simulation finished');
toc;

analObject=analyzer(paramAnal,model);

%% preview
%disp('Preview of simulation');
%analObject.preview();

%% save workspace and make videos
analObject.makevideo(filename);

%% end!
disp('End!');

%With function
%161.85s
%159.34s
%152.58s
%149.87s

%Embedded
%143.03s
%140.56s
%143.10s
%

%parallel
%231.17s

%serial
%229.07s
