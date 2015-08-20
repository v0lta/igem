close all;clear all;clc;
filename=['gaussian_no_interaction_large_radius'];
framerate=10;

%% Simulation parameters
dt=1;				%time step
tend=50;			%end time of simulation
N=tend/dt;			%time steps
%N=200;				%time steps
nBacteriaA=10;	%number of bacteria A

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
bandwidth=.8;

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
paramA.k=1;			%cell radius
paramA.f=1;			%cell radius
paramA.modulo=5;		%cell radius

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

%% Analytic parameters
paramAnal.framerate=framerate;	%framerate

model=cellmodel(paramA,bacteriaA,domain);

%% run simulation
disp('Running simulation');
for i=1:N
	model.update(dt);
end

beep on;beep;beep off;

disp('Simulation finished');

analObject=analyzer(paramAnal,model);

%% preview
%disp('Preview of simulation');
%analObject.preview();

%% save workspace and make videos
analObject.makevideo(filename);

%% end!
disp('End!');
