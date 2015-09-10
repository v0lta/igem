%% Written bij KU Leuven iGEM team
%% Model I
close all;clear all;

%% realistic set of parameters %%

numSimulation=1;
for simulationCounter=1:numSimulation
	poolsize=4;

	filename='10_sep_real_test';
	list=ls([filename '*_data.mat']);
	[k,~]=size(list);

	filename=[filename num2str(k+1)];

	%% Simulation parameters
	dtPDE=1;       %time step PDE (h)
    dtBact=0.1;   %time step Bact (h)
	tend=1;			%end time of simulation (h)
	%tend=2;			%end time of simulation
	%tend=30;			%end time of simulation
	N=round(tend/dtPDE);			%time steps
	%N=200;				%time steps
	%nBacteriaA=1000;	%number of bacteria A
	nBacteriaA=5e4;		%number of bacteria A
	%nBacteriaA=0;		%number of bacteria A
	%nBacteriaA=1;		%number of bacteria A
	%nBacteriaB=1000;	%number of bacteria B
	nBacteriaB=5e4;		%number of bacteria B
	%nBacteriaB=0;		%number of bacteria B
	%nBacteriaB=1;		%number of bacteria B
	%initialpattern='gaussian';
	initialpattern='uniform_random';
	%initialpattern='spot';

	%% define domain
	XLength=1000;			%Length of domain um
	YLength=1000;			%Length of domain um
	Jx=1000;				%# of subdivisions
	Jy=1000;				%# of subdivisions

	%% define kernel functions and bandwidth
	%bandwidth=1e-4;
	%bandwidth=1;
	bandwidth=5;

	%% define constants
	%Bacteria A and B
	kappa=2;
	r0=0.5;		%um
	k1=4e-4;			%mN/um
	k2=2e-4;			%mN/um
	k3=5e-4;			%mN/um
	%gamma=2.62e-8;	%mN*h/cm
	gamma=1e-5;	%mN*h/um
	modulo=10;		%#
	muHighA=0.072e5;	%high diffusion constant of bacteria A um^2/h
	muLowA=0;			%low diffusion constant of bacteria A
	muHighB=2.376e5;	%high diffusion constant of bacteria B um^2/h
	muLowB=0;	%low diffusion constant of bacteria B
	VthA=0.2;
	VthB=0.1;

	%AHL and leucine
	alpha=17.9e-4;		%nmol/h
	%alpha=17.9e-2;		%nmol/h
	beta=5.4199e-4;		%nmol/h
	%beta=5.4199e-0;		%nmol/h
	kAHL=1/48;			%1/h
	kleucine=1/80;		%1/h
	DAHL=50e5;			%um^2/h
	Dleucine=26.46e5;	%um^2/h

	runsimulation(filename,simulationCounter,poolsize,...
		dtPDE,dtBact,tend,nBacteriaA,nBacteriaB,initialpattern,...
		XLength,YLength,Jx,Jy,...
		bandwidth,...
		kappa,r0,k1,k2,k3,gamma,modulo,...
		muHighA,muLowA,muHighB,muLowB,VthA,VthB,...
		alpha,beta,kAHL,kleucine,DAHL,Dleucine);
end

% toy parameters %%
%numSimulation=1;
%for simulationCounter=1:numSimulation
%	poolsize=4;
%
%	filename='optimization_test';
%
%	list=ls([filename '*_data.mat']);
%	[k,~]=size(list);
%	filename=[filename num2str(k+1)];
%
%	%% Simulation parameters
%	dtPDE=10;       %time step PDE
%    dtBact = 1;   %time step Bact
%	tend=10;			 %end time of simulation
%	%tend=1;			%end time of simulation
%	%tend=30;			%end time of simulation
%	N=tend/dtPDE;			%time steps
%	%N=200;				%time steps
%	%nBacteriaA=1000;	%number of bacteria A
%	%nBacteriaA=300;		%number of bacteria A
%	%nBacteriaA=0;		%number of bacteria A
%	nBacteriaA=1000;		%number of bacteria A
%	%nBacteriaB=1000;	%number of bacteria B
%	%nBacteriaB=300;		%number of bacteria B
%	%nBacteriaB=0;		%number of bacteria B
%	nBacteriaB=1000;		%number of bacteria B
%	%initialpattern='gaussian';
%	initialpattern='uniform_random';
%	%initialpattern='spot';
%	%initialpattern='corner_spot';
%	%initialpattern='center_spot';
%
%	%% define domain
%	%square
%	XLength=10;			%Length of domain um
%	YLength=10;			%Length of domain um
%	Jx=101;				%# of subdivisions
%	Jy=101;				%# of subdivisions
%
%	%narrow strip
%	%XLength=100;			%Length of domain um
%	%YLength=3;			%Length of domain um
%	%Jx=1001;				%# of subdivisions
%	%Jy=30;				%# of subdivisions
%
%	%% define kernel functions and bandwidth
%	%bandwidth=1e-4;
%	bandwidth=1;
%
%	%% define constants
%	%Bacteria A and B
%	kappa=2;
%	r0=0.5;		%um
%	k1=10;			%mN/um
%	k2=5;			%mN/um
%	k3=20;			%mN/um	%original parameter
%	%k3=0;			%mN/um
%	%gamma=2.62e-8;	%mN*h/cm
%	gamma=60;	%mN*h/um
%	modulo=10;		%#
%	muHighA=1/30;	%high diffusion constant of bacteria A um^2/h
%	muLowA=1/300;			%low diffusion constant of bacteria A
%	muHighB=1/30;	%high diffusion constant of bacteria B um^2/h
%	muLowB=1/300;	%low diffusion constant of bacteria B
%	VthA=0.2;
%	VthB=0.1;
%
%	%AHL and leucine
%	%alpha=17.9e-4;		%nmol/h
%	alpha=1e-3;		%nmol/h
%	%beta=5.4199e-4;		%nmol/h
%	beta=2e-3;		%nmol/h
%	kAHL=2e-1;			%1/h
%	kleucine=1e-1;		%1/h
%	%kAHL=0;			%1/h
%	%kleucine=0;		%1/h
%	DAHL=1/30;			%cm^2/h
%	Dleucine=1/20;	%um^2/h
%
%	runsimulation(filename,simulationCounter,poolsize,...
%		dtPDE,dtBact,tend,nBacteriaA,nBacteriaB,initialpattern,...
%		XLength,YLength,Jx,Jy,...
%		bandwidth,...
%		kappa,r0,k1,k2,k3,gamma,modulo,...
%		muHighA,muLowA,muHighB,muLowB,VthA,VthB,...
%		alpha,beta,kAHL,kleucine,DAHL,Dleucine);
%
%	beep on;beep;beep off;
%end

