%% Written bij KU Leuven iGEM team
%% Model I
close all;clear all;clc;

%numSimulation=1;
%for simulationCounter=1:numSimulation
%	poolsize=4;
%
%	filename=['functiontest' num2str(simulationCounter)];
%	framerate=25;
%	scaling=20;
%
%	%% Simulation parameters
%	dt=1e-2;				%time step (h)
%	tend=1e-1;			%end time of simulation
%	%tend=2;			%end time of simulation
%	%tend=30;			%end time of simulation
%	N=tend/dt;			%time steps
%	%N=200;				%time steps
%	%nBacteriaA=1000;	%number of bacteria A
%	nBacteriaA=10;		%number of bacteria A
%	%nBacteriaA=0;		%number of bacteria A
%	%nBacteriaA=1;		%number of bacteria A
%	%nBacteriaB=1000;	%number of bacteria B
%	nBacteriaB=10;		%number of bacteria B
%	%nBacteriaB=0;		%number of bacteria B
%	%nBacteriaB=1;		%number of bacteria B
%	%initialpattern='gaussian';
%	initialpattern='uniform_random';
%	%initialpattern='spot';
%
%	%% define domain
%	XLength=10;			%Length of domain um
%	YLength=10;			%Length of domain um
%	Jx=101;				%# of subdivisions
%	Jy=101;				%# of subdivisions
%
%	%% define kernel functions and bandwidth
%	%bandwidth=1e-4;
%	bandwidth=1;
%
%	%% define constants
%	%Bacteria A and B
%	kappa=2;
%	r0=0.5;		%um
%	k1=4e-4;			%mN/um
%	k2=2e-4;			%mN/um
%	k3=5e-4;			%mN/um
%	%gamma=2.62e-8;	%mN*h/cm
%	gamma=1e-5;	%mN*h/um
%	modulo=10;		%#
%	muHighA=0.072e5;	%high diffusion constant of bacteria A um^2/h
%	muLowA=0;			%low diffusion constant of bacteria A
%	muHighB=2.376e5;	%high diffusion constant of bacteria B um^2/h
%	muLowB=0;	%low diffusion constant of bacteria B
%	VthA=0.2;
%	VthB=0.1;
%
%	%AHL and leucine
%	%alpha=17.9e-4;		%nmol/h
%	alpha=17.9e-2;		%nmol/h
%	%beta=5.4199e-4;		%nmol/h
%	beta=5.4199e-0;		%nmol/h
%	kAHL=1/48;			%1/h
%	kleucine=1/80;		%1/h
%	DAHL=50e6;			%cm^2/h
%	Dleucine=26.46e5;	%um^2/h
%
%	runsimulation(filename,simulationCounter,poolsize,...
%		dt,tend,nBacteriaA,nBacteriaB,initialpattern,...
%		XLength,YLength,Jx,Jy,...
%		bandwidth,...
%		kappa,r0,k1,k2,k3,gamma,modulo,...
%		muHighA,muLowA,muHighB,muLowB,VthA,VthB,...
%		alpha,beta,kAHL,kleucine,DAHL,Dleucine);
%end

numSimulation=1;
for simulationCounter=1:numSimulation
	poolsize=4;

	filename=['toytest' num2str(simulationCounter)];
	framerate=25;
	scaling=20;

	%% Simulation parameters
	dt=0.1;				%time step
	tend=10;			%end time of simulation
	%tend=2;			%end time of simulation
	%tend=30;			%end time of simulation
	N=tend/dt;			%time steps
	%N=200;				%time steps
	%nBacteriaA=1000;	%number of bacteria A
	nBacteriaA=300;		%number of bacteria A
	%nBacteriaA=0;		%number of bacteria A
	%nBacteriaA=1;		%number of bacteria A
	%nBacteriaB=1000;	%number of bacteria B
	nBacteriaB=300;		%number of bacteria B
	%nBacteriaB=0;		%number of bacteria B
	%nBacteriaB=1;		%number of bacteria B
	%initialpattern='gaussian';
	initialpattern='uniform_random';
	%initialpattern='spot';

	%% define domain
	XLength=10;			%Length of domain um
	YLength=10;			%Length of domain um
	Jx=101;				%# of subdivisions
	Jy=101;				%# of subdivisions

	%% define kernel functions and bandwidth
	%bandwidth=1e-4;
	bandwidth=.5;

	%% define constants
	%Bacteria A and B
	kappa=2;
	r0=0.5;		%um
	k1=10;			%mN/um
	k2=5;			%mN/um
	k3=20;			%mN/um
	%gamma=2.62e-8;	%mN*h/cm
	gamma=60;	%mN*h/um
	modulo=10;		%#
	muHighA=1/30;	%high diffusion constant of bacteria A um^2/h
	muLowA=1/300;			%low diffusion constant of bacteria A
	muHighB=1/30;	%high diffusion constant of bacteria B um^2/h
	muLowB=1/300;	%low diffusion constant of bacteria B
	VthA=0.2;
	VthB=0.1;

	%AHL and leucine
	%alpha=17.9e-4;		%nmol/h
	alpha=1e-3;		%nmol/h
	%beta=5.4199e-4;		%nmol/h
	beta=2e-3;		%nmol/h
	kAHL=5e-1;			%1/h
	kleucine=3e-1;		%1/h
	DAHL=1/30;			%cm^2/h
	Dleucine=1/20;	%um^2/h

	runsimulation(filename,simulationCounter,poolsize,...
		dt,tend,nBacteriaA,nBacteriaB,initialpattern,...
		XLength,YLength,Jx,Jy,...
		bandwidth,...
		kappa,r0,k1,k2,k3,gamma,modulo,...
		muHighA,muLowA,muHighB,muLowB,VthA,VthB,...
		alpha,beta,kAHL,kleucine,DAHL,Dleucine);
end
