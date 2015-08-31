%% Written bij KU Leuven iGEM team
%% Model I
close all;clear all;clc;

fid=fopen('simulation.pbs','w');
header='';
header=[header '#!/bin/bash -l\r\n'];
header=[header '#PBS -l nodes=1:ppn=20'];
header=[header '#PBS -l walltime=01:00:00'];
header=[header 'cd $hybridfast'];
header=[header 'module load matlab/R2014a'];
header=[header 'module load matlab/R2014a'];

numSimulation=2;
for simulationCounter=1:numSimulation
	poolsize=4;

	filename=['functiontest' num2str(simulationCounter)];
	framerate=25;
	scaling=20;

	%% Simulation parameters
	dt=1e-3;				%time step (h)
	tend=1e-2;			%end time of simulation
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
	XLength=1;			%Length of domain cm
	YLength=1;			%Length of domain cm
	Jx=101;				%# of subdivisions
	Jy=101;				%# of subdivisions

	%% define kernel functions and bandwidth
	bandwidth=1e-4;
	%bandwidth=10;

	%% define constants
	%Bacteria A and B
	kappa=2;
	r0=0.5e-4;		%cm
	k1=4;			%mN/cm
	k2=2;			%mN/cm
	k3=0;			%mN/cm
	%gamma=2.62e-8;	%mN*h/cm
	gamma=1e-1;	%mN*h/cm
	modulo=10;		%#
	muHighA=0.072e-3;	%high diffusion constant of bacteria A
	muLowA=0;			%low diffusion constant of bacteria A
	muHighB=2.376e-3;	%high diffusion constant of bacteria B
	muLowB=0;	%low diffusion constant of bacteria B
	VthA=0.2;
	VthB=0.1;

	%AHL and leucine
	alpha=17.9e-4;
	beta=5.4199e-4;
	kAHL=1/48;
	kleucine=1/80;
	DAHL=50e-2;
	Dleucine=26.46e-3;

end
