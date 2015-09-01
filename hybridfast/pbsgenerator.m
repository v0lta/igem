%% Written bij KU Leuven iGEM team
%% Model I
close all;clear all;clc;

outputstr1='';
outputstr2='';

outputstr1=[outputstr1 '#!/bin/bash -l\n'];
outputstr1=[outputstr1 '#PBS -l nodes=1:ppn=20\n'];
outputstr1=[outputstr1 '#PBS -l walltime=01:00:00\n'];
outputstr1=[outputstr1 'cd ~/hybridfast\n'];
outputstr1=[outputstr1 'module load matlab/R2014a\n'];

outputstr2=[outputstr2 'cd ~/hybridfast\n'];
outputstr2=[outputstr2 'module load matlab/R2014a\n'];

filename='realparameters';
poolsize=20;
numSimulation=1;

framerate=25;
scaling=20;

%% Simulation parameters
dt=1e-2;			%time step (h)
tend=24;			%end time of simulation
%tend=2;			%end time of simulation
%tend=30;			%end time of simulation
%nBacteriaA=1000;	%number of bacteria A
%nBacteriaAArray=300;		%number of bacteria A
nBacteriaAArray=1000;		%number of bacteria A
%nBacteriaA=0;		%number of bacteria A
%nBacteriaA=1;		%number of bacteria A
%nBacteriaB=1000;	%number of bacteria B
%nBacteriaBArray=300;		%number of bacteria B
nBacteriaBArray=1000;		%number of bacteria B
%nBacteriaB=0;		%number of bacteria B
%nBacteriaB=1;		%number of bacteria B
%initialpattern='gaussian';
initialpattern='uniform_random';
%initialpattern='spot';

%% define domain
XLength=1;			%Length of domain cm
YLength=1;			%Length of domain cm
Jx=101;			%# of subdivisions
Jy=101;			%# of subdivisions

%% define kernel functions and bandwidth
bandwidthArray=1e-4;
%bandwidth=10;

%% define constants
%Bacteria A and B
%kappaArray=[2,3];
kappaArray=2;
r0Array=0.5e-4;		%cm
k1Array=4;			%mN/cm
k2Array=2;			%mN/cm
k3Array=0;			%mN/cm
%gamma=2.62e-8;	%mN*h/cm
gammaArray=1e-1;	%mN*h/cm
moduloArray=10;		%#
muHighAArray=0.072e-3;	%high diffusion constant of bacteria A
muLowAArray=0;			%low diffusion constant of bacteria A
muHighBArray=2.376e-3;	%high diffusion constant of bacteria B
muLowBArray=0;	%low diffusion constant of bacteria B
VthAArray=0.2;
VthBArray=0.1;

%AHL and leucine
alphaArray=17.9e-4;
betaArray=5.4199e-4;
kAHLArray=1/48;
kleucineArray=1/80;
DAHLArray=50e-2;
DleucineArray=26.46e-3;

for nBacteriaA=nBacteriaAArray
for nBacteriaB=nBacteriaBArray
for bandwidth=bandwidthArray
for kappa=kappaArray
for r0=r0Array
for k1=k1Array
for k2=k2Array
for k3=k3Array
for gamma=gammaArray
for modulo=moduloArray
for muHighA=muHighAArray
for muLowA=muLowAArray
for muHighB=muHighBArray
for muLowB=muLowBArray
for VthA=VthAArray
for VthB=VthBArray
for alpha=alphaArray
for beta=betaArray
for kAHL=kAHLArray
for kleucine=kleucineArray
for DAHL=DAHLArray
for Dleucine=DleucineArray
	for simulationCounter=1:numSimulation
		temp=['runsimulation ',...
			filename,' ',...
			num2str(simulationCounter),' ',...
			num2str(poolsize),' ',...
			num2str(dt),' ',...
			num2str(tend),' ',...
			num2str(nBacteriaA),' ',...
			num2str(nBacteriaB),' ',...
			initialpattern,' ',...
			num2str(XLength),' ',...
			num2str(YLength),' ',...
			num2str(Jx),' ',...
			num2str(Jy),' ',...
			num2str(bandwidth),' ',...
			num2str(kappa),' ',...
			num2str(r0),' ',...
			num2str(k1),' ',...
			num2str(k2),' ',...
			num2str(k3),' ',...
			num2str(gamma),' ',...
			num2str(modulo),' ',...
			num2str(muHighA),' ',...
			num2str(muLowA),' ',...
			num2str(muHighB),' ',...
			num2str(muLowB),' ',...
			num2str(VthA),' ',...
			num2str(VthB),' ',...
			num2str(alpha),' ',...
			num2str(beta),' ',...
			num2str(kAHL),' ',...
			num2str(kleucine),' ',...
			num2str(DAHL),' ',...
			num2str(Dleucine),'\n'];
		outputstr1=[outputstr1 temp];
		outputstr2=[outputstr2 temp];
	end
end;end;end;end;end;end;end;end;end;end;end;end;end;end;end;end;end;end;end;end;end;end;

fid=fopen('simulation.pbs','w');
fprintf(fid,outputstr1);
fclose(fid);

fid=fopen('runsimulation.pbs','w');
fprintf(fid,outputstr2);
fclose(fid);
