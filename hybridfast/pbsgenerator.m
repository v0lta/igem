%% Written bij KU Leuven iGEM team
%% Model I
close all;clear all;clc;

outputstr1='';
outputstr2='';

outputstr1=[outputstr1 '#!/bin/bash -l\n'];
outputstr1=[outputstr1 '#PBS -l nodes=1:ppn=20\n'];
outputstr1=[outputstr1 '#PBS -l walltime=02:00:00\n'];
outputstr1=[outputstr1 'cd ~/hybridfast\n'];
outputstr1=[outputstr1 'module load matlab/R2014a\n'];

outputstr2=[outputstr2 'cd ~/hybridfast\n'];
outputstr2=[outputstr2 'module load matlab/R2014a\n'];

filename='real_large_test5';
poolsize=20;
numSimulation=1;

%% - realistic parameters - %%
%% Simulation parameters
dt=1e-2;			%time step (h)
%tend=10;			%end time of simulation
tend=1;			%end time of simulation
%tend=2;			%end time of simulation
%tend=30;			%end time of simulation
%nBacteriaA=1000;	%number of bacteria A
%nBacteriaAArray=300;		%number of bacteria A
nBacteriaAArray=5e3;		%number of bacteria A
%nBacteriaA=0;		%number of bacteria A
%nBacteriaA=1;		%number of bacteria A
%nBacteriaB=1000;	%number of bacteria B
%nBacteriaBArray=300;		%number of bacteria B
nBacteriaBArray=5e3;		%number of bacteria B
%nBacteriaB=0;		%number of bacteria B
%nBacteriaB=1;		%number of bacteria B
%initialpattern='gaussian';
initialpattern='uniform_random';
%initialpattern='spot';

%% define domain
XLength=10000;			%Length of domain um
YLength=10000;			%Length of domain um
Jx=2000;			%# of subdivisions
Jy=2000;			%# of subdivisions

%% define kernel functions and bandwidth
bandwidthArray=100;
%bandwidth=10;

%% define constants
%Bacteria A and B
%kappaArray=[2,3];
kappaArray=2;
r0Array=0.5;		%um
k1Array=4e-4;			%mN/um
k2Array=2e-4;			%mN/um
k3Array=5e-4;			%mN/um
%gamma=2.62e-8;	%mN*h/um
gammaArray=1e-5;	%mN*h/um
moduloArray=10;		%#
muHighAArray=0.072e5;	%high diffusion constant of bacteria A
muLowAArray=muHighAArray*1e-2;			%low diffusion constant of bacteria A
muHighBArray=2.376e5;	%high diffusion constant of bacteria B
muLowBArray=muHighBArray*1e-2;	%low diffusion constant of bacteria B
VthAArray=0.2e-3;		%nmol/cl
%VthBArray=[0.1 100];%nmol/cl
VthBArray=0.1e-3;%nmol/cl

%AHL and leucine
alphaArray=17.9e-4;	%nmol/h
betaArray=5.4199e-4;	%nmol/h
kAHLArray=1/48;		%1/h
kleucineArray=1/80;	%1/h
DAHLArray=50e5;		%um^2/h
DleucineArray=26.46e5;	%um^2/h


%%% -- toy parameters -- %%
%%% Simulation parameters
%dt=0.1;				%time step
%tend=300;			%end time of simulation
%%tend=2;			%end time of simulation
%%tend=30;			%end time of simulation
%N=tend/dt;			%time steps
%%N=200;				%time steps
%%nBacteriaA=1000;	%number of bacteria A
%nBacteriaA=300;		%number of bacteria A
%%nBacteriaA=0;		%number of bacteria A
%%nBacteriaA=1;		%number of bacteria A
%%nBacteriaB=1000;	%number of bacteria B
%nBacteriaB=300;		%number of bacteria B
%%nBacteriaB=0;		%number of bacteria B
%%nBacteriaB=1;		%number of bacteria B
%%initialpattern='gaussian';
%initialpattern='uniform_random';
%%initialpattern='spot';
%
%%% define domain
%XLength=10;			%Length of domain um
%YLength=10;			%Length of domain um
%Jx=101;				%# of subdivisions
%Jy=101;				%# of subdivisions
%
%%% define kernel functions and bandwidth
%%bandwidth=1e-4;
%bandwidthArray=.5;
%
%%% define constants
%%Bacteria A and B
%kappaArray=2;
%r0Array=0.5;		%um
%k1Array=10;			%mN/um
%k2Array=5;			%mN/um
%k3Array=20;			%mN/um
%%gammaArray=2.62e-8;	%mN*h/cm
%gammaArray=60;	%mN*h/um
%moduloArray=10;		%#
%muHighAArray=1/30;	%high diffusion constant of bacteria A um^2/h
%muLowAArray=1/300;			%low diffusion constant of bacteria A
%muHighBArray=1/30;	%high diffusion constant of bacteria B um^2/h
%muLowBArray=1/300;	%low diffusion constant of bacteria B
%VthAArray=0.2;
%VthBArray=0.1;
%
%%AHL and leucine
%%alphaArray=17.9e-4;		%nmol/h
%alphaArray=1e-3;		%nmol/h
%%betaArray=5.4199e-4;		%nmol/h
%betaArray=2e-3;		%nmol/h
%kAHLArray=5e-1;			%1/h
%kleucineArray=3e-1;		%1/h
%DAHLArray=1/30;			%cm^2/h
%DleucineArray=1/20;	%um^2/h


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
