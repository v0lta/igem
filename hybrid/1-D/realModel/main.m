%Model I
close all;clear all;clc;
filename='test.avi';
framerate=10;

N=100;			%time steps
nBacteriaA=100;	%number of bacteria A
nBacteriaB=100;	%number of bacteria B
%k=20;			%plot kth iteration
L=10;			%Length of domain
dx=0.1;			%grid spacing
domain=[0:dx:L];%domain
%domain=[0:5];%domain

%define constants
alpha=1e-3;		%bacterial production rate of AHL
beta=1e-3;		%bacterial production rate of leucine
k1=5e-3;		%degradation rate of AHL
k2=5e-3;		%degradation rate of leucine
DAHL=1/300;	%Diffusion constant of AHL
Dleucine=2/300;%Diffusion constant of leucine
%muHighA=1/30;	%high diffusion constant of bacteria A
%muLowA=1/3000;	%low diffusion constant of bacteria A
muHighA=1e-5;	%high diffusion constant of bacteria A
muLowA=1e-5;	%low diffusion constant of bacteria A
muHighB=1/30;	%high diffusion constant of bacteria B
muLowB=1/300;	%low diffusion constant of bacteria B
muA=[muLowA muHighA];
muB=[muLowB muHighB];
kappaA=2;		%chemotactic sensitivity constant bacteria A
kappaB=2;		%chemotactic sensitivity constant bacteria B
VthA=1.2;		%threshold concentration of AHL for bacteria A
VthB=1.0;		%threshold concentration of AHL for bacteria B

%initialize bacteria A
bacteriaA=[];

for i=1:nBacteriaA
	%normal distribution
	x=normrnd(L/2,1);

	%uniform random distribution
	%x=rand*L-L/2;

	%uniform distribution
	%x=L/(nBacteriaA+1)*i;

	%single peak
	%x=5;
	%x=10;

	b=bacterium(x);
	bacteriaA=[bacteriaA b];
end

bacteriaPopA=bacteriaPopulationA(bacteriaA,domain);
%bacteriaPopA=bacteriaPopulationStill(bacteriaA,domain);

%initialize bacteria B
bacteriaB=[];

for i=1:nBacteriaB
	%normal distribution
	%x=normrnd(0,1);

	%uniform random distribution
	%x=rand*L-L/2;

	%uniform distribution
	x=L/(nBacteriaB+1)*i;

	%single peak
	%x=5;
	%x=10;

	b=bacterium(x);
	bacteriaB=[bacteriaB b];
end

bacteriaPopB=bacteriaPopulationB(bacteriaB,domain);

%initialize AHL field
%uniform
n=length(domain);
concentration=zeros(1,n)+1;
%sinusoidal
%lambda=2;
%k=2*pi/lambda;
%%concentration=cos(domain)+1;
%concentration=sin(k*domain)+1; 

%boundary conditions
boundaries=[1 1];
AHLField=AHL(domain,concentration,boundaries);

%initialize leucine field
n=length(domain);
concentration=zeros(1,n)+1;
%concentration=cos(domain)+1;
%concentration=sin(domain)+1;

%boundary conditions
boundaries=[1 1];
leucineField=leucine(domain,concentration,boundaries);

%define kernel function and bandwidth
addpath ..\..\..\kernel;
%kernelfun=@uni;
%kernelfun=@tri;
%kernelfun=@(x) normpdf(x,0,1);
kernelfun=@epanechnikov;
bandwidth=0.25;

%define timestep
timestep=1;

%define scaling for plotting concentrations
scaling=50;

%initialize model
model=realModel(bacteriaPopA,bacteriaPopB,AHLField,leucineField,...
		alpha,beta,k1,k2,DAHL,Dleucine,muA,muB,kappaA,kappaB,VthA,VthB,...
		kernelfun,bandwidth,timestep,scaling);

runSimulation;
makeVideo;

beep on;
beep;
beep off;

%% Junk code

%% run for N (extra) time steps
%N=300;
%for i=1:N
%	model.update()
%end
%
%% beep on;
%% beep;
%% beep off;
%
%%plot kth iteration
%%fig=figure(1);
%%model.plot(k,fig);
%%% 
%%vidObj=VideoWriter('simulation2.avi');
%vidObj=VideoWriter(filename);
%set(vidObj,'FrameRate',framerate);
%open(vidObj);
%
%nFrames=model.getlength();
%fig=figure('Position',[100 100 1000 1000]);
%for i=1:nFrames
%	model.plot(i,fig);
%	ylim([-5 75]);
%	writeVideo(vidObj,getframe());
%	clf;
%end
%
%close(fig);
%vidObj.close();

%Double peak
%for i=1:nBacteria/2
%	x=-5;
%	b=bacterium(x);
%	bacteria=[bacteria b];
%end
%for i=(nBacteria/2+1):nBacteria
%	x=5;
%	b=bacterium(x);
%	bacteria=[bacteria b];
%end
