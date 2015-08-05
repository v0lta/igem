%Model II
close all;clear all;clc;
filename='turning_frequency_modulated_diffusion_abrupt_spot.avi';
framerate=10;

N=300;			%time steps
nBacteriaA=100;	%number of bacteria A
nBacteriaB=100;	%number of bacteria B
%k=20;			%plot kth iteration
L=15;			%Length of domain
dx=0.1;			%grid spacing
domain=[0:dx:L];%domain
%domain=[0:5];%domain

%define constants
alpha=1e-3;		%bacterial production rate of AHL
beta=1e-3;		%bacterial production rate of leucine
%alpha=0;		%bacterial production rate of AHL
%beta=0;		%bacterial production rate of leucine
k1=5e-3;		%degradation rate of AHL
k2=5e-3;		%degradation rate of leucine
DAHL=1/300;	%Diffusion constant of AHL
Dleucine=2/300;%Diffusion constant of leucine

%variable speed
%speedHighA=1e-2;	%high constant speed of bacteria A
%speedLowA=.5e-2;		%low constant speed of bacteria A
%speedHighB=1e-2;	%high constant speed of bacteria B
%speedLowB=0.7e-2;		%low constant speed of bacteria B
%%speedLowB=5e-2;		%low constant speed of bacteria B
%speedA=[speedLowA speedHighA];
%speedB=[speedLowB speedHighB];

%fixed speed
speedA=1e-2;
speedB=1e-2;

%variable turning frequency
lambda0HighA=10e-3;	%high base turning frequency of bacteria A
lambda0LowA=10e-3;		%low base turning frequency of bacteria A
lambda0HighB=1.5e-3;	%high base turning frequency of bacteria B
lambda0LowB=5e-3;		%low base turning frequency of bacteria B
lambda0A=[lambda0LowA lambda0HighA];
lambda0B=[lambda0LowB lambda0HighB];

%fixed turning frequency
%lambda0A=1.5e-3;
%lambda0B=1.5e-3;

kappaA=2;		%chemotactic sensitivity constant bacteria A
kappaB=2;		%chemotactic sensitivity constant bacteria B
VthA=1.5;		%threshold concentration of AHL for bacteria A
VthB=1.0;		%threshold concentration of AHL for bacteria B

%initialize bacteria A
bacteriaA=[];

for i=1:nBacteriaA
	%normal distribution
	%x=normrnd(L/2,1);

	%uniform random distribution
	%x=rand*L-L/2;

	%uniform distribution
	%x=L/(nBacteriaA)*(i-1);

	%single peak
	%x=L*0.01;
	x=L/2;

	%two peaks
	%if i<nBacteriaA/2
	%	%x=0.25*L;
	%	x=normrnd(0.25*L,1);
	%else
	%	%x=0.75*L;
	%	x=normrnd(0.75*L,1);
	%end

	b=bacterium(x);
	bacteriaA=[bacteriaA b];
end

bacteriaPopA=bacteriaPopulationA(bacteriaA,domain);
%bacteriaPopA=bacteriaPopulationStill(bacteriaA,domain);

%initialize bacteria B
bacteriaB=[];

for i=1:nBacteriaB
	%normal distribution
	%x=normrnd(L/2,1);

	%uniform random distribution
	%x=rand*L-L/2;

	%uniform distribution
	%x=L/(nBacteriaB)*(i-1);

	%single peak
	%x=5;
	%x=10;
	%x=L*(.99);
	x=L/2;

	%two peaks
	%if i<nBacteriaB/2
	%	%x=.35*L;
	%	x=normrnd(0.35*L,1);
	%else
	%	%x=.65*L;
	%	x=normrnd(0.65*L,1);
	%end

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
%bandwidth=1;

%define timestep
timestep=1;

%define scaling for plotting concentrations
scaling=30;

%initialize model
model=realModel_II(bacteriaPopA,bacteriaPopB,AHLField,leucineField,...
		alpha,beta,k1,k2,DAHL,Dleucine,lambda0A,lambda0B,speedA,speedB,kappaA,kappaB,VthA,VthB,...
		kernelfun,bandwidth,timestep,scaling);

%%
runSimulation;
makeVideo;
%%

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
