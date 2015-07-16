%Model I
close all;clear all;clc;
filename='simulation_periodic.avi';
framerate=10;

N=200;			%time steps
nBacteria=1;	%number of bacteria
%nBacteria=1;	%number of bacteria
%k=20;			%plot kth iteration
L=1;			%Length of domain
dx=0.1;			%grid spacing
domain=[0:dx:L];%domain
%domain=[0:5];%domain

%define constants
muHigh=1/30;	%high diffusion constant of bacteria
muLow=1/300;	%low diffusion constant of bacteria
mu=[muLow muHigh];
Vth=1.2;		%threshold concentration of attractant
kappa=2;		%chemotactic sensitivity constant
eta=1e-3;		%bacterial production rate of attractant
%eta=0;			%bacterial production rate of attractant
Da=1/300;		%Diffusion constant of attractant
alpha=5e-3;		%degradation rate of attractant
Ds=1/300;		%Diffusion constant of nutrient
beta=1e-3;		%bacterial consumption rate of nutrient

%initialize bacteria
bacteria=[];

for i=1:nBacteria
	%normal distribution
	%x=normrnd(0,1);

	%uniform random distribution
	%x=rand*L-L/2;

	%uniform distribution
	%x=L/(nBacteria+1)*i;

	%single peak
	%x=5;
	x=1;

	b=bacterium(x);
	bacteria=[bacteria b];
end

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

bacteriaPop=bacteriaPopulation(bacteria,domain);

%initialize attractant field
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
attractantField=attractant(domain,concentration,boundaries);

%initialize nutrient field
n=length(domain);
concentration=zeros(1,n)+1;
%concentration=cos(domain)+1;
%concentration=sin(domain)+1;

%boundary conditions
boundaries=[1 1];
nutrientField=nutrient(domain,concentration,boundaries);

%define kernel function and bandwidth
addpath ..\..\..\kernel;
%kernelfun=@uni;
%kernelfun=@tri;
%kernelfun=@(x) normpdf(x,0,1);
kernelfun=@epanechnikov;
bandwidth=0.25;

%define timestep
timestep=1;

%initialize model
model=chemotaxisModel(bacteriaPop,attractantField,nutrientField,mu,Vth,kappa,eta,Da,alpha,...
	Ds,beta,kernelfun,bandwidth,timestep);

%% run for (extra) N time steps
%N=800;
for i=1:N
	model.update()
end

%plot kth iteration
%fig=figure(1);
%model.plot(k,fig);
%% 
%vidObj=VideoWriter('simulation2.avi');
vidObj=VideoWriter(filename);
set(vidObj,'FrameRate',framerate);
open(vidObj);

nFrames=model.getlength();
fig=figure(1);
for i=1:nFrames
	model.plot(i,fig);
	ylim([-5 75]);
	writeVideo(vidObj,getframe());
	clf;
end

close(fig);
vidObj.close();

beep on;
beep;
beep off;
%% 
