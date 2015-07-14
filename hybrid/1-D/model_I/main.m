%Model I
close all;clear all;clc;

N=200;			%time steps
nBacteria=100;	%number of bacteria
%nBacteria=1;	%number of bacteria
%k=20;			%plot kth iteration
domain=[-10:0.1:10];%domain
%domain=[0:5];%domain


%initialize bacteria
bacteria=[];
for i=1:nBacteria
	%x=normrnd(0,1);
	%b=bacterium(5);
	b=bacterium(0);
	%b=bacterium(x);
	bacteria=[bacteria b];
end
bacteriaPop=bacteriaPopulation(bacteria,domain);

%initialize nutrient field
n=length(domain);
concentration=zeros(1,n)+1;
%concentration=cos(domain)+1;
%concentration=sin(domain)+1;

%boundary conditions
boundaries=[1 1];
nutrientField=nutrient(domain,concentration,boundaries);

%define constants
mu=1/30;	%diffusion constant of bacteria
kappa=2;	%chemotactic sensitivity constant
d=1e-3;		%consumption rate
Ds=1/300;	%Diffusion constant of nutrient

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
model=chemotaxisModel(bacteriaPop,nutrientField,mu,kappa,d,Ds,kernelfun,bandwidth,timestep);

%run for N time steps
for i=1:N
	model.update()
end

%plot kth iteration
%fig=figure(1);
%model.plot(k,fig);

%vidObj=VideoWriter('simulation2.avi');
vidObj=VideoWriter('simulation_sin.avi');
set(vidObj,'FrameRate',10);
open(vidObj);

fig=figure(1);
for i=1:N+1
	model.plot(i,fig);
	ylim([0 75]);
	writeVideo(vidObj,getframe());
	clf;
end

close(fig);
vidObj.close();

%densityfun=bacteriaPop.bacteriadensity(d,bandwidth);
%
%y=densityfun(x);
%
%fig1=figure(1);
%plot(x,y);
%hold on;
%x=bacteriaPop.coordinates();
%y=x*0;
%plot(x,y,'*');
