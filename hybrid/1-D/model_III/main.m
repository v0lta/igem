%Model II
close all;clear all;clc;

N=200;			%time steps
nBacteria=100;	%number of bacteria
%nBacteria=1;	%number of bacteria
%k=20;			%plot kth iteration
%domain=-10:0.01:10;%domain
%domain=-5:0.1:5;%domain
domain=-5:0.01:5;%domain
%timestep=1e-5;	%time stepsize
timestep=1;


%initialize bacteria
bacteria=[];
for i=1:nBacteria
	%x=normrnd(0,1);
	%b=bacterium(5);
	b=bacterium(0);
	bacteria=[bacteria b];
end
bacteriaPop=bacteriaPopulation(bacteria,domain);

%initialize nutrient field
n=length(domain);
concentration=zeros(1,n)+1;
%concentration=zeros(1,n);
%concentration(8:12)=1;
%concentration=domain;
nutrientField=nutrient(domain,concentration);

%define constants
%mu=1/30;		%diffusion constant of bacteria
kappa=2;		%chemotactic sensitivity constant
d=1e-3;			%consumption rate
%d=0;			%consumption rate
lambda0=1.5e-3;	%base turning frequency
speed=1e-2;		%constant speed
%Ds=1;		%diffusion constant of nutrient
Ds=0;		%diffusion constant of nutrient

%define kernel function and bandwidth
addpath ..\..\..\kernel;
%kernelfun=@uni;
%kernelfun=@tri;
%kernelfun=@(x) normpdf(x,0,1);
kernelfun=@epanechnikov;
bandwidth=0.25;

%initialize model
model=chemotaxisModel(bacteriaPop,nutrientField,kappa,d,lambda0,speed,Ds,kernelfun,bandwidth,timestep);

%run for N time steps
for i=1:N
	model.update()
end

%plot kth iteration
%fig=figure(1);
%model.plot(k,fig);

vidObj=VideoWriter('test5.avi');
set(vidObj,'FrameRate',30);
open(vidObj);

fig=figure(1);
for i=1:N+1
	model.plot(i,fig);
	ylim([0 75]);
	%ylim([0 1.2]);
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
