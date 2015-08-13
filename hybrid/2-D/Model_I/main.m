%Model I
close all;clear all;clc;
filename=['zero_flux_',...
	'zero_initial_concentration_',...
	'high_degradation_',...
	'spot_A_',...
	'spot_B_',...
	'small_simulation'];
framerate=10;

%Simulation parameters
N=100;				%time steps
%nBacteriaA=1000;	%number of bacteria A
nBacteriaA=300;		%number of bacteria A
%nBacteriaA=1;		%number of bacteria A
%nBacteriaB=1000;	%number of bacteria B
nBacteriaB=300;		%number of bacteria B
%nBacteriaB=1;		%number of bacteria B
XLength=25;			%Length of domain
YLength=25;			%Length of domain
Jx=101;				%# of subdivisions
Jy=101;				%# of subdivisions

%define domain
domain.x=linspace(0,XLength,Jx);
domain.y=linspace(0,YLength,Jy);
%domain=[domainx',domainy'];

%define kernel functions and bandwidth
kernelfun=@epanechnikov2DNorm;
bandwidth=.8;

%define timestep
dt=1;

%define constants
alpha=3e-3;		%production rate of AHL
%alpha=-3e-3;	%production rate of AHL
beta=4e-3;		%production rate of leucine
%k1=5e-3;		%degradation rate of AHL
%k1=0;			%degradation rate of AHL
k1=5e-1;		%degradation rate of AHL
k2=3e-1;		%degradation rate of leucine
%DAHL=1/300;	%Diffusion constant of AHL
DAHL=1/30;		%Diffusion constant of AHL
Dleucine=1/20;	%Diffusion constant of leucine
%muA=1/30;
muHighA=1/30;	%high diffusion constant of bacteria A
muLowA=1/300;	%low diffusion constant of bacteria A
muHighB=1/30;	%high diffusion constant of bacteria B
muLowB=1/300;	%low diffusion constant of bacteria B
muA.low=muLowA;
muA.high=muHighA;
muB.low=muLowB;
muB.high=muHighB;
%alpha=0;		%production rate of AHL
kappaA=2;		%chemotactic sensitivity constant of bacteria A
kappaB=2;		%chemotactic sensitivity constant of bacteria B
VthA=1.2;		%threshold concentration of AHL for bacteria A
VthB=1.0;		%threshold concentration of AHL for bacteria B

%Initialize bacteria A
bacteriaA=[];

for i=1:nBacteriaA
	%x=XLength/2;
	%y=YLength/2;

	%gaussian
	x=normrnd(1/2*XLength,1);
	y=normrnd(1/2*YLength,1);
	%x=normrnd(9/10*XLength,1);
	%y=normrnd(9/10*YLength,1);

	%uniform random
	%x=rand*XLength;
	%y=rand*YLength;

	if x>XLength
		x=XLength;
	end

	if y>YLength
		y=YLength;
	end

	b=bacterium(x,y);
	bacteriaA=[bacteriaA b];
end

bacteriaPopA=bacteriaPopulationA(bacteriaA,domain);

%Initialize bacteria B
bacteriaB=[];

for i=1:nBacteriaB
	%x=XLength/2;
	%y=YLength/2;

	%gaussian
	x=normrnd(XLength/2,1);
	y=normrnd(YLength/2,1);

	%uniform random
	%x=rand*XLength;
	%y=rand*YLength;

	if x>XLength
		x=XLength;
	end

	if y>YLength
		y=YLength;
	end

	b=bacterium(x,y);
	bacteriaB=[bacteriaB b];
end

bacteriaPopB=bacteriaPopulationB(bacteriaB,domain);

%initialize AHL field
m=length(domain.y);
n=length(domain.x);
[X,Y]=meshgrid(domain.x,domain.y);

%uniform
%concentration=zeros(m,n)+1;
concentration=zeros(m,n)+1e-5;

%block
%a=floor(J/3);
%idy=a:J-a;
%idx=a:J-a;
%concentration(idy,idx)=ones(length(idy),length(idx));

%boundary conditions
west=concentration(:,1);
east=concentration(:,end);
south=concentration(1,:);
north=concentration(end,:);
boundaries=[west,east,south',north'];

AHLField=AHL(domain,concentration,boundaries);

%initialize leucine field
m=length(domain.y);
n=length(domain.x);
[X,Y]=meshgrid(domain.x,domain.y);

%uniform
%concentration=zeros(m,n)+1;
concentration=zeros(m,n)+1e-5;

%block
%a=floor(J/3);
%idy=a:J-a;
%idx=a:J-a;
%concentration(idy,idx)=ones(length(idy),length(idx));

%boundary conditions
west=concentration(:,1);
east=concentration(:,end);
south=concentration(1,:);
north=concentration(end,:);
boundaries=[west,east,south',north'];

leucineField=leucine(domain,concentration,boundaries);

%define scaling for plotting concentrations
scaling=20;

model=model1(bacteriaPopA,bacteriaPopB,AHLField,leucineField,...
	alpha,beta,k1,k2,muA,muB,DAHL,Dleucine,kappaA,kappaB,VthA,VthB,...
	kernelfun,bandwidth,dt,scaling);

%% run simulation
disp('Running simulation');
for i=1:N
	model.update();
end

beep on;
beep;
beep off;

disp('Simulation finished');

%% preview
disp('Preview of simulation');
preview;

%% save workspace and make videos
%makevideo;

disp('End!');
