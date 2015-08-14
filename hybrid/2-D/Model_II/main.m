%Model I
close all;clear all;clc;
filename=['zero_flux_',...
	'zero_initial_concentration_',...
	'high_degradation_',...
	'spot_A_',...
	'spot_B_',...
	'square_domain_',...
	'small_simulation'];
framerate=10;

%Simulation parameters
N=600;				%time steps
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

%variable speed
%speedHighA=1e-2;	%high constant speed of bacteria A
%speedLowA=.5e-2;		%low constant speed of bacteria A
%speedHighB=1e-2;	%high constant speed of bacteria B
%speedLowB=0.7e-2;		%low constant speed of bacteria B
%%speedLowB=5e-2;		%low constant speed of bacteria B
%speedA=[speedLowA speedHighA];
%speedB=[speedLowB speedHighB];

%fixed speed
speedA=1.2e-2;
speedB=1.2e-2;

%variable turning frequency
lambda0HighA=1.5e-2;	%high base turning frequency of bacteria A
lambda0LowA=1.5e-3;		%low base turning frequency of bacteria A
lambda0HighB=1.5e-2;	%high base turning frequency of bacteria B
lambda0LowB=1.5e-3;		%low base turning frequency of bacteria B
lambda0A.low=lambda0LowA;
lambda0A.high=lambda0HighA;
lambda0B.low=lambda0LowB;
lambda0B.high=lambda0HighB;

%fixed turning frequency
%lambda0A=1.5e-3;
%lambda0B=1.5e-3;

%alpha=0;		%production rate of AHL
kappaA=2;		%chemotactic sensitivity constant of bacteria A
kappaB=2;		%chemotactic sensitivity constant of bacteria B
%kappaA=0;		%chemotactic sensitivity constant of bacteria A
%kappaB=0;		%chemotactic sensitivity constant of bacteria B
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
	elseif x<0
		x=0
	end

	if y>YLength
		y=YLength;
	elseif y<0
		y=0;
	end

	direction=rand*2*pi;
	b=bacterium(x,y,direction);
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
	elseif x<0
		x=0
	end

	if y>YLength
		y=YLength;
	elseif y<0
		y=0
	end

	direction=rand*2*pi;
	b=bacterium(x,y,direction);
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
if Jx<Jy
	foo=0;
	bar=Jy-Jx;
else
	foo=Jx-Jy;
	bar=0;
end
boundaries=[[west;zeros(foo,1)],[east;zeros(foo,1)],[south';zeros(bar,1)],[north';zeros(bar,1)]];

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
if Jx<Jy
	foo=0;
	bar=Jy-Jx;
else
	foo=Jx-Jy;
	bar=0;
end
boundaries=[[west;zeros(foo,1)],[east;zeros(foo,1)],[south';zeros(bar,1)],[north';zeros(bar,1)]];

leucineField=leucine(domain,concentration,boundaries);

%define scaling for plotting concentrations
scaling=20;

model=model2(bacteriaPopA,bacteriaPopB,AHLField,leucineField,...
	alpha,beta,k1,k2,lambda0A,lambda0B,speedA,speedB,DAHL,Dleucine,kappaA,kappaB,VthA,VthB,...
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
%disp('Preview of simulation');
%preview;

%% save workspace and make videos
makevideo;

disp('End!');
