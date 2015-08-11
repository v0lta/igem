%Model I
close all;clear all;clc;
filename='diffusion_consumption_of_AHL';
framerate=10;

%Simulation parameters
N=100;			%time steps
nBacteriaA=300;	%number of bacteria A
%nBacteriaA=1;	%number of bacteria A
XLength=30;		%Length of domain
YLength=30;		%Length of domain
J=100;			%# of subdivisions

%define domain
domain.x=linspace(0,XLength,J);
domain.y=linspace(0,YLength,J);
%domain=[domainx',domainy'];

%define kernel functions and bandwidth
kernelfun=@epanechnikov2DNorm;
bandwidth=.5;

%define timestep
dt=1;

%define constants
%muA=1/30;
muA=1/50;		%diffusion constant of bacteria A
alpha=-3e-3;	%production rate of AHL
%alpha=0;		%production rate of AHL
kappaA=2;		%chemotactic sensitivity constant
DAHL=1/300;		%Diffusion constant of AHL

%Initialize bacteria A
bacteriaA=[];

for i=1:nBacteriaA
	%x=XLength/2;
	%y=YLength/2;

	x=normrnd(XLength/2,1);
	y=normrnd(YLength/2,1);

	b=bacterium(x,y);
	bacteriaA=[bacteriaA b];
end

bacteriaPopA=bacteriaPopulationA(bacteriaA,domain);

%initialize AHL field
%uniform
m=length(domain.y);
n=length(domain.x);
[X,Y]=meshgrid(domain.x,domain.y);
concentration=zeros(m,n)+1;
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

%define scaling for plotting concentrations
scaling=20;

model=model1(bacteriaPopA,AHLField,alpha,muA,DAHL,kappaA,kernelfun,bandwidth,dt,scaling);

%% run simulation
for i=1:N
	model.update();
end

%% save workspace
save(filename);

%% Make videos
vidObj3D=VideoWriter([filename '_3D.avi']);
set(vidObj3D,'FrameRate',framerate);
open(vidObj3D);

vidObj2D=VideoWriter([filename '_2D.avi']);
set(vidObj2D,'FrameRate',framerate);
open(vidObj2D);

nFrames=model.getlength();
fig=figure();
set(fig,'units','normalized','outerposition',[0 0 1 1]);
%fig=figure(1);
for i=1:nFrames
	%3D
	model.plot3D(i,fig);
	%ylim([-5 75]);
	writeVideo(vidObj3D,getframe(fig));
	clf;
end

for i=1:nFrames
	%2D
	model.plot2D(i,fig);
	%ylim([-5 75]);
	writeVideo(vidObj2D,getframe(fig));
	clf;
end

close(fig);
vidObj3D.close();
vidObj2D.close();

beep on;
beep;
beep off;
