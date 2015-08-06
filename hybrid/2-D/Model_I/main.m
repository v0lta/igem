%Model I
close all;clear all;clc;
filename='test.avi';
framerate=10;

%Simulation parameters
N=100;			%time steps
nBacteriaA=100;	%number of bacteria A
XLength=10;		%Length of domain
YLength=10;		%Length of domain
J=100;			%# of subdivisions

%define domain
domainx=linspace(0,XLength,J);
domainy=linspace(0,YLength,J);
domain=[domainx',domainy'];

%define kernel functions and bandwidth
kernelfun=@epanechnikov2DNorm;
bandwidth=.5;

%define timestep
dt=1;

%define constants
%muA=1/30;
muA=1/300;

%Initialize bacteria A
bacteriaA=[];

for i=1:nBacteriaA
	x=XLength/2;
	y=YLength/2;

	b=bacterium(x,y);
	bacteriaA=[bacteriaA b];
end

bacteriaPopA=bacteriaPopulationA(bacteriaA,domain);

model=model1(bacteriaPopA,muA,kernelfun,bandwidth,dt);

%% run simulation
for i=1:N
	model.update();
end

%% Make video
vidObj=VideoWriter(filename);
set(vidObj,'FrameRate',framerate);
open(vidObj);

nFrames=model.getlength();
fig=figure('Position',[100 100 1000 1000]);
set(fig,'units','normalized','outerposition',[0 0 1 1]);
%fig=figure(1);
for i=1:nFrames
	model.plot(i,fig);
	%ylim([-5 75]);
	writeVideo(vidObj,getframe());
	clf;
end

close(fig);
vidObj.close();

beep on;
beep;
beep off;
