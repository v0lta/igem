close all;clear all;clc;

N=100;			%time steps
nBacteria=100;	%number of bacteria
%k=20;			%plot kth iteration
domain=-10:0.01:10;%domain


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
concentration=domain*0+1;
nutrientField=nutrient(domain,concentration);

%define constants
mu=1/30;	%diffusion constant
kappa=2;	%chemotactic sensitivity constant
d=1e-3;		%consumption rate

%define kernel function and bandwidth
addpath ..\..\..\kernel;
%kernelfun=@uni;
%kernelfun=@tri;
%kernelfun=@(x) normpdf(x,0,1);
kernelfun=@epanechnikov;
bandwidth=0.25;

%initialize model
model=chemotaxisModel(bacteriaPop,nutrientField,mu,kappa,d,kernelfun,bandwidth);

%run for N time steps
for i=1:N
	model.update()
end

%plot kth iteration
%fig=figure(1);
%model.plot(k,fig);

vidObj=VideoWriter('test.avi');
set(vidObj,'FrameRate',5);
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
