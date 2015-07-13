clear all;close all;clc;

%Parameters
n=20;			%number of cells
h=0.9;			%bandwidth
x=-4:0.1:4;		%domain

%sample from normal probability
pointArray=normrnd(zeros(1,n),1);
%pointArray=rand(1,n);

%define kernel function
%d=@uni;
d=@tri;
%d=@(x) normpdf(x,0,1);
%d=@epanechnikov;

%calculate kernel function
f=KDE(pointArray,d,h);
y=f(x);

%plot kernel
fig1=figure(1);
%y=[];
%for ix=x
%	y(end+1)=K_h(ix);
%end
plot(x,y,'-b');
hold on;
plot(pointArray,zeros(1,n),'k*');
