clear all;close all;clc;

%% Calculation of R_cutoff for given mu, dt, r0 and p
%mu=1/300;
%dt=3.3248e-4;
%r0=1/100;
%
%lambda=(2*r0)^2/(4*mu*dt);
%
%p=0.99;	% 5 percent chance of detaching
%v=2;	% 2 degrees of freedom
%
%R_star2=ncx2inv(p,v,lambda);
%
%R_cutoff=sqrt(4*mu*dt*R_star2);
%
%disp(['R_cutoff: ' num2str(R_cutoff)]);

%for
%mu=1/300;
%dt=1;
%r0=2;
%p=0.95; 

%R cutoff is 2.1931

%% Calculation of k for given mu, dt, gamma, r0, rcut
%mu=1/300;
%dt=1;
%r0=2;
%gamma=2;
%
%rcut=2.25;
%
%dx=(r0-rcut)/2;
%
%k=dx/(mu*dt*gamma*(r0-rcut));
%
%disp(['k: ' num2str(k)]);

%% Determination of suitable dt for given r
%mu=1/30;
mu=1/300;
%dt=1;
%r0=1/100;
r0=0.05;

%lambda2=(2*r0)^2/(4*mu*dt)

p=0.99;	% 5 percent chance of detaching
v=2;	% 2 degrees of freedom

f=@(dt) ncx2inv(p,v,(2*r0)^2/(4*mu*dt))*(4*mu*dt)-(2*r0*1.25)^2;

dt=fsolve(f,0.1)

%dt = 3.3248e-4
