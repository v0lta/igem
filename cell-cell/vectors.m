close all;clear all;clc;

%Interaction parameters
r0=2;
k=1;

%Initial positions
%Cell 1
x1=2.5;
y1=2.5;

%Cell 2
x2=3.5;
y2=2.5;

%Calculate r and absolute force
r=sqrt((x2-x1)^2+(y2-y1)^2);
Fr=k*(r0-r);

%Calculate force vectors
F1=[x1-x2;y1-y2]*1/r*Fr;
F2=[x2-x1;y2-y1]*1/r*Fr;

%plot results

fig=figure(1);
hold on;

plot(x1,y1,'b*');
plot(x2,y2,'c*');

quiver(x1,y1,F1(1),F1(2));
quiver(x2,y2,F2(1),F2(2));

xlim([0 5]);
ylim([0 5]);
