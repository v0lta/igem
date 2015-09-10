clear all;close all;
%N=10;
%
%%save('test.mat','-v7.3');
%save('test.mat');
%m=matfile('test.mat','Writable',true,'-v7.3');
%m.A=zeros(N,N,2);
%
%B=magic(N);
%
%disp('piep');
%for i=1:N
%	%m.A=cat(3,m.A,magic(N));
%	m.A(:,:,i)=B;
%end
%
%parfor i=1:N
%	disp(['i: ' num2str(i)]);
%	m.A(:,:,i);
%end

XLength=10;
rsearch=1;
x=4.5;
y=1;

cellid=determinecellidtest(x,y,XLength,rsearch);
disp(x);
disp(y);
disp(cellid);
