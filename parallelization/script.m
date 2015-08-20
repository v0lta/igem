close all;clear all;clc;

k=6144;
%A=distributed(rand(k));
A=gpuArray.rand(k);
%A=rand(k);

%B=distributed(rand(k));
B=gpuArray.rand(k);
%B=rand(k);

%serial
Fserial=A\B;
%Fserial=A*B;
%Fserial=gather(Fserial);

%parallel
%[m,n]=size(B);
%
%Fparallel=[];
%parfor i=1:n
%	Fparallel=[Fparallel,A*B(:,i)];
%end

%if sum(sum(Fserial==Fparallel)) ~= m*n
%	disp('Matrices are not equal');
%else
%	disp('Matrices are equal');
%end

beep on;
beep;
beep off;

%mldivide
%k=6000
%parallel: 59.413s
%serial: 41.405s
%gpu: 29.968s

%serial: 19.209s
%gpu: 26.995s
%A gpu, B serial: 27.084s

%serial: 23.866s
%gpu: 26.330s

%k=6144
%serial
%gpu

%k=1000
%parallel: 5.211s
%serial: 0.545s
%gpu: 0.408s
%A gpu, B serial: 0.267s

%matrix multiplication
%k=6000
%serial: 13.555s
%gpu: 5.331s
