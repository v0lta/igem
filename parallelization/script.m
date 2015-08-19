close all;clear all;clc;
n=16;
alpha=20;

a=0;
array=(1:n)+alpha;

%% parallel for loop to compute sum of squares of numbers 1 to 5
parfor i=1:n
%for i=1:n
	a=a+fib(array(i));
    %A(i)=fib(array(i));
end

%a=sum(A)
a

beep on; beep; beep off;