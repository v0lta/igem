function output=fib(n)
if n==1
	output=0;
elseif n==2
	output=1;
else
	output=fib(n-1)+fib(n-2);
end
