clear all;close all;clc;
%define stripe pattern
p=pattern;
stripes=[];

%stripes=[stripes stripe(0,100)];
%border=circleBorder(150);

stripes=[stripes stripe(0,0.8e4)];
stripes=[stripes stripe(1e4,1.3e4)];
stripes=[stripes stripe(1.5e4,1.8e4)];
stripes=[stripes stripe(2e4,2.3e4)];
stripes=[stripes stripe(2.5e4,2.8e4)];
border=circleBorder(3e4);


%stripes=[stripes stripe(0,0.8e3)];
%stripes=[stripes stripe(1e3,1.3e3)];
%stripes=[stripes stripe(1.5e3,1.8e3)];
%stripes=[stripes stripe(2e3,2.3e3)];
%stripes=[stripes stripe(2.5e3,2.8e3)];
%border=circleBorder(1.5e3);
%border=circleBorder(3e3);

%Add patterns and border
for stripe=stripes
	p.addComp(stripe);
end
p.addBorder(border);

fig=figure(1);
p.plot(fig);
axis equal;

%timesteps
for N=1e8
	disp([num2str(N) ' timesteps']);
	%disp('test');
	tic;

	b=bacterium(0,0);

	%path=[0 0];

	for i=1:N
		b.updateCoordinates(p);
	end

	disp('Ratio of high density vs total:');
	disp(num2str(b.highCounter/N));

	b.plotPath(fig);
	disp(['N: ' num2str(N)]);
	toc;
end

