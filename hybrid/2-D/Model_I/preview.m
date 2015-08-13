close all;

nFrames=model.getlength();
fig=figure();
set(fig,'units','normalized','outerposition',[0 0 1 1]);
%fig=figure(1);
for i=1:nFrames
	%3D
	model.plot3D(i,fig);
	%ylim([-5 75]);
	pause(0.1);
	clf;
end

for i=1:nFrames
	%2D
	model.plot2D(i,fig);
	pause(0.1);
	%ylim([-5 75]);
	clf;
end

